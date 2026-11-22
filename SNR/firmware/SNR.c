/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    main.c
  * @brief   SNR measurement firmware
  *
  * Project : Implementation and Performance Analysis of a Baseband Phase
  *           Modulation System with STM32 for Low-Power Wireless Communication
  * Board   : STM32 Nucleo-F303RE
  * Output  : DAC1_OUT1 on PA4 -> Sallen-Key Butterworth LPF -> oscilloscope
  *
  * ---------------------------------------------------------------------------
  * HOW TO USE
  * ---------------------------------------------------------------------------
  * Change ONE line to select the build:
  *
  *     #define SNR_MODE   MODE_CARRIER
  *
  * Then Project -> Clean, then Project -> Build.
  *
  *   Mode              LED flashes   Measurement case
  *   ----------------- ------------- --------------------------------------
  *   MODE_DC_DMA            1        N1 - instrument noise floor
  *   MODE_CARRIER           2        N2 / N3 - unmodulated 5 kHz carrier
  *   MODE_PM_1              3        N4 - analog PM, kp = 1.0
  *   MODE_PM_BESSEL         4        N5 - analog PM, kp = 2.4048
  *   MODE_BPSK              5        N6 - BPSK, 500 bps
  *   MODE_DC_STATIC         6        optional - DMA stopped, DAC at 2048
  *
  * ---------------------------------------------------------------------------
  * NOTE ON PROJECT STRUCTURE
  * ---------------------------------------------------------------------------
  * This file assumes CubeMX was set to generate peripheral initialisation as
  * separate .c/.h file pairs, so MX_GPIO_Init, MX_DMA_Init, MX_DAC1_Init and
  * MX_TIM6_Init live in gpio.c, dma.c, dac.c and tim.c.
  *
  * If instead your project keeps all the MX_..._Init functions inside main.c,
  * delete the four #include lines for "dac.h", "dma.h", "gpio.h" and "tim.h",
  * and paste your existing MX_..._Init function bodies back in at the marked
  * place near the bottom of this file.
  *
  * ---------------------------------------------------------------------------
  * NOTE ON THE CLOCK
  * ---------------------------------------------------------------------------
  * SystemClock_Config below uses HSI directly at 8 MHz with NO PLL, because
  * the HSI PLL path on the F303 is ambiguous between HSI/2 and HSI/PREDIV.
  *
  * If your existing SystemClock_Config already works and was the one used for
  * the output frequency accuracy measurement, KEEP YOURS and delete the one
  * here. Method B of the SNR analysis cross-checks the fitted carrier
  * frequency against that measurement, and the check is only meaningful if
  * both used the same clock.
  *
  * Do not let CubeMX regenerate this function.
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "dac.h"
#include "dma.h"
#include "gpio.h"
#include "tim.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <math.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* ---------- mode numbers ---------- */
#define MODE_DC_DMA      0   /* N1 floor  - DMA runs, buffer is all midcode */
#define MODE_CARRIER     1   /* N2 / N3   - unmodulated 5 kHz carrier       */
#define MODE_PM_1        2   /* N4        - analog PM, kp = 1.0             */
#define MODE_PM_BESSEL   3   /* N5        - analog PM, kp = 2.4048          */
#define MODE_BPSK        4   /* N6        - BPSK, 500 bps                   */
#define MODE_DC_STATIC   5   /* optional  - DMA stopped, DAC held at 2048   */

/* =========================================================== */
/* ===== THE ONLY LINE YOU CHANGE BETWEEN BUILDS ============== */
#define SNR_MODE   MODE_CARRIER
/* =========================================================== */

/* ---------- fixed system parameters ---------- */
/* These match every other methodology document in the thesis. */
#define FS_HZ        100000UL   /* DAC sample rate                */
#define FC_HZ          5000UL   /* carrier frequency              */
#define FM_HZ           500UL   /* modulating frequency           */
#define RS_BPS          500UL   /* BPSK symbol rate               */
#define N_SAMPLES      1000UL   /* DMA buffer length              */
#define DAC_MID         2048    /* 12-bit midpoint                */
#define DAC_AMP         1200    /* amplitude in codes             */

/* ---------- phase deviation for the PM modes ---------- */
#if   (SNR_MODE == MODE_PM_1)
  #define KP   1.0f
#elif (SNR_MODE == MODE_PM_BESSEL)
  #define KP   2.4048f          /* first zero of J0 - carrier null */
#else
  #define KP   0.0f
#endif

/* ---------- BPSK bit pattern ---------- */
/* One buffer is 10 ms long. At 500 bps that is exactly 5 symbols of
   200 samples each. The pattern below repeats every buffer and gives
   three phase reversals: 0, 1, 1, 0, 1.                             */
#define BPSK_PATTERN     0x16u  /* 0b10110 */
#define BPSK_NBITS       5u

#define TWO_PI   6.283185307179586f
#define PI_F     3.141592653589793f

/* ---------- status LED ---------- */
/* LD2 on the Nucleo-F303RE is PA5. These are normally defined in
   main.h by CubeMX; the fallback covers the case where they are not. */
#ifndef LD2_Pin
  #define LD2_Pin        GPIO_PIN_5
  #define LD2_GPIO_Port  GPIOA
#endif

/* ---------- safety: MCO must be off ---------- */
/* Left over from frequency accuracy test T4. If it is still enabled it
   puts a clock signal on an output pin and injects a tone into the
   supply rails, which would land straight in the noise floor.        */
#ifdef MCO_ENABLE
  #if (MCO_ENABLE != 0)
    #error "MCO_ENABLE must be 0 for the SNR measurement."
  #endif
#endif

/* ---------- compile-time checks on buffer coherence ---------- */
/* The buffer must hold a whole number of carrier cycles, modulating
   cycles and BPSK symbols. If not, the waveform jumps at every DMA
   wrap. The buffer repeats at 100 Hz, so any jump there produces a
   comb of spurious lines 100 Hz apart across the whole spectrum -
   counted as noise, and almost invisible on the scope display.
   The rule is N = k * lcm(fs/fc, fs/fm).                            */
#if (FS_HZ % FC_HZ) != 0
  #error "fs must be a whole multiple of fc"
#endif
#if (FS_HZ % FM_HZ) != 0
  #error "fs must be a whole multiple of fm"
#endif
#if (FS_HZ % RS_BPS) != 0
  #error "fs must be a whole multiple of the symbol rate"
#endif
#if (N_SAMPLES % (FS_HZ / FC_HZ)) != 0
  #error "N_SAMPLES is not a whole number of carrier cycles"
#endif
#if (N_SAMPLES % (FS_HZ / FM_HZ)) != 0
  #error "N_SAMPLES is not a whole number of modulating cycles"
#endif
#if (N_SAMPLES % (FS_HZ / RS_BPS)) != 0
  #error "N_SAMPLES is not a whole number of BPSK symbols"
#endif
#if (DAC_MID + DAC_AMP) > 4095
  #error "amplitude too large - the DAC would clip"
#endif
#if (DAC_MID - DAC_AMP) < 0
  #error "amplitude too large - the DAC would clip"
#endif

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */
/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
/* USER CODE BEGIN PV */

/* Defined in dac.c, dma.c and tim.c. The externs are harmless if the
   headers already declare them.                                      */
extern DAC_HandleTypeDef hdac1;
extern DMA_HandleTypeDef hdma_dac1_ch1;
extern TIM_HandleTypeDef htim6;

/* The DMA buffer. uint16_t because the DAC is used in 12-bit
   right-aligned mode (DHR12R1).                                      */
static uint16_t dac_buf[N_SAMPLES];

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */
static void led_init(void);
static void blink_mode(void);
static void build_waveform(void);
static void force_dac_config(void);
static void force_dma_config(void);
static void force_sample_rate(void);
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* ---------------------------------------------------------------------
 * led_init()
 *
 * Configures PA5 as a push-pull output. Harmless if MX_GPIO_Init has
 * already done it; included so the status LED works even if PA5 was
 * never set up in CubeMX.
 * ------------------------------------------------------------------- */
static void led_init(void)
{
  GPIO_InitTypeDef g = {0};

  __HAL_RCC_GPIOA_CLK_ENABLE();

  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  g.Pin   = LD2_Pin;
  g.Mode  = GPIO_MODE_OUTPUT_PP;
  g.Pull  = GPIO_NOPULL;
  g.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &g);
}

/* ---------------------------------------------------------------------
 * blink_mode()
 *
 * Flashes LD2 (SNR_MODE + 1) times, then leaves it off. Runs BEFORE the
 * DAC starts, so it cannot add noise to any capture.
 *
 * In the lab this tells you at a glance which build is on the board.
 * Without it you will eventually record 32 files of the wrong mode and
 * only find out during analysis, at home, after the lab has closed.
 * ------------------------------------------------------------------- */
static void blink_mode(void)
{
  for (uint32_t i = 0; i <= (uint32_t)SNR_MODE; i++)
  {
    HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_SET);
    HAL_Delay(150);
    HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);
    HAL_Delay(250);
  }
  HAL_Delay(500);
}

/* ---------------------------------------------------------------------
 * build_waveform()
 *
 * Fills dac_buf once at startup. Nothing is computed while the waveform
 * is running - the DMA just replays the buffer forever. This is the
 * same precomputation approach measured in the power chapter.
 * ------------------------------------------------------------------- */
static void build_waveform(void)
{
  for (uint32_t n = 0; n < N_SAMPLES; n++)
  {
#if (SNR_MODE == MODE_DC_DMA) || (SNR_MODE == MODE_DC_STATIC)

    /* Flat midcode.
       In MODE_DC_DMA every DMA transfer, every DAC write and every TIM6
       trigger still happens exactly as in a real measurement, but no
       signal is produced. That is what makes it a fair noise floor: the
       switching activity, the supply current pattern and the source
       impedance seen by the probe are all unchanged. Only the data is
       different.

       Simply disabling the DAC would give a quieter floor than reality
       and would make the SNR look better than it is.                  */
    dac_buf[n] = (uint16_t)DAC_MID;

#else

    /* Carrier phase. fc/fs = 5000/100000 = 1/20 exactly. */
    float carrier = TWO_PI * (float)FC_HZ * (float)n / (float)FS_HZ;

    /* Modulating phase term. */
    float mod = 0.0f;

  #if (SNR_MODE == MODE_PM_1) || (SNR_MODE == MODE_PM_BESSEL)
    mod = KP * sinf(TWO_PI * (float)FM_HZ * (float)n / (float)FS_HZ);

  #elif (SNR_MODE == MODE_BPSK)
    /* Symbol index. fs/Rs = 200 divides N exactly, so integer division
       is exact here and cannot drift.                                 */
    uint32_t sym = n / (FS_HZ / RS_BPS);
    uint32_t bit = (BPSK_PATTERN >> (sym % BPSK_NBITS)) & 1u;
    mod = bit ? PI_F : 0.0f;
  #endif

    float s = (float)DAC_MID + (float)DAC_AMP * sinf(carrier + mod);

    /* Round to nearest, then clamp. The clamp should never trigger -
       the compile-time checks guarantee it - but leaving it in means a
       mistyped amplitude shows as a flat top on the scope instead of
       wrapping around silently.                                       */
    int32_t code = (int32_t)(s + 0.5f);
    if (code < 0)    { code = 0;    }
    if (code > 4095) { code = 4095; }

    dac_buf[n] = (uint16_t)code;

#endif
  }
}

/* ---------------------------------------------------------------------
 * force_dac_config()
 *
 * Guarantees the two DAC settings the measurement depends on, whatever
 * CubeMX may have generated:
 *   - output buffer ENABLED, so the DAC can drive the LPF input
 *   - trigger source TIM6 TRGO
 * ------------------------------------------------------------------- */
static void force_dac_config(void)
{
  DAC_ChannelConfTypeDef sConfig = {0};

  sConfig.DAC_Trigger      = DAC_TRIGGER_T6_TRGO;
  sConfig.DAC_OutputBuffer = DAC_OUTPUTBUFFER_ENABLE;

  if (HAL_DAC_ConfigChannel(&hdac1, &sConfig, DAC_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }
}

/* ---------------------------------------------------------------------
 * force_dma_config()
 *
 * Guarantees circular mode and half-word transfers on both sides.
 *
 * Half-word matters: the buffer is uint16_t and the DAC is used in
 * 12-bit right-aligned mode. If either side is set to Word, every
 * second sample is dropped and the output frequency doubles.
 * ------------------------------------------------------------------- */
static void force_dma_config(void)
{
  HAL_DMA_DeInit(&hdma_dac1_ch1);

  hdma_dac1_ch1.Init.Direction           = DMA_MEMORY_TO_PERIPH;
  hdma_dac1_ch1.Init.PeriphInc           = DMA_PINC_DISABLE;
  hdma_dac1_ch1.Init.MemInc              = DMA_MINC_ENABLE;
  hdma_dac1_ch1.Init.PeriphDataAlignment = DMA_PDATAALIGN_HALFWORD;
  hdma_dac1_ch1.Init.MemDataAlignment    = DMA_MDATAALIGN_HALFWORD;
  hdma_dac1_ch1.Init.Mode                = DMA_CIRCULAR;
  hdma_dac1_ch1.Init.Priority            = DMA_PRIORITY_HIGH;

  if (HAL_DMA_Init(&hdma_dac1_ch1) != HAL_OK)
  {
    Error_Handler();
  }

  __HAL_LINKDMA(&hdac1, DMA_Handle1, hdma_dac1_ch1);
}

/* ---------------------------------------------------------------------
 * force_sample_rate()
 *
 * Sets TIM6 to overflow at exactly FS_HZ, working the ARR value out
 * from the clock that TIM6 is actually running on.
 *
 * This removes the most common firmware mistake on this board. TIM6
 * sits on APB1: if the APB1 prescaler is 1 the timer clock equals
 * PCLK1, but if it is anything else the timer clock is TWICE PCLK1.
 * Hard-coding ARR = 79 is only correct for one particular clock setup,
 * and getting it wrong puts the carrier at 2.5 kHz or 10 kHz.
 *
 * It also forces TRGO to Update Event, which is what triggers the DAC.
 * ------------------------------------------------------------------- */
static void force_sample_rate(void)
{
  RCC_ClkInitTypeDef clk = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};
  uint32_t flash_latency = 0;
  uint32_t tim_clk;
  uint32_t arr;

  HAL_RCC_GetClockConfig(&clk, &flash_latency);

  tim_clk = HAL_RCC_GetPCLK1Freq();
  if (clk.APB1CLKDivider != RCC_HCLK_DIV1)
  {
    tim_clk *= 2U;              /* APB1 prescaler > 1 doubles the timer clock */
  }

  /* The sample rate must come out exact. If the clock does not divide
     evenly by 100 kHz, stop rather than run at the wrong rate.        */
  if ((tim_clk % FS_HZ) != 0U)
  {
    Error_Handler();
  }

  arr = tim_clk / FS_HZ;        /* e.g. 8 MHz / 100 kHz = 80 */
  if ((arr < 2U) || (arr > 65536U))
  {
    Error_Handler();
  }

  htim6.Instance               = TIM6;
  htim6.Init.Prescaler         = 0U;
  htim6.Init.CounterMode       = TIM_COUNTERMODE_UP;
  htim6.Init.Period            = arr - 1U;   /* 8 MHz -> ARR = 79 */
  htim6.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;

  if (HAL_TIM_Base_Init(&htim6) != HAL_OK)
  {
    Error_Handler();
  }

  sMasterConfig.MasterOutputTrigger = TIM_TRGO_UPDATE;
  sMasterConfig.MasterSlaveMode     = TIM_MASTERSLAVEMODE_DISABLE;

  if (HAL_TIMEx_MasterConfigSynchronization(&htim6, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */
  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  HAL_Init();

  /* USER CODE BEGIN Init */
  /* USER CODE END Init */

  SystemClock_Config();

  /* USER CODE BEGIN SysInit */
  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_DAC1_Init();
  MX_TIM6_Init();

  /* USER CODE BEGIN 2 */

  /* Show which build is loaded, before anything else starts. */
  led_init();
  blink_mode();

  /* Force the settings the measurement depends on, so a stray CubeMX
     regeneration cannot silently change them.                        */
  force_dac_config();
  force_dma_config();
  force_sample_rate();

  /* Precompute the whole waveform. */
  build_waveform();

#if (SNR_MODE == MODE_DC_STATIC)

  /* DMA and timer stay stopped. The DAC output buffer is still on, so
     the op-amp sees the same source impedance as always.             */
  if (HAL_DAC_Start(&hdac1, DAC_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_1,
                       DAC_ALIGN_12B_R, DAC_MID) != HAL_OK)
  {
    Error_Handler();
  }

#else

  /* Normal path: circular DMA fed by TIM6 TRGO.
     The cast to uint32_t* is what the HAL prototype expects. The real
     transfer width comes from the DMA configuration set above.       */
  if (HAL_DAC_Start_DMA(&hdac1, DAC_CHANNEL_1,
                        (uint32_t *)dac_buf,
                        (uint32_t)N_SAMPLES,
                        DAC_ALIGN_12B_R) != HAL_OK)
  {
    Error_Handler();
  }

  if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
  {
    Error_Handler();
  }

#endif

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* Leave this loop empty.

       No printf, no LED toggling, no ADC reads. Anything the CPU does
       here injects switching noise into the supply rails, and on an SNR
       measurement that noise lands directly in the result.

       Do not use sleep mode either, even though the power chapter
       showed it saves current. Waking and sleeping creates current
       steps at a regular rate, and regular current steps appear as
       spurs in the spectrum. A plain idle loop keeps conditions
       identical across every capture, including the noise floor.     */

    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  *
  * HSI at 8 MHz, no PLL.
  *
  * The HSI PLL path on the STM32F303 is ambiguous between HSI/2 and
  * HSI/PREDIV, so the PLL is avoided entirely. At 8 MHz, FLASH_LATENCY_0
  * is correct and APB1 runs undivided, giving TIM6 an 8 MHz clock and
  * ARR = 79 for a 100 kSa/s sample rate.
  *
  * If you use HSE instead, remember the prediv lives inside
  * osc.PLL.PREDIV as RCC_PREDIV_DIV1, NOT in osc.HSEPredivValue.
  *
  * Do not let CubeMX regenerate this function.
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /* Configure HSI. No PLL. */
  RCC_OscInitStruct.OscillatorType      = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState            = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState        = RCC_PLL_NONE;

  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /* SYSCLK = HSI = 8 MHz, AHB /1, APB1 /1, APB2 /1 */
  RCC_ClkInitStruct.ClockType      = RCC_CLOCKTYPE_HCLK
                                   | RCC_CLOCKTYPE_SYSCLK
                                   | RCC_CLOCKTYPE_PCLK1
                                   | RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_HSI;
  RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_0) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */

/* If your project keeps the MX_..._Init functions inside main.c rather
   than in separate files, paste them here and remove the includes for
   "dac.h", "dma.h", "gpio.h" and "tim.h" at the top.                  */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */

  /* Fast continuous blink means something failed during setup - most
     likely the clock does not divide evenly to 100 kSa/s. Visible
     without a debugger, which matters in the lab.                    */
  __disable_irq();

  while (1)
  {
    HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
    for (volatile uint32_t d = 0; d < 100000U; d++)
    {
      __NOP();
    }
  }

  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  UNUSED(file);
  UNUSED(line);
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
