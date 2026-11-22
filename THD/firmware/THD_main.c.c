/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    main.c
  * @brief   Pure sine-wave generator for THD measurement
  *
  *          Thesis: Implementation and Performance Analysis of a Baseband
  *                  Phase Modulation System with STM32 for Low-Power
  *                  Wireless Communication
  *
  *          Target : STM32 Nucleo-F303RE
  *          Output : DAC1_OUT1 on PA4  -> Sallen-Key Butterworth LPF
  *          Engine : TIM6 TRGO -> DAC -> DMA1 Channel 3, circular mode
  *
  *          THIS BUILD PRODUCES A PURE SINE WAVE ONLY.
  *          No phase modulation. THD must never be computed from a PM
  *          waveform, because the intended sidebands cannot be separated
  *          from genuine harmonic distortion.
  *
  *          Change only the CONFIGURATION BLOCK below to move between the
  *          builds B1 to B7 in the THD test matrix.
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "dac.h"
#include "dma.h"
#include "tim.h"
#include "gpio.h"

/* USER CODE BEGIN Includes */
#include <math.h>
#include <stdint.h>
/* USER CODE END Includes */

/* USER CODE BEGIN PD */
/* ==========================================================================
 *  CONFIGURATION BLOCK  --  THIS IS THE ONLY PART YOU EDIT
 * ==========================================================================
 *
 *  Build table for the THD document
 *
 *   Build | OUTPUT_MODE      | F0_HZ | DAC_AMPLITUDE | Test cases
 *   ------|------------------|-------|---------------|-------------
 *    B1   | MODE_SINE        |  1000 |     1200      | T1, T2
 *    B2   | MODE_SINE        |  1000 |      500      | T5
 *    B3   | MODE_SINE        |  1000 |      800      | T6
 *    B4   | MODE_SINE        |  1000 |     1000      | T7
 *    B5   | MODE_SINE        |  1000 |     1500      | T8
 *    B6   | MODE_SINE        |  2000 |     1200      | T10 (optional)
 *    --   | MODE_SINE        |  5000 |     1200      | T3, T4 (reused
 *         |                  |       |               | from SNR N2/N3)
 *    B7a  | MODE_DC_MIDPOINT |   --  |      --       | T9 system floor
 *    B7b  | MODE_OFF         |   --  |      --       | T9 instrument floor
 *
 *  After editing, do a full rebuild (Project > Clean, then Build) so the
 *  lookup table is regenerated. Then flash.
 * ========================================================================== */

#define MODE_SINE           0   /* normal sine output                       */
#define MODE_DC_MIDPOINT    1   /* DAC on, constant code 2048, no DMA       */
#define MODE_OFF            2   /* DAC and DMA off, PA4 analog high-Z       */

#define OUTPUT_MODE         MODE_SINE      /* <-- EDIT for B7               */

#define F0_HZ               1000U          /* <-- EDIT: fundamental, Hz     */
#define DAC_AMPLITUDE       1200U          /* <-- EDIT: peak, in DAC codes  */

/* --- Fixed system parameters. Do not change between builds. -------------- */
#define SAMPLE_RATE_HZ      100000U        /* fs = 100 kSa/s                */
#define BUF_LEN             1000U          /* N = 1000 samples              */
#define DAC_MIDPOINT        2048U          /* mid-scale of the 12-bit DAC   */
#define TIMER_CLK_HZ        8000000U       /* TIM6 clock, see note below    */

/* TIMER_CLK_HZ note:
 * SystemClock_Config below runs SYSCLK from HSI directly at 8 MHz with no
 * PLL, AHB /1 and APB1 /1. TIM6 therefore sees 8 MHz. This matches the
 * configuration used for the power and frequency-accuracy measurements and
 * keeps digital switching noise low, which helps a THD measurement.
 * If you change the clock tree, change TIMER_CLK_HZ to match.
 */

/* ==========================================================================
 *  DERIVED CONSTANTS AND COMPILE-TIME CHECKS  --  do not edit
 * ========================================================================== */

/* Whole number of sine cycles held in the buffer. Must be an integer, or the
 * waveform will have a discontinuity every time the circular DMA wraps, and
 * that step generates its own harmonics that have nothing to do with the DAC. */
#define CYCLES_PER_BUFFER   ((F0_HZ * BUF_LEN) / SAMPLE_RATE_HZ)

/* TIM6 auto-reload for the required sample rate. Prescaler stays at 0. */
#define TIM6_PRESCALER      0U
#define TIM6_PERIOD         ((TIMER_CLK_HZ / SAMPLE_RATE_HZ) - 1U)

#if ((F0_HZ * BUF_LEN) % SAMPLE_RATE_HZ) != 0
#error "F0_HZ * BUF_LEN must divide exactly by SAMPLE_RATE_HZ. Choose an f0 that gives a whole number of cycles in the 1000-sample buffer (1, 2, 5, 10 kHz all work)."
#endif

#if (TIMER_CLK_HZ % SAMPLE_RATE_HZ) != 0
#error "TIMER_CLK_HZ must divide exactly by SAMPLE_RATE_HZ, otherwise fs is not exactly 100 kSa/s."
#endif

#if (DAC_MIDPOINT + DAC_AMPLITUDE) > 4095
#error "Amplitude too large: the sine would exceed DAC code 4095 and clip."
#endif

#if (DAC_AMPLITUDE > DAC_MIDPOINT)
#error "Amplitude too large: the sine would go below DAC code 0 and clip."
#endif

#if (DAC_AMPLITUDE > 1500)
#warning "Amplitude above 1500 codes drives the DAC output buffer close to the rails. Section 7.1 of the THD document says do not use 1800."
#endif
/* USER CODE END PD */

/* USER CODE BEGIN PV */
/* Lookup table. uint16_t because the DMA is configured for half-word
 * transfers with DAC_ALIGN_12B_R. See the DMA note in the header comment
 * of MX_DMA_Init in dma.c. */
static uint16_t sine_table[BUF_LEN];

/* These handles are created by CubeMX in dac.c, dma.c and tim.c. */
extern DAC_HandleTypeDef  hdac1;
extern DMA_HandleTypeDef  hdma_dac1_ch1;
extern TIM_HandleTypeDef  htim6;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */
static void Build_Sine_Table(void);
static void Start_Output(void);
/* USER CODE END PFP */

/* USER CODE BEGIN 0 */

/**
  * @brief  Fill the lookup table with one or more whole sine cycles.
  *
  *         sample[n] = 2048 + A * sin(2*pi * CYCLES * n / N)
  *
  *         The table is centred on code 2048 and holds exactly
  *         CYCLES_PER_BUFFER complete cycles, so the circular DMA wraps
  *         with no phase step.
  */
static void Build_Sine_Table(void)
{
    const double two_pi = 6.283185307179586476925286766559;
    double phase_step   = two_pi * (double)CYCLES_PER_BUFFER / (double)BUF_LEN;

    for (uint32_t n = 0U; n < BUF_LEN; n++)
    {
        double s = sin(phase_step * (double)n);

        /* Round to nearest code rather than truncating. Truncation adds a
         * small systematic error to every sample and shows up as extra
         * harmonic content in the very measurement we are trying to make. */
        double v = (double)DAC_MIDPOINT + (double)DAC_AMPLITUDE * s;
        int32_t code = (int32_t)(v + 0.5);

        if (code < 0)    { code = 0;    }
        if (code > 4095) { code = 4095; }

        sine_table[n] = (uint16_t)code;
    }
}

/**
  * @brief  Configure TIM6 for exactly fs, then start the DAC in the mode
  *         selected by OUTPUT_MODE.
  */
static void Start_Output(void)
{
#if (OUTPUT_MODE == MODE_OFF)

    /* ---- T9, instrument floor -------------------------------------------
     * DAC and DMA are left stopped. PA4 stays in analog mode from
     * MX_GPIO_Init / MX_DAC1_Init, so the pin is high impedance and drives
     * nothing. The LPF input therefore floats, which is what you want if
     * you are measuring the oscilloscope and pickup floor only.
     * -------------------------------------------------------------------- */
    (void)sine_table;

#elif (OUTPUT_MODE == MODE_DC_MIDPOINT)

    /* ---- T9, system floor -----------------------------------------------
     * The DAC drives a steady code 2048. The output buffer and the op-amp
     * run at the same operating point as during a real measurement, but no
     * signal is present. With AC coupling on the scope the DC disappears
     * and what remains is the noise floor of the whole chain.
     * This is the more useful floor for the 6 dB harmonic decision rule.
     * -------------------------------------------------------------------- */
    HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_1, DAC_ALIGN_12B_R, DAC_MIDPOINT);
    HAL_DAC_Start(&hdac1, DAC_CHANNEL_1);

#else   /* MODE_SINE */

    /* ---- Normal sine output --------------------------------------------- */

    /* Force the timer registers to the values derived above, so the sample
     * rate is correct even if CubeMX regenerates tim.c with different
     * numbers. Set before the timer is started. */
    __HAL_TIM_SET_PRESCALER(&htim6, TIM6_PRESCALER);
    __HAL_TIM_SET_AUTORELOAD(&htim6, TIM6_PERIOD);
    __HAL_TIM_SET_COUNTER(&htim6, 0U);

    /* Start the DAC on channel 1 with circular DMA from the lookup table.
     * The cast is required because the HAL prototype takes uint32_t*; the
     * DMA itself is configured for half-word transfers. */
    HAL_DAC_Start_DMA(&hdac1,
                      DAC_CHANNEL_1,
                      (uint32_t *)sine_table,
                      BUF_LEN,
                      DAC_ALIGN_12B_R);

    /* Mask the DMA half-transfer and transfer-complete interrupts.
     *
     * Circular DMA needs no callbacks. The 1000-sample buffer wraps every
     * 10 ms, so those interrupts would fire at 100 Hz, and harmonics of
     * 100 Hz land on 1 kHz, 2 kHz and 3 kHz - exactly the frequencies
     * being measured. Masking them keeps the core asleep and keeps those
     * current bursts out of the spectrum.
     *
     * This must come AFTER HAL_DAC_Start_DMA, because that function
     * enables the interrupts internally. */
    __HAL_DMA_DISABLE_IT(&hdma_dac1_ch1, DMA_IT_TC | DMA_IT_HT);

    /* Start the trigger source last, so the first sample goes out on a
     * clean timer update. */
    HAL_TIM_Base_Start(&htim6);

#endif
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  */
int main(void)
{
    /* MCU Configuration ---------------------------------------------------*/
    HAL_Init();

    SystemClock_Config();

    /* Initialise all configured peripherals.
     * DMA must be initialised before the DAC, or the DAC DMA request is
     * never connected. */
    MX_GPIO_Init();
    MX_DMA_Init();
    MX_DAC1_Init();
    MX_TIM6_Init();

    /* USER CODE BEGIN 2 */

    Build_Sine_Table();
    Start_Output();

    /* Stop the 1 ms SysTick interrupt.
     *
     * This matters for THD. SysTick fires at exactly 1 kHz, which is the
     * fundamental frequency of the primary test. Each interrupt is a small
     * burst of current that couples into the supply rails and can appear in
     * the spectrum at 1 kHz and its multiples, exactly where the harmonics
     * being measured live. With the DAC running entirely on hardware, no
     * tick is needed after start-up.
     *
     * Note: HAL_Delay will no longer work after this call. Nothing below
     * uses it. */
    HAL_SuspendTick();

    /* USER CODE END 2 */

    /* Infinite loop */
    while (1)
    {
        /* USER CODE BEGIN 3 */

        /* Put the core to sleep. The DAC, DMA and TIM6 keep running on
         * their own. A sleeping core switches far less, which lowers the
         * digital noise floor of the measurement. */
        __WFI();

        /* USER CODE END 3 */
    }
}

/**
  * @brief  System Clock Configuration
  *
  *         SYSCLK  = HSI = 8 MHz, no PLL
  *         AHB     = 8 MHz  (/1)
  *         APB1    = 8 MHz  (/1)  -> TIM6 clock = 8 MHz
  *         APB2    = 8 MHz  (/1)
  *         Flash latency 0 (valid up to 24 MHz)
  *
  *  WARNING: do not let CubeMX regenerate this function. The HSI PLL path
  *  on the F303 is ambiguous between HSI/2 and HSI/PREDIV, so the PLL is
  *  avoided entirely and HSI is used directly.
  */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef RCC_OscInitStruct = {0};
    RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

    /* Enable HSI and turn the PLL off. */
    RCC_OscInitStruct.OscillatorType      = RCC_OSCILLATORTYPE_HSI;
    RCC_OscInitStruct.HSIState            = RCC_HSI_ON;
    RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
    RCC_OscInitStruct.PLL.PLLState        = RCC_PLL_OFF;
    if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
    {
        Error_Handler();
    }

    /* Take SYSCLK straight from HSI. All bus prescalers at /1. */
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

/**
  * @brief  This function is executed in case of error occurrence.
  */
void Error_Handler(void)
{
    __disable_irq();
    while (1)
    {
    }
}

#ifdef  USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
    /* User can add his own implementation to report the file name and line
       number, e.g. printf("Wrong parameters value: file %s on line %d\r\n",
       file, line) */
}
#endif /* USE_FULL_ASSERT */
