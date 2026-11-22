/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Frequency Accuracy Test Mode
  *                   STM32 Nucleo-F303RE, DAC1_OUT1 (PA4), TIM6 TRGO, DMA
  *                   plus MCO clock output on PA8 for T4
  *
  *  Preset for F4 (MCO, HSI). Change the CONFIG block for the other tests.
  *
  *  F4 routes the raw oscillator to PA8 so the clock can be measured on its
  *  own, without the DAC, DMA, timer or buffer in the path. Comparing the
  *  PA8 result with the DAC-output result from F1 and F2 turns the argument
  *  that the generation chain adds no frequency error into a measured fact.
  *
  *  Clock plan
  *  ----------
  *  HSI : HSI 8 MHz straight to SYSCLK, PLL OFF.
  *        APB1 /1  ->  PCLK1 8 MHz  ->  TIM6 clock 8 MHz
  *        Fs = 8 MHz / (0+1)(79+1) = 100 000 Sa/s exactly
  *
  *        The PLL is deliberately unused here. On the F303xE the HAL
  *        constant RCC_PLLSOURCE_HSI may resolve to HSI/2 or to HSI/PREDIV
  *        depending on library version. With HSI/PREDIV the PLL input is
  *        8 MHz and x16 gives 128 MHz, far past the 72 MHz limit; the MCU
  *        then runs unstably and the DAC output collapses. Core speed does
  *        not matter for this metric - only the sample rate does, and that
  *        is exactly 100 kSa/s either way.
  *
  *  HSE : HSE 8 MHz bypass (ST-LINK MCO), PREDIV /1, PLL x9 -> 72 MHz.
  *        APB1 /2  ->  PCLK1 36 MHz  ->  TIM6 clock 72 MHz (timer doubler)
  *        Fs = 72 MHz / (0+1)(719+1) = 100 000 Sa/s exactly
  *
  *  MCO note
  *  --------
  *  The F303 MCO has no prescaler, so PA8 carries the oscillator at its
  *  full 8 MHz in both cases. Measure it at 500 ns/div with the scope
  *  bandwidth limit OFF, and keep the probe ground lead short.
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
/* USER CODE END Includes */

/* USER CODE BEGIN PD */

/* ==========================================================================
 * CONFIG   -   the only lines that change between measurement runs
 * ========================================================================== */

#define CLOCK_SOURCE    CLOCK_HSI      /* CLOCK_HSI  or  CLOCK_HSE           */
#define TEST_FREQ_HZ    5000.0f        /* 1000 / 2000 / 5000 / 10000 only    */
#define KP              0.0f           /* 0.0f normally, 1.0f for F6         */
#define MCO_ENABLE      1              /* 1 for F4, 0 for every other test   */

/* ========================================================================== */

#define CLOCK_HSI       0
#define CLOCK_HSE       1

#define PI              3.14159265358979323846f

#define FS              100000UL       /* DAC sample rate, both clock paths  */
#define FM              500.0f         /* message frequency, used if KP > 0  */

#define DAC_MID         2048           /* 12-bit midpoint                    */
#define DAC_AMP         1200           /* 0.68 V to 2.62 V swing             */

#define DAC_BUF_SIZE    1000           /* N = 1000, buffer repeats at 100 Hz */

/*
 * Whole-cycle rule.
 *
 * The buffer holds 1000 samples at 100 kSa/s, so it repeats every 10 ms and
 * only integer multiples of 100 Hz fit a whole number of cycles.
 *
 *    1 kHz  ->  10 cycles, 100 samples per cycle
 *    2 kHz  ->  20 cycles,  50 samples per cycle
 *    5 kHz  ->  50 cycles,  20 samples per cycle
 *   10 kHz  -> 100 cycles,  10 samples per cycle
 */
#define N_CYCLES  ((uint32_t)((TEST_FREQ_HZ * (float)DAC_BUF_SIZE) / (float)FS))

/* USER CODE END PD */

/* USER CODE BEGIN PV */

uint16_t dac_buffer[DAC_BUF_SIZE];

/* Watch these in the debugger to confirm the clock actually switched. */
volatile uint32_t dbg_sysclk;      /* HSI: 8000000    HSE: 72000000 */
volatile uint32_t dbg_timclk;      /* HSI: 8000000    HSE: 72000000 */
volatile uint32_t dbg_arr;         /* HSI:      79    HSE:      719 */
volatile uint32_t dbg_ncycles;     /* 50 at 5 kHz                   */
volatile uint32_t dbg_mco;         /* frequency expected on PA8, Hz */

/* USER CODE END PV */

void SystemClock_Config(void);

/* USER CODE BEGIN PFP */
void Generate_Test_Buffer(void);
static uint32_t Get_Timer_Clock(void);
static void MCO_Setup(void);
/* USER CODE END PFP */

/* USER CODE BEGIN 0 */

/**
  * @brief Actual clock feeding TIM6.
  *
  * On the F3, if the APB1 prescaler is 1 the timer clock equals PCLK1;
  * otherwise the timer clock is twice PCLK1. Reading it back means the
  * auto-reload value is derived rather than hard-coded, so it can never
  * be left out of step with the clock source.
  */
static uint32_t Get_Timer_Clock(void)
{
    uint32_t pclk1 = HAL_RCC_GetPCLK1Freq();

    if ((RCC->CFGR & RCC_CFGR_PPRE1) != RCC_CFGR_PPRE1_DIV1)
    {
        pclk1 *= 2U;
    }
    return pclk1;
}

/**
  * @brief Route the selected oscillator to PA8.
  *
  * HAL_RCC_MCOConfig configures the GPIO itself, so nothing needs adding
  * to the .ioc. If the compiler does not recognise the source constant,
  * type RCC_MCO1SOURCE_ and let autocomplete list what this HAL version
  * provides; older trees use RCC_MCOSOURCE_ instead.
  */
static void MCO_Setup(void)
{
#if MCO_ENABLE
  #if (CLOCK_SOURCE == CLOCK_HSE)
    HAL_RCC_MCOConfig(RCC_MCO1, RCC_MCO1SOURCE_HSE, RCC_MCODIV_1);
  #else
    HAL_RCC_MCOConfig(RCC_MCO1, RCC_MCO1SOURCE_HSI, RCC_MCODIV_1);
  #endif
    dbg_mco = 8000000U;        /* no MCO prescaler on the F303 */
#else
    dbg_mco = 0U;
#endif
}

/**
  * @brief Fill the DAC buffer.
  *
  * With KP = 0 this is a pure cosine, so the zero crossings are undisplaced,
  * which is what the zero-crossing regression needs.
  *
  * The carrier phase is built from the integer cycle count rather than from
  * elapsed time, which guarantees the last sample joins the first with no
  * phase step regardless of floating point rounding.
  */
void Generate_Test_Buffer(void)
{
    uint32_t n;

    for (n = 0; n < DAC_BUF_SIZE; n++)
    {
        float carrier_phase =
            2.0f * PI * (float)((N_CYCLES * n) % DAC_BUF_SIZE)
                      / (float)DAC_BUF_SIZE;

        float phase = carrier_phase;

        if (KP != 0.0f)                       /* T6 only */
        {
            float t = (float)n / (float)FS;
            phase += KP * sinf(2.0f * PI * FM * t);
        }

        float y = cosf(phase);
        float dac_value = (float)DAC_MID + ((float)DAC_AMP * y);

        if (dac_value < 0.0f)    { dac_value = 0.0f;    }
        if (dac_value > 4095.0f) { dac_value = 4095.0f; }

        dac_buffer[n] = (uint16_t)(dac_value + 0.5f);
    }
}

/* USER CODE END 0 */

int main(void)
{
    /* USER CODE BEGIN 1 */
    uint32_t timclk;
    uint32_t arr;
    /* USER CODE END 1 */

    HAL_Init();

    SystemClock_Config();

    /* USER CODE BEGIN SysInit */
    MCO_Setup();
    /* USER CODE END SysInit */

    MX_DMA_Init();
    MX_GPIO_Init();
    MX_DAC1_Init();
    MX_TIM6_Init();

    /* USER CODE BEGIN 2 */

    /*
     * Derive PSC and ARR from the clock TIM6 is really running on. The
     * correct ARR differs between the two clock paths, and computing it
     * here removes the possibility of flashing a mismatched pair.
     */
    timclk = Get_Timer_Clock();

    if ((timclk % FS) != 0U)          /* must divide exactly */
    {
        Error_Handler();
    }

    arr = (timclk / FS) - 1U;

    __HAL_TIM_SET_PRESCALER(&htim6, 0U);
    __HAL_TIM_SET_AUTORELOAD(&htim6, arr);
    __HAL_TIM_SET_COUNTER(&htim6, 0U);

    dbg_sysclk  = HAL_RCC_GetSysClockFreq();
    dbg_timclk  = timclk;
    dbg_arr     = arr;
    dbg_ncycles = N_CYCLES;

    Generate_Test_Buffer();

    if (HAL_DAC_Start_DMA(&hdac1,
                          DAC_CHANNEL_1,
                          (uint32_t *)dac_buffer,
                          DAC_BUF_SIZE,
                          DAC_ALIGN_12B_R) != HAL_OK)
    {
        Error_Handler();
    }

    if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
    {
        Error_Handler();
    }

    /* USER CODE END 2 */

    while (1)
    {
        /*
         * Nothing to do. TIM6 TRGO and DMA drive PA4, and the MCO hardware
         * drives PA8, both without CPU involvement. Probing either pin
         * measures the same oscillator.
         */
    }
}

/**
  * @brief System Clock Configuration
  */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef RCC_OscInitStruct = {0};
    RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

#if (CLOCK_SOURCE == CLOCK_HSE)

    /* ---- HSE bypass 8 MHz, PREDIV /1, PLL x9 -> 72 MHz ---- */
    RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
    RCC_OscInitStruct.HSEState       = RCC_HSE_BYPASS;
    RCC_OscInitStruct.PLL.PLLState   = RCC_PLL_ON;
    RCC_OscInitStruct.PLL.PLLSource  = RCC_PLLSOURCE_HSE;
    RCC_OscInitStruct.PLL.PREDIV     = RCC_PREDIV_DIV1;
    RCC_OscInitStruct.PLL.PLLMUL     = RCC_PLL_MUL9;

#else

    /* ---- HSI 8 MHz used directly, PLL off ---- */
    RCC_OscInitStruct.OscillatorType      = RCC_OSCILLATORTYPE_HSI;
    RCC_OscInitStruct.HSIState            = RCC_HSI_ON;
    RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
    RCC_OscInitStruct.PLL.PLLState        = RCC_PLL_NONE;

#endif

    if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
    {
        Error_Handler();
    }

    RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK   |
                                  RCC_CLOCKTYPE_SYSCLK |
                                  RCC_CLOCKTYPE_PCLK1  |
                                  RCC_CLOCKTYPE_PCLK2;

#if (CLOCK_SOURCE == CLOCK_HSE)

    RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_PLLCLK;
    RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;
    RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;   /* PCLK1 max 36 MHz */
    RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

    if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)

#else

    RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_HSI;
    RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;
    RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;   /* 8 MHz to TIM6 */
    RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

    if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_0) != HAL_OK)

#endif
    {
        Error_Handler();
    }
}

void Error_Handler(void)
{
    __disable_irq();
    while (1)
    {
    }
}

#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
}
#endif
