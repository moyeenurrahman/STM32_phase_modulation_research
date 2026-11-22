/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Frequency Accuracy Test Mode
  *                   STM32 Nucleo-F303RE, DAC1_OUT1 (PA4) + TIM6 TRGO + DMA
  *
  *  Covers test matrix entries:
  *    F7  chain check   : CLOCK_HSE, KP 0, 5 kHz, probe PA4 and LPF output
  *    F1  HSI sweep     : CLOCK_HSI, KP 0, 1 / 2 / 5 / 10 kHz
  *    F2  HSE sweep     : CLOCK_HSE, KP 0, 1 / 2 / 5 / 10 kHz
  *    F3  ratiometric   : alternate CLOCK_HSI and CLOCK_HSE at 5 kHz
  *    F6  PM check      : CLOCK_HSE, KP 1.0f, 5 kHz
  *
  *  Change ONLY the two defines in the CONFIG block between runs.
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
 * CONFIG  -  the only two lines you change between measurement runs
 * ========================================================================== */

#define CLOCK_SOURCE    CLOCK_HSE      /* CLOCK_HSI  or  CLOCK_HSE           */
#define TEST_FREQ_HZ    5000.0f        /* 1000 / 2000 / 5000 / 10000 only    */
#define KP              0.0f           /* 0.0f for F1,F2,F3,F7 - 1.0f for F6 */

/* ========================================================================== */

#define CLOCK_HSI       0
#define CLOCK_HSE       1

#define PI              3.14159265358979323846f

#define FS              100000.0f      /* DAC sample rate, both clock configs */
#define FM              500.0f         /* message frequency, used only if KP>0 */

#define DAC_MID         2048           /* 12-bit midpoint                     */
#define DAC_AMP         1200           /* 0.68 V to 2.62 V swing              */

#define DAC_BUF_SIZE    1000           /* N = 1000, buffer repeats at 100 Hz  */

/*
 * Clock plan.
 *
 *  HSI path : HSI 8 MHz used directly as SYSCLK, PLL OFF.
 *             APB1 /1 -> PCLK1 8 MHz -> TIM6 clock 8 MHz
 *             Fs = 8e6 / (0+1)(79+1) = 100 000 Sa/s  exactly
 *
 *             The PLL is deliberately not used on the HSI path. On the
 *             F303xE the HAL constant RCC_PLLSOURCE_HSI can select either
 *             HSI/2 or HSI/PREDIV depending on the library version. If it
 *             selects HSI/PREDIV the PLL input is 8 MHz, and x16 gives
 *             128 MHz, far beyond the 72 MHz limit. The chip then runs
 *             unstably and the DAC output collapses. Running HSI straight
 *             through avoids the ambiguity entirely, and system clock speed
 *             is irrelevant to this metric: only the sample rate matters,
 *             and that is exactly 100 kSa/s either way.
 *
 *  HSE path : HSE 8 MHz bypass (ST-LINK MCO) -> PREDIV /1 -> PLL x9 -> 72 MHz
 *             APB1 /2 -> PCLK1 36 MHz -> TIM6 clock 72 MHz
 *             Fs = 72e6 / (0+1)(719+1) = 100 000 Sa/s  exactly
 *
 *  Both paths are exact integer chains, so the fractional error of the DAC
 *  output equals the fractional error of the oscillator, and the ratiometric
 *  comparison in T3 stays valid despite the different multiplier and ARR.
 */
#if (CLOCK_SOURCE == CLOCK_HSE)
  #define TIM6_ARR_VALUE  719      /* 72 MHz timer clock */
#else
  #define TIM6_ARR_VALUE  79       /*  8 MHz timer clock */
#endif

/*
 * Whole-cycle rule.
 *
 * The buffer holds N = 1000 samples at 100 kSa/s, so it repeats every 10 ms.
 * Only integer multiples of 100 Hz fit a whole number of cycles. Anything
 * else produces a phase jump at every DMA wrap and corrupts the spectrum.
 *
 *   1 kHz  ->  10 cycles, 100 samples per cycle
 *   2 kHz  ->  20 cycles,  50 samples per cycle
 *   5 kHz  ->  50 cycles,  20 samples per cycle
 *  10 kHz  -> 100 cycles,  10 samples per cycle
 */
#define N_CYCLES  ((uint32_t)((TEST_FREQ_HZ * (float)DAC_BUF_SIZE) / FS))

/* USER CODE END PD */

/* USER CODE BEGIN PV */
uint16_t dac_buffer[DAC_BUF_SIZE];
/* USER CODE END PV */

void SystemClock_Config(void);

/* USER CODE BEGIN PFP */
void Generate_Test_Buffer(void);
/* USER CODE END PFP */

/* USER CODE BEGIN 0 */

/**
  * @brief Fill the DAC buffer.
  *
  * With KP = 0 this is a pure cosine and the zero crossings are undisplaced,
  * which is what the zero-crossing regression needs.
  *
  * The carrier phase is built from the integer cycle count rather than from
  * elapsed time. This guarantees the last sample joins the first sample with
  * no phase step, independently of floating point rounding.
  */
void Generate_Test_Buffer(void)
{
    uint32_t n;

    for (n = 0; n < DAC_BUF_SIZE; n++)
    {
        /* carrier phase, exact at the buffer boundary */
        float carrier_phase =
            2.0f * PI * (float)((N_CYCLES * n) % DAC_BUF_SIZE)
                      / (float)DAC_BUF_SIZE;

        float phase = carrier_phase;

        /* modulation, active only for T6 */
        if (KP != 0.0f)
        {
            float t = (float)n / FS;
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
    HAL_Init();

    SystemClock_Config();

    MX_DMA_Init();
    MX_GPIO_Init();
    MX_DAC1_Init();
    MX_TIM6_Init();

    /* USER CODE BEGIN 2 */

    /*
     * Force PSC and ARR here rather than trusting the CubeMX value.
     * The correct ARR differs between the 64 MHz and 72 MHz clock paths,
     * and setting it in code keeps the two configurations self-consistent.
     */
    __HAL_TIM_SET_PRESCALER(&htim6, 0);
    __HAL_TIM_SET_AUTORELOAD(&htim6, TIM6_ARR_VALUE);
    __HAL_TIM_SET_COUNTER(&htim6, 0);

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
        /* TIM6 TRGO and DMA run the output with no CPU involvement */
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

    /* ---- HSE bypass, 8 MHz from the ST-LINK MCO, PLL x9 -> 72 MHz ---- */
    RCC_OscInitStruct.OscillatorType   = RCC_OSCILLATORTYPE_HSE;
    RCC_OscInitStruct.HSEState         = RCC_HSE_BYPASS;
    RCC_OscInitStruct.PLL.PLLState     = RCC_PLL_ON;
    RCC_OscInitStruct.PLL.PLLSource    = RCC_PLLSOURCE_HSE;
    RCC_OscInitStruct.PLL.PREDIV       = RCC_PREDIV_DIV1;
    RCC_OscInitStruct.PLL.PLLMUL       = RCC_PLL_MUL9;

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

    RCC_ClkInitStruct.ClockType      = RCC_CLOCKTYPE_HCLK   |
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
