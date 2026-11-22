/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    : main.c
  * @brief   : Mode M3 - PM active, waveform recomputed at run time
  *
  *            Thesis: Implementation and Performance Analysis of a Baseband
  *            Phase Modulation (PM) System with STM32
  *
  *            PURPOSE
  *            -------
  *            M3 produces exactly the same output waveform as M1, but the
  *            samples are recomputed continuously at run time inside the DMA
  *            half-transfer and transfer-complete callbacks instead of being
  *            precomputed once at startup.
  *
  *            The difference M3 - M1 quantifies the cost of real-time
  *            floating-point waveform generation, and therefore justifies
  *            the design decision to use a precomputed buffer with circular
  *            DMA.
  *
  *            HOW IT WORKS
  *            ------------
  *            The DMA operates on the buffer in circular mode. While it is
  *            reading one half, the CPU rewrites the other half:
  *
  *              Half-transfer interrupt -> DMA now reads the second half,
  *                                         so the CPU refills the first half
  *              Transfer-complete       -> DMA wraps to the first half,
  *                                         so the CPU refills the second half
  *
  *            Because the waveform is periodic over exactly DAC_BUF_SIZE
  *            samples, each block is refilled with the same sample indices
  *            it originally held. The recomputed values are numerically
  *            identical to the precomputed ones, so the output is identical
  *            to M1 and the comparison isolates CPU work alone.
  *
  *            CLOCK CONFIGURATION
  *            -------------------
  *            Clock source : HSE (8 MHz, bypass, from ST-LINK MCO)
  *            PLL          : 8 MHz x 9 = 72 MHz
  *            HCLK         : 72 MHz
  *            PCLK1        : 36 MHz   (APB1 timer clocks = 72 MHz)
  *            PCLK2        : 72 MHz
  *            Flash latency: 2 wait states
  *            CPU state    : Running, computing samples in DMA callbacks
  *
  *            TIMER CONFIGURATION (set in CubeMX, not here)
  *            --------------------------------------------
  *            TIM6 clock = 72 MHz, PSC = 0, ARR = 719
  *            Fs = 72 MHz / (0 + 1)(719 + 1) = 100.000 kSa/s
  *
  *            REQUIRED: the DMA interrupt must be enabled in the NVIC, or
  *            the callbacks never fire. Check the NVIC Settings tab for the
  *            DMA1 channel assigned to DAC1_CH1.
  *
  *            TIMING BUDGET
  *            -------------
  *            One half-buffer = 1000 samples = 10 ms at 100 kSa/s.
  *            Each sample needs one sinf(), one cosf() and float scaling.
  *            See the note in the source below on build optimisation.
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

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

#define PI              3.14159265358979323846f

#define FS              100000.0f   /* DAC sampling frequency = 100 kSa/s */
#define FC              5000.0f     /* Carrier frequency      = 5 kHz     */
#define FM              500.0f      /* Message frequency      = 500 Hz    */
#define KP              1.0f        /* Phase sensitivity      = 1.0 rad   */

#define DAC_MID         2048        /* Midpoint for 12-bit DAC            */
#define DAC_AMP         1000        /* Amplitude, keep below clipping     */

#define DAC_BUF_SIZE    2000        /* 100 carrier cycles, 10 message cycles */
#define HALF_BUF        (DAC_BUF_SIZE / 2)

/* USER CODE END PD */

/* Private variables ---------------------------------------------------------*/
/* USER CODE BEGIN PV */

uint16_t dac_buffer[DAC_BUF_SIZE];

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);

/* USER CODE BEGIN PFP */
static void Compute_Block(uint32_t start, uint32_t count);
/* USER CODE END PFP */

/* USER CODE BEGIN 0 */

/**
  * @brief  Compute a block of phase-modulated samples at run time.
  *
  *         s(t) = Ac * cos(2*pi*fc*t + kp*m(t))
  *         m(t) = sin(2*pi*fm*t)
  *
  * @param  start  index of the first sample of the block
  * @param  count  number of samples to compute
  */
static void Compute_Block(uint32_t start, uint32_t count)
{
    for (uint32_t i = 0; i < count; i++)
    {
        uint32_t n = start + i;

        /* Sample time. At Fs = 100 kSa/s each sample is 10 us apart. */
        float t = (float)n / FS;

        /* Sinusoidal message signal */
        float message = sinf(2.0f * PI * FM * t);

        /* Instantaneous phase: carrier term plus phase deviation */
        float phase = (2.0f * PI * FC * t) + (KP * message);

        /* Phase-modulated output, -1 to +1 */
        float y = cosf(phase);

        /* Map to the 12-bit DAC range, centred on the midpoint */
        float dac_value = (float)DAC_MID + ((float)DAC_AMP * y);

        /* Safety clipping */
        if (dac_value < 0.0f)    { dac_value = 0.0f; }
        if (dac_value > 4095.0f) { dac_value = 4095.0f; }

        dac_buffer[n] = (uint16_t)dac_value;
    }
}

/**
  * @brief  DMA half-transfer callback.
  *         The DMA has finished the first half and is now reading the
  *         second half, so the first half is free to be rewritten.
  */
void HAL_DAC_ConvHalfCpltCallbackCh1(DAC_HandleTypeDef *hdac)
{
    Compute_Block(0, HALF_BUF);
}

/**
  * @brief  DMA transfer-complete callback.
  *         The DMA has wrapped back to the first half, so the second half
  *         is free to be rewritten.
  */
void HAL_DAC_ConvCpltCallbackCh1(DAC_HandleTypeDef *hdac)
{
    Compute_Block(HALF_BUF, HALF_BUF);
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
    /* Reset of all peripherals, initialise the Flash interface and SysTick. */
    HAL_Init();

    /* Configure the system clock: HSE + PLL -> 72 MHz */
    SystemClock_Config();

    /* Initialise all configured peripherals */
    MX_DMA_Init();
    MX_GPIO_Init();
    MX_DAC1_Init();
    MX_TIM6_Init();

    /* USER CODE BEGIN 2 */

    /* Prime the whole buffer once so the very first playback is valid.
     * From this point on the CPU refills it continuously in the callbacks.
     */
    Compute_Block(0, DAC_BUF_SIZE);

    /* Start the DAC in circular DMA mode. */
    if (HAL_DAC_Start_DMA(&hdac1,
                          DAC_CHANNEL_1,
                          (uint32_t *)dac_buffer,
                          DAC_BUF_SIZE,
                          DAC_ALIGN_12B_R) != HAL_OK)
    {
        Error_Handler();
    }

    /* Start TIM6. Its TRGO event triggers the DAC at 100 kSa/s. */
    if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
    {
        Error_Handler();
    }

    /* USER CODE END 2 */

    /* Infinite loop */
    while (1)
    {
        /* No sleep. The CPU is doing real work in the DMA callbacks; between
         * callbacks it spins here. This is the condition being measured.
         */
    }
}

/**
  * @brief  System Clock Configuration
  *
  *         HSE  = 8 MHz, supplied as a clock signal by the on-board
  *                ST-LINK MCO output. BYPASS mode is required because
  *                no crystal is fitted for HSE on the Nucleo board.
  *
  *         PLL  = HSE / PREDIV(1) x PLLMUL(9) = 72 MHz
  *
  *         Note: on the STM32F3 HAL, PREDIV is a member of the PLL
  *         sub-structure, not of RCC_OscInitTypeDef directly.
  *
  * @retval None
  */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef RCC_OscInitStruct = {0};
    RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

    RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
    RCC_OscInitStruct.HSEState       = RCC_HSE_BYPASS;
    RCC_OscInitStruct.PLL.PLLState   = RCC_PLL_ON;
    RCC_OscInitStruct.PLL.PLLSource  = RCC_PLLSOURCE_HSE;
    RCC_OscInitStruct.PLL.PREDIV     = RCC_PREDIV_DIV1;
    RCC_OscInitStruct.PLL.PLLMUL     = RCC_PLL_MUL9;

    if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
    {
        Error_Handler();
    }

    RCC_ClkInitStruct.ClockType      = RCC_CLOCKTYPE_HCLK
                                     | RCC_CLOCKTYPE_SYSCLK
                                     | RCC_CLOCKTYPE_PCLK1
                                     | RCC_CLOCKTYPE_PCLK2;
    RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_PLLCLK;
    RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;   /* HCLK  = 72 MHz */
    RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;     /* PCLK1 = 36 MHz */
    RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;     /* PCLK2 = 72 MHz */

    /* FLASH_LATENCY_2 is mandatory at 72 MHz. */
    if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
    {
        Error_Handler();
    }
}

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
    __disable_irq();
    while (1)
    {
    }
}

#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
}
#endif /* USE_FULL_ASSERT */
