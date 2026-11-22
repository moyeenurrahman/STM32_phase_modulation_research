/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    : main.c
  * @brief   : Mode P0 - Baseline power measurement condition
  *
  *            Thesis: Implementation and Performance Analysis of a Baseband
  *            Phase Modulation (PM) System with STM32
  *
  *            PURPOSE
  *            -------
  *            M0 establishes the static power floor of the microcontroller:
  *            a powered and clocked STM32F303RE with NO signal generation.
  *            All other power results are interpreted relative to this value.
  *
  *            The DAC, DMA and TIM6 peripherals are deliberately NOT
  *            initialised, so their peripheral clocks remain disabled and
  *            they contribute no current.
  *
  *            CONFIGURATION
  *            -------------
  *            Clock source : HSE (8 MHz, bypass, from ST-LINK MCO)
  *            PLL          : 8 MHz x 9 = 72 MHz
  *            HCLK         : 72 MHz
  *            PCLK1        : 36 MHz   (APB1 timer clocks = 72 MHz)
  *            PCLK2        : 72 MHz
  *            Flash latency: 2 wait states
  *            CPU state    : Sleep, via __WFI()
  *
  *            EXPECTED RESULT
  *            ---------------
  *            Oscilloscope: flat line. PA4 is not driven.
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "gpio.h"

/* NOTE: dac.h, dma.h and tim.h are intentionally NOT included.
 *       Those peripherals must stay unclocked in the baseline condition.
 */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* Reset of all peripherals, initialise the Flash interface and the SysTick. */
  HAL_Init();

  /* Configure the system clock: HSE + PLL -> 72 MHz */
  SystemClock_Config();

  /* Initialise GPIO only.
   *
   * MX_DMA_Init(), MX_DAC1_Init() and MX_TIM6_Init() are deliberately
   * omitted. The RCC peripheral-clock enable for each of those blocks
   * lives inside its own HAL_MspInit function, which is called from
   * those init functions. Skipping the calls therefore leaves DMA1,
   * DAC1 and TIM6 completely unclocked.
   */
  MX_GPIO_Init();

  /* Infinite loop */
  while (1)
  {
    /* Sleep mode. The CPU clock stops until an interrupt occurs.
     *
     * SysTick remains enabled (started by HAL_Init) and wakes the core
     * once per millisecond. This is identical in M1, M2 and M3, so the
     * comparison between modes is unaffected.
     */
    __WFI();
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
  *         TIM6 clock = 72 MHz, because the APB1 prescaler is greater
  *         than 1 and the timer clock is therefore doubled from PCLK1.
  *
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /* Configure HSE and the PLL */
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

  /* Configure the CPU, AHB and APB bus clocks */
  RCC_ClkInitStruct.ClockType      = RCC_CLOCKTYPE_HCLK
                                   | RCC_CLOCKTYPE_SYSCLK
                                   | RCC_CLOCKTYPE_PCLK1
                                   | RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource   = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider  = RCC_SYSCLK_DIV1;   /* HCLK  = 72 MHz */
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;     /* PCLK1 = 36 MHz */
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;     /* PCLK2 = 72 MHz */

  /* FLASH_LATENCY_2 is mandatory at 72 MHz. A lower setting causes a
   * hard fault as soon as the clock is switched to the PLL.
   */
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
