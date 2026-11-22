/**
  ******************************************************************************
  * @file    main.c
  * @brief   R3 NOISE FLOOR reference capture.
  *
  *          Baseband Phase Modulation - phase deviation and phase accuracy
  *          measurement, reference capture R3.
  *
  *          Target : STM32 Nucleo-F303RE
  *          Clock  : HSE bypass 8 MHz (ST-LINK MCO) -> PLL x9 -> 72 MHz
  *          Output : DAC1_OUT1 on PA4, held STATIC at code 2048
  *
  * ---------------------------------------------------------------------------
  * WHAT THIS BUILD DOES DIFFERENTLY
  * ---------------------------------------------------------------------------
  * The DAC is driven by a single HAL_DAC_SetValue call and then left alone:
  *
  *     - no DMA          (no buffer transfers, so no transfer noise)
  *     - no TIM6         (timer left stopped, so no switching noise at fs)
  *     - no sync pulse   (PB0 is initialised low and never toggled)
  *     - no button       (nothing to step through)
  *
  * The captured record therefore contains only the noise of the DAC output
  * stage, the Sallen-Key filter, the wiring, and the oscilloscope front end.
  * Its RMS value in codes is the measurement floor quoted in the results.
  *
  * ---------------------------------------------------------------------------
  * BENCH PROCEDURE
  * ---------------------------------------------------------------------------
  *   1. Flash this build FIRST, before the main phase firmware, so the whole
  *      reference set needs only two flashes.
  *   2. Probes exactly as for R1: CH1 on the LPF output, ground clip on
  *      LPF GND. CH2 may stay on PB0 but will show a flat low line.
  *   3. DO NOT TOUCH THE VERTICAL SCALE. The noise floor is only meaningful
  *      at the same volts/div used for the measurement captures. If you
  *      rescale, the floor you compute will not apply to your data.
  *   4. Timebase 1 ms/div, same as every other phase capture.
  *   5. There is no sync pulse, so:
  *          trigger mode  = Auto
  *          acquisition   = Sample   (NOT Average - nothing to average to)
  *   6. CH1 coupling AC, as for the analysis captures.
  *   7. Capture and export as ref_noise.csv
  *
  * Expected on screen: a nearly flat trace with visible noise. If you see a
  * waveform, the DMA did not stop - check that HAL_TIM_Base_Start is still
  * commented out below.
  *
  * ---------------------------------------------------------------------------
  * PROJECT STRUCTURE
  * ---------------------------------------------------------------------------
  * Peripheral init lives in the CubeMX files. This file defines none of:
  *     hdac1, hdma_dac1_ch1, HAL_DAC_MspInit, MX_DAC1_Init   -> dac.c
  *     MX_DMA_Init                                           -> dma.c
  *     MX_GPIO_Init                                          -> gpio.c
  *     htim6, HAL_TIM_Base_MspInit, MX_TIM6_Init             -> tim.c
  *     DMA1_Channel3_IRQHandler, TIM6_DAC_IRQHandler         -> stm32f3xx_it.c
  *     HAL_MspInit                                           -> stm32f3xx_hal_msp.c
  *
  * USART2 is initialised but nothing is printed. The virtual COM port on this
  * board produces garbled output, so verification is done on the oscilloscope
  * instead. USART2 is kept only so the LPF sees the same board activity as in
  * the other captures.
  ******************************************************************************
  */

#include "main.h"
#include "dac.h"
#include "dma.h"
#include "gpio.h"
#include "tim.h"

#include <math.h>

/* ==========================================================================
 *                                CONFIG
 * ========================================================================== */
#define MID_CODE         2048       /* DAC held here for the whole capture   */

#define SYNC_PORT        GPIOB
#define SYNC_PIN         GPIO_PIN_0
/* ========================================================================== */

/* --------------------------------------------------------------------------
 * Handles defined by CubeMX elsewhere
 * -------------------------------------------------------------------------- */
extern DAC_HandleTypeDef hdac1;           /* dac.c */
extern DMA_HandleTypeDef hdma_dac1_ch1;   /* dac.c */
extern TIM_HandleTypeDef htim6;           /* tim.c */

UART_HandleTypeDef huart2;                /* defined here */

/* --------------------------------------------------------------------------
 * Prototypes
 * -------------------------------------------------------------------------- */
void SystemClock_Config(void);
static void MX_USART2_UART_Init(void);
static void MX_ExtraGPIO_Init(void);
void Error_Handler(void);

/* ==========================================================================
 *                                 main
 * ========================================================================== */
int main(void)
{
    HAL_Init();
    SystemClock_Config();

    MX_GPIO_Init();          /* gpio.c - CubeMX pins                        */
    MX_DMA_Init();           /* dma.c  - clock only, no transfer started    */
    MX_DAC1_Init();          /* dac.c                                       */
    MX_TIM6_Init();          /* tim.c  - configured but NOT started         */
    MX_USART2_UART_Init();   /* here                                        */
    MX_ExtraGPIO_Init();     /* here - last                                 */

    /* ------------------------------------------------------------------
     * TIM6 is deliberately NOT started. Leaving it stopped keeps its
     * switching activity out of the noise record. Do not uncomment.
     * ------------------------------------------------------------------ */
    /* HAL_TIM_Base_Start(&htim6); */

    /* ------------------------------------------------------------------
     * Static DAC output. No DMA, no trigger, no repetition.
     *
     * The DAC channel was configured in dac.c with DAC_TRIGGER_T6_TRGO.
     * With TIM6 stopped there are no trigger events, so the value written
     * below would sit in the holding register and never transfer. The
     * trigger is therefore disabled here so the value reaches the output
     * immediately.
     * ------------------------------------------------------------------ */
    {
        DAC_ChannelConfTypeDef ch = {0};
        ch.DAC_Trigger      = DAC_TRIGGER_NONE;
        ch.DAC_OutputBuffer = DAC_OUTPUTBUFFER_ENABLE;
        if (HAL_DAC_ConfigChannel(&hdac1, &ch, DAC_CHANNEL_1) != HAL_OK) {
            Error_Handler();
        }
    }

    if (HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_1,
                         DAC_ALIGN_12B_R, MID_CODE) != HAL_OK) {
        Error_Handler();
    }

    if (HAL_DAC_Start(&hdac1, DAC_CHANNEL_1) != HAL_OK) {
        Error_Handler();
    }

    /* PB0 stays low for the whole capture - no sync pulse in this build. */
    HAL_GPIO_WritePin(SYNC_PORT, SYNC_PIN, GPIO_PIN_RESET);

    /* LD2 on steadily, so it is obvious at a glance that this is the
     * noise-floor build and not the phase firmware.                     */
    HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, GPIO_PIN_SET);

    while (1) {
        /* Nothing. Any activity here would add noise to the record. */
        __WFI();
    }
}

/* ==========================================================================
 *                           Clock configuration
 * ==========================================================================
 * Identical to the phase firmware. The clock configuration must match the
 * measurement captures, because the switching activity of the PLL and the
 * core is part of what sets the noise floor.
 *
 * HSE on the Nucleo-64 is the 8 MHz MCO from the ST-LINK, not a crystal, so
 * HSE_BYPASS is correct. On the STM32F3 the prediv lives inside the PLL
 * struct, not beside HSEState.
 * ========================================================================== */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef       osc = {0};
    RCC_ClkInitTypeDef       clk = {0};
    RCC_PeriphCLKInitTypeDef per = {0};

    osc.OscillatorType = RCC_OSCILLATORTYPE_HSE;
    osc.HSEState       = RCC_HSE_BYPASS;      /* 8 MHz from ST-LINK MCO */
    osc.PLL.PLLState   = RCC_PLL_ON;
    osc.PLL.PLLSource  = RCC_PLLSOURCE_HSE;
    osc.PLL.PREDIV     = RCC_PREDIV_DIV1;     /* 8 MHz / 1 */
    osc.PLL.PLLMUL     = RCC_PLL_MUL9;        /* x 9 = 72 MHz */
    if (HAL_RCC_OscConfig(&osc) != HAL_OK) {
        Error_Handler();
    }

    clk.ClockType      = RCC_CLOCKTYPE_HCLK   | RCC_CLOCKTYPE_SYSCLK
                       | RCC_CLOCKTYPE_PCLK1  | RCC_CLOCKTYPE_PCLK2;
    clk.SYSCLKSource   = RCC_SYSCLKSOURCE_PLLCLK;
    clk.AHBCLKDivider  = RCC_SYSCLK_DIV1;     /* HCLK  = 72 MHz */
    clk.APB1CLKDivider = RCC_HCLK_DIV2;       /* PCLK1 = 36 MHz */
    clk.APB2CLKDivider = RCC_HCLK_DIV1;       /* PCLK2 = 72 MHz */
    if (HAL_RCC_ClockConfig(&clk, FLASH_LATENCY_2) != HAL_OK) {
        Error_Handler();
    }

    per.PeriphClockSelection = RCC_PERIPHCLK_USART2;
    per.Usart2ClockSelection = RCC_USART2CLKSOURCE_PCLK1;
    if (HAL_RCCEx_PeriphCLKConfig(&per) != HAL_OK) {
        Error_Handler();
    }
}

/* ==========================================================================
 *                        USART2 - not generated by CubeMX
 * ========================================================================== */
static void MX_USART2_UART_Init(void)
{
    huart2.Instance                    = USART2;
    huart2.Init.BaudRate               = 115200;
    huart2.Init.WordLength             = UART_WORDLENGTH_8B;
    huart2.Init.StopBits               = UART_STOPBITS_1;
    huart2.Init.Parity                 = UART_PARITY_NONE;
    huart2.Init.Mode                   = UART_MODE_TX_RX;
    huart2.Init.HwFlowCtl              = UART_HWCONTROL_NONE;
    huart2.Init.OverSampling           = UART_OVERSAMPLING_16;
    huart2.Init.OneBitSampling         = UART_ONE_BIT_SAMPLE_DISABLED;
    huart2.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
    if (HAL_UART_Init(&huart2) != HAL_OK) {
        Error_Handler();
    }
}

void HAL_UART_MspInit(UART_HandleTypeDef *huart)
{
    GPIO_InitTypeDef g = {0};

    if (huart->Instance != USART2) {
        return;
    }

    __HAL_RCC_USART2_CLK_ENABLE();
    __HAL_RCC_GPIOA_CLK_ENABLE();

    /* PA2 = TX, PA3 = RX, routed to the ST-LINK virtual COM port */
    g.Pin       = GPIO_PIN_2 | GPIO_PIN_3;
    g.Mode      = GPIO_MODE_AF_PP;
    g.Pull      = GPIO_NOPULL;
    g.Speed     = GPIO_SPEED_FREQ_HIGH;
    g.Alternate = GPIO_AF7_USART2;
    HAL_GPIO_Init(GPIOA, &g);
}

/* ==========================================================================
 *      Extra GPIO - pins CubeMX does not know about. Call this LAST.
 * ==========================================================================
 * PC13 is configured as a plain input here, not an interrupt. There is no
 * kp list to step through in this build, and an unused EXTI would only add
 * a way for the noise record to be disturbed by an accidental button press.
 * ========================================================================== */
static void MX_ExtraGPIO_Init(void)
{
    GPIO_InitTypeDef g = {0};

    __HAL_RCC_GPIOA_CLK_ENABLE();
    __HAL_RCC_GPIOB_CLK_ENABLE();
    __HAL_RCC_GPIOC_CLK_ENABLE();

    /* PB0 - held low, no sync pulse in this build */
    HAL_GPIO_WritePin(SYNC_PORT, SYNC_PIN, GPIO_PIN_RESET);
    g.Pin   = SYNC_PIN;
    g.Mode  = GPIO_MODE_OUTPUT_PP;
    g.Pull  = GPIO_NOPULL;
    g.Speed = GPIO_SPEED_FREQ_HIGH;
    HAL_GPIO_Init(SYNC_PORT, &g);

    /* PA5 - LD2, on steadily to identify this build */
    HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, GPIO_PIN_RESET);
    g.Pin   = GPIO_PIN_5;
    g.Mode  = GPIO_MODE_OUTPUT_PP;
    g.Pull  = GPIO_NOPULL;
    g.Speed = GPIO_SPEED_FREQ_LOW;
    HAL_GPIO_Init(GPIOA, &g);

    /* PC13 - plain input, no interrupt */
    g.Pin  = GPIO_PIN_13;
    g.Mode = GPIO_MODE_INPUT;
    g.Pull = GPIO_PULLUP;
    HAL_GPIO_Init(GPIOC, &g);
}

/* ==========================================================================
 *                            Error handler
 * ========================================================================== */
void Error_Handler(void)
{
    __disable_irq();
    while (1) {
        /* LD2 stays off. If the board is dark, HSE most likely failed. */
    }
}

#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
    (void)file;
    (void)line;
}
#endif
