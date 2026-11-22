/**
  ******************************************************************************
  * @file    main.c
  * @brief   Baseband Phase Modulation - phase deviation and phase accuracy
  *          measurement firmware.
  *
  *          Target : STM32 Nucleo-F303RE
  *          Clock  : HSE bypass 8 MHz (ST-LINK MCO) -> PLL x9 -> 72 MHz
  *          Output : DAC1_OUT1 on PA4, TIM6 TRGO, circular DMA
  *          Sync   : PB0 pulse at every DMA buffer wrap, 40 us wide
  *          Button : PC13 (B1) steps through the kp list
  *          UART   : USART2 (PA2/PA3) -> ST-LINK VCP, 115200 8N1
  *                   (VCP output is garbled on this board - verification is
  *                    done on the oscilloscope instead, see notes below)
  * ---------------------------------------------------------------------------
  * R1 REFERENCE CAPTURE
  * ---------------------------------------------------------------------------
  * This file boots at kp_idx = 0, which is kp = 0.0 rad, an unmodulated 5 kHz
  * carrier. The virtual COM port on this board produces garbled output, so the
  * kp index cannot be read back over UART. Either:
  *   - set kp_idx below and reflash for each sweep point (recommended: no
  *     ambiguity about which kp is loaded), or
  *   - press B1 to step, and count the LD2 blinks: one blink = index 0.
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
  * This file DOES define, because CubeMX did not generate them:
  *     huart2, MX_USART2_UART_Init, HAL_UART_MspInit
  *     MX_ExtraGPIO_Init   (PB0 sync, PA5 LD2, PC13 button interrupt)
  *     EXTI15_10_IRQHandler
  ******************************************************************************
  */

#include "main.h"
#include "dac.h"
#include "dma.h"
#include "gpio.h"
#include "tim.h"

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/* ==========================================================================
 *                                CONFIG
 * ========================================================================== */
#define MODE_ANALOG_PM   1          /* 1 = analog PM (PH1, PH2)                */
                                    /* 0 = BPSK (PH3)                         */

#define FC_HZ            5000.0f    /* baseband carrier                      */
#define FM_HZ            500.0f     /* message. Use 1000.0f for M2 only      */
#define FS_HZ            100000.0f  /* DAC update rate                       */
#define NSAMP            1000u      /* Use 4000 for BPSK at Rs = 250         */
#define RS_SYM           250.0f     /* BPSK symbol rate                      */

#define AMP_CODES        1200       /* peak amplitude in DAC codes           */
#define MID_CODE         2048       /* 2048 +/- 1200 = 848..3248 = 0.68..2.62 V */

#define SYNC_PULSE_US    40u        /* sync pulse width. Raise to 60 or 80   */
                                    /* if the trigger is still not stable.   */

#define TIM6_ARR_EXPECT  719u       /* checked against tim.c in the banner   */
#define TIM6_PSC_EXPECT  0u

#define SYNC_PORT        GPIOB
#define SYNC_PIN         GPIO_PIN_0
/* ========================================================================== */

/* kp sweep list. Index 0 is kp = 0, which is the R1 reference capture.
 * Index 4 is kp = 1.0, which is the R4 envelope check.
 * For the PH2 Bessel null, replace this list with the coarse or fine grid
 * and set FM_HZ to 1000.0f.                                                 */
static const float kp_list[] = {
    0.00f, 0.25f, 0.50f, 0.75f, 1.00f,
    1.50f, 2.00f, 2.50f, 3.00f, 3.50f, 4.00f
};
#define KP_COUNT  (sizeof(kp_list) / sizeof(kp_list[0]))

/* Set this to select the sweep point, then rebuild and flash.
 * 0 = 0.00 (R1)   4 = 1.00 (R4)   10 = 4.00                                 */
static volatile uint8_t kp_idx = 0;

/* BPSK pattern. 10 bits, so NSAMP must be 4000 at Rs = 250.
 * Unused while MODE_ANALOG_PM is 1 - the compiler warning is expected.      */
static const uint8_t bits[10] = { 0, 1, 0, 1, 1, 0, 1, 0, 0, 1 };
#define NBITS  (sizeof(bits) / sizeof(bits[0]))

/* --------------------------------------------------------------------------
 * Handles defined by CubeMX elsewhere
 * -------------------------------------------------------------------------- */
extern DAC_HandleTypeDef hdac1;           /* dac.c */
extern DMA_HandleTypeDef hdma_dac1_ch1;   /* dac.c */
extern TIM_HandleTypeDef htim6;           /* tim.c */

/* Defined here - CubeMX did not generate a usart.c */
UART_HandleTypeDef huart2;

/* --------------------------------------------------------------------------
 * Globals
 * -------------------------------------------------------------------------- */
static uint16_t dac_buf[NSAMP];

static volatile uint8_t  rebuild_req   = 0;
static volatile uint32_t last_btn_tick = 0;
static volatile uint32_t wrap_count    = 0;
static volatile uint8_t  clamp_flag    = 0;

/* --------------------------------------------------------------------------
 * Prototypes
 * -------------------------------------------------------------------------- */
void SystemClock_Config(void);
static void MX_USART2_UART_Init(void);
static void MX_ExtraGPIO_Init(void);
static void DWT_Init(void);
static void build_buffer(float kp);
static void blink_index(uint8_t n);
static void uprintf(const char *fmt, ...);
void Error_Handler(void);

/* ==========================================================================
 *                    Cycle-accurate microsecond delay
 * ==========================================================================
 * Uses the Cortex-M4 DWT cycle counter rather than a NOP loop, so the pulse
 * width is exact and repeatable regardless of optimisation level. The
 * subtraction is unsigned, so it remains correct across counter wrap.
 * ========================================================================== */
static void DWT_Init(void)
{
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL  |= DWT_CTRL_CYCCNTENA_Msk;
}

static inline void delay_us(uint32_t us)
{
    uint32_t start = DWT->CYCCNT;
    uint32_t ticks = us * (SystemCoreClock / 1000000UL);
    while ((DWT->CYCCNT - start) < ticks) { }
}

/* ==========================================================================
 *                              UART helper
 * ========================================================================== */
static void uprintf(const char *fmt, ...)
{
    char    buf[160];
    va_list ap;

    va_start(ap, fmt);
    int n = vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);

    if (n > 0) {
        if (n > (int)sizeof(buf) - 1) {
            n = (int)sizeof(buf) - 1;
        }
        HAL_UART_Transmit(&huart2, (uint8_t *)buf, (uint16_t)n, HAL_MAX_DELAY);
    }
}

/* ==========================================================================
 *                    LD2 blink, so kp index is readable
 * ==========================================================================
 * The VCP is unusable on this board, so the index is signalled by blinking
 * LD2: one blink = index 0, two blinks = index 1, and so on.
 * ========================================================================== */
static void blink_index(uint8_t n)
{
    for (uint8_t i = 0; i <= n; i++) {
        HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, GPIO_PIN_SET);
        HAL_Delay(120);
        HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, GPIO_PIN_RESET);
        HAL_Delay(180);
    }
}

/* ==========================================================================
 *                          Waveform generation
 * ==========================================================================
 * The DMA is stopped before the buffer is touched. Writing into a live
 * circular buffer would create a phase discontinuity in the middle of a
 * record, which shows up as a large spurious phase error.
 * ========================================================================== */
static void build_buffer(float kp)
{
    HAL_DAC_Stop_DMA(&hdac1, DAC_CHANNEL_1);

    clamp_flag = 0;

    for (uint32_t n = 0; n < NSAMP; n++) {
        float t  = (float)n / FS_HZ;
        float ph = 2.0f * (float)M_PI * FC_HZ * t;

#if MODE_ANALOG_PM
        ph += kp * sinf(2.0f * (float)M_PI * FM_HZ * t);
#else
        uint32_t k = (uint32_t)(t * RS_SYM) % NBITS;
        ph += (bits[k] ? (float)M_PI : 0.0f);
        (void)kp;
#endif

        int32_t code = MID_CODE + (int32_t)lrintf((float)AMP_CODES * cosf(ph));

        if (code < 0)    { code = 0;    clamp_flag = 1; }
        if (code > 4095) { code = 4095; clamp_flag = 1; }

        dac_buf[n] = (uint16_t)code;
    }

    if (HAL_DAC_Start_DMA(&hdac1, DAC_CHANNEL_1,
                          (uint32_t *)dac_buf, NSAMP,
                          DAC_ALIGN_12B_R) != HAL_OK) {
        Error_Handler();
    }
}

/* ==========================================================================
 *                              Callbacks
 * ==========================================================================
 * Sync pulse on PB0 at every buffer wrap, SYNC_PULSE_US wide. Because the
 * buffer holds an exact whole number of carrier cycles, message periods and
 * symbols, this edge is a known time origin and the only trigger that
 * anchors the same point of the message envelope on every acquisition.
 *
 * The DAC timing is not disturbed by the blocking delay: TIM6 triggers the
 * DMA in hardware, so interrupt duration does not reach the output.
 * ========================================================================== */
void HAL_DAC_ConvCpltCallbackCh1(DAC_HandleTypeDef *hdac)
{
    (void)hdac;

    SYNC_PORT->BSRR = SYNC_PIN;                          /* high */
    delay_us(SYNC_PULSE_US);
    SYNC_PORT->BSRR = (uint32_t)SYNC_PIN << 16u;         /* low  */

    wrap_count++;
}

void HAL_GPIO_EXTI_Callback(uint16_t pin)
{
    if (pin != GPIO_PIN_13) {
        return;
    }

    uint32_t now = HAL_GetTick();
    if ((now - last_btn_tick) < 250u) {                  /* debounce */
        return;
    }
    last_btn_tick = now;

    kp_idx      = (uint8_t)((kp_idx + 1u) % KP_COUNT);
    rebuild_req = 1;                                     /* do it in main() */
}

/* ==========================================================================
 *                                 main
 * ========================================================================== */
int main(void)
{
    HAL_Init();
    SystemClock_Config();
    DWT_Init();              /* must precede any delay_us call              */

    MX_GPIO_Init();          /* gpio.c - CubeMX pins                        */
    MX_DMA_Init();           /* dma.c  - must precede the DAC init          */
    MX_DAC1_Init();          /* dac.c                                       */
    MX_TIM6_Init();          /* tim.c                                       */
    MX_USART2_UART_Init();   /* here                                        */
    MX_ExtraGPIO_Init();     /* here - last, so PC13 EXTI is not overwritten */

    /* TIM6 clock is PCLK1 x2 because APB1 is prescaled. */
    uint32_t tim6_clk = HAL_RCC_GetPCLK1Freq() * 2u;
    float    fs_meas  = (float)tim6_clk
                        / (float)((htim6.Init.Prescaler + 1u)
                                  * (htim6.Init.Period + 1u));

    uprintf("\r\n=== PM phase metric firmware ===\r\n");
    uprintf("SYSCLK    : %lu Hz\r\n", (unsigned long)HAL_RCC_GetSysClockFreq());
    uprintf("TIM6 ARR  : %lu   (expect %u)\r\n",
            (unsigned long)htim6.Init.Period, (unsigned)TIM6_ARR_EXPECT);
    uprintf("fs actual : %.3f Sa/s\r\n", (double)fs_meas);
    uprintf("fc        : %.1f Hz  (%.1f cycles/buffer)\r\n",
            (double)FC_HZ, (double)(FC_HZ * NSAMP / FS_HZ));
#if MODE_ANALOG_PM
    uprintf("mode      : ANALOG PM, fm = %.1f Hz (%.1f periods/buffer)\r\n",
            (double)FM_HZ, (double)(FM_HZ * NSAMP / FS_HZ));
#else
    uprintf("mode      : BPSK, Rs = %.1f sym/s (%.1f symbols/buffer)\r\n",
            (double)RS_SYM, (double)(RS_SYM * NSAMP / FS_HZ));
#endif
    uprintf("sync pulse: %lu us wide, every %.3f ms\r\n",
            (unsigned long)SYNC_PULSE_US, (double)(1000.0f * NSAMP / FS_HZ));
    uprintf("amplitude : %d codes about %d  (%d..%d)\r\n",
            AMP_CODES, MID_CODE, MID_CODE - AMP_CODES, MID_CODE + AMP_CODES);

    if (htim6.Init.Period != TIM6_ARR_EXPECT ||
        htim6.Init.Prescaler != TIM6_PSC_EXPECT) {
        uprintf("*** WARNING: TIM6 config does not match. ***\r\n");
    }
    if (HAL_RCC_GetSysClockFreq() != 72000000u) {
        uprintf("*** WARNING: SYSCLK is not 72 MHz. HSE may not have "
                "started. ***\r\n");
    }
    if (hdma_dac1_ch1.Init.Mode != DMA_CIRCULAR) {
        uprintf("*** WARNING: DMA is not CIRCULAR. ***\r\n");
    }
    uprintf("--------------------------------\r\n");

    if (HAL_TIM_Base_Start(&htim6) != HAL_OK) {
        Error_Handler();
    }

    build_buffer(kp_list[kp_idx]);
    uprintf("kp = %.4f rad   [index %u/%u]\r\n",
            (double)kp_list[kp_idx], (unsigned)kp_idx,
            (unsigned)(KP_COUNT - 1u));

    /* Signal the loaded index on LD2, since the VCP is unreadable. */
    blink_index(kp_idx);

    while (1) {
        if (rebuild_req) {
            rebuild_req = 0;

            build_buffer(kp_list[kp_idx]);

            uprintf("kp = %.4f rad   [index %u/%u]%s\r\n",
                    (double)kp_list[kp_idx],
                    (unsigned)kp_idx,
                    (unsigned)(KP_COUNT - 1u),
                    clamp_flag ? "   *** CLAMPED ***" : "");

            blink_index(kp_idx);
        }
        HAL_Delay(10);
    }
}

/* ==========================================================================
 *                           Clock configuration
 * ==========================================================================
 * HSE on the Nucleo-64 is the 8 MHz MCO from the ST-LINK, not a crystal, so
 * HSE_BYPASS is correct. SB54 and SB55 must be fitted and SB16/SB50 removed.
 * On the STM32F3 the prediv lives inside the PLL struct, not beside HSEState.
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
    clk.APB1CLKDivider = RCC_HCLK_DIV2;       /* PCLK1 = 36 MHz, TIM6 = 72 MHz */
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
 * ========================================================================== */
static void MX_ExtraGPIO_Init(void)
{
    GPIO_InitTypeDef g = {0};

    __HAL_RCC_GPIOA_CLK_ENABLE();
    __HAL_RCC_GPIOB_CLK_ENABLE();
    __HAL_RCC_GPIOC_CLK_ENABLE();

    /* PB0 - sync marker to oscilloscope CH2 */
    HAL_GPIO_WritePin(SYNC_PORT, SYNC_PIN, GPIO_PIN_RESET);
    g.Pin   = SYNC_PIN;
    g.Mode  = GPIO_MODE_OUTPUT_PP;
    g.Pull  = GPIO_NOPULL;
    g.Speed = GPIO_SPEED_FREQ_HIGH;
    HAL_GPIO_Init(SYNC_PORT, &g);

    /* PA5 - LD2, blinks the kp index */
    HAL_GPIO_WritePin(GPIOA, GPIO_PIN_5, GPIO_PIN_RESET);
    g.Pin   = GPIO_PIN_5;
    g.Mode  = GPIO_MODE_OUTPUT_PP;
    g.Pull  = GPIO_NOPULL;
    g.Speed = GPIO_SPEED_FREQ_LOW;
    HAL_GPIO_Init(GPIOA, &g);

    /* PC13 - B1 user button, pulled up, active low */
    g.Pin  = GPIO_PIN_13;
    g.Mode = GPIO_MODE_IT_FALLING;
    g.Pull = GPIO_PULLUP;
    HAL_GPIO_Init(GPIOC, &g);

    HAL_NVIC_SetPriority(EXTI15_10_IRQn, 5, 0);
    HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);
}

/* ==========================================================================
 *   EXTI handler - CubeMX did not generate one, since PC13 is not in the .ioc
 * ========================================================================== */
void EXTI15_10_IRQHandler(void)
{
    HAL_GPIO_EXTI_IRQHandler(GPIO_PIN_13);
}

/* ==========================================================================
 *                            Error handler
 * ========================================================================== */
void Error_Handler(void)
{
    __disable_irq();
    while (1) {
        /* Board dark and silent: HSE most likely failed to start. */
    }
}

#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
    (void)file;
    (void)line;
}
#endif
