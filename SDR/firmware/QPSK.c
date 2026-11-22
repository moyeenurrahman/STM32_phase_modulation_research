/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : QPSK baseband / passband generator
  *
  *  Project : Implementation and Performance Analysis of a Baseband Phase
  *            Modulation System with STM32 for Low-Power Wireless Communication
  *  Board   : NUCLEO-F303RE
  *
  *  This file replaces the BPSK demo main.c.  Everything that changes between
  *  runs lives in the CONFIG block below.  No other file needs editing.
  *
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

/* USER CODE BEGIN PV */

/* ==========================================================================
 *  CONFIG BLOCK  -  edit only this section between runs
 * ========================================================================== */

/* Output mode:
 *   0 = PASSBAND  : composite QPSK on DAC1_OUT1 (PA4) only.
 *                   s(t) = I*cos(2*pi*fc*t) - Q*sin(2*pi*fc*t)
 *                   PA5 is not driven.  Safe on Nucleo (LD2 stays off).
 *
 *   1 = BASEBAND  : I on DAC1_OUT1 (PA4), Q on DAC1_OUT2 (PA5).
 *                   WARNING: on NUCLEO-F303RE, PA5 also drives the green
 *                   user LED LD2.  The LED clamps the waveform above about
 *                   1.8 V.  Either remove the LD2 solder bridge, or move Q
 *                   to PA6 (DAC2_OUT1), before trusting any measurement
 *                   taken in this mode.
 */
#define QPSK_MODE            0

/* Pulse shaping:
 *   0 = rectangular NRZ  (simple, matches the BPSK demo)
 *   1 = raised cosine    (band-limited, ISI-free at symbol instants)
 *
 * Note: with PULSE_SHAPE = 1 at an 8 MHz system clock, buffer generation
 * takes a few seconds at start-up.  This is normal, not a hang.
 */
#define PULSE_SHAPE          0
#define RC_BETA              0.35f   /* roll-off factor, 0 < beta <= 1        */
#define RC_HALF_SPAN_SYM     3       /* pulse truncated at +/- this many symbols */

/* Signal parameters (integers, so the compile-time checks below can run) */
#define FS_HZ_INT            100000U /* DAC sample rate,  samples per second  */
#define FC_HZ_INT            1000U   /* carrier frequency, Hz                 */
#define RS_HZ_INT            200U    /* symbol rate, symbols per second       */

/* Amplitude */
#define DAC_MID              2048U   /* 12-bit midpoint, about 1.65 V         */
#define DAC_AMP              800.0f  /* peak deviation in DAC codes           */

/* Pattern length in symbols */
#define PATTERN_LEN_SYM      8U

/* ========================================================================== */

#define PI                   3.14159265358979323846f

#define FS_HZ                ((float)FS_HZ_INT)
#define FC_HZ                ((float)FC_HZ_INT)
#define RS_HZ                ((float)RS_HZ_INT)

/* Samples per symbol: 100000 / 200 = 500 */
#define SPS                  (FS_HZ_INT / RS_HZ_INT)

/* Buffer length: 500 * 8 = 4000 samples = 40 ms */
#define BUF_LEN              (SPS * PATTERN_LEN_SYM)

/* --- Compile-time sanity checks -------------------------------------------
 * 1) The sample rate must divide evenly into symbols, or symbol boundaries
 *    drift away from sample boundaries.
 * 2) The buffer must hold a whole number of carrier cycles, or the circular
 *    DMA produces a phase step every time it wraps.
 */
_Static_assert((FS_HZ_INT % RS_HZ_INT) == 0U,
               "FS must be an integer multiple of RS");
_Static_assert(((FC_HZ_INT * PATTERN_LEN_SYM) % RS_HZ_INT) == 0U,
               "Buffer must contain a whole number of carrier cycles");

/* --- Constellation --------------------------------------------------------
 * Gray-coded QPSK.  Dibit d = (b1 b0):
 *
 *   d = 0  (00)  ->  I = +1, Q = +1  ->   45 degrees
 *   d = 1  (01)  ->  I = -1, Q = +1  ->  135 degrees
 *   d = 2  (10)  ->  I = +1, Q = -1  ->  315 degrees
 *   d = 3  (11)  ->  I = -1, Q = -1  ->  225 degrees
 *
 * Neighbouring points differ by exactly one bit, so a phase error that
 * lands in the wrong quadrant costs one bit, not two.
 */
static const float qpsk_i[4] = { +1.0f, -1.0f, +1.0f, -1.0f };
static const float qpsk_q[4] = { +1.0f, +1.0f, -1.0f, -1.0f };

/* Symbol pattern.
 * Visits all four constellation points, and contains both 90 degree
 * transitions and 180 degree transitions (0->3 and 1->2).  The 180 degree
 * transitions are the interesting ones: they pass through the origin.
 *
 * Bit stream: 00 01 11 10 00 11 01 10   (16 bits, 400 bps at Rs = 200)
 */
static const uint8_t symbol_pattern[PATTERN_LEN_SYM] =
    { 0U, 1U, 3U, 2U, 0U, 3U, 1U, 2U };

/* DMA buffers.  4000 samples * 2 bytes * 2 channels = 16 kB. */
static uint16_t dac_buf_i[BUF_LEN];
static uint16_t dac_buf_q[BUF_LEN];

/* Debug variables - inspect these in the CubeIDE Live Expressions window
 * to confirm the clock tree and the ARR actually loaded at run time.
 */
volatile uint32_t dbg_pclk1;
volatile uint32_t dbg_timclk;
volatile uint32_t dbg_arr;
volatile uint32_t dbg_buflen;
volatile float    dbg_peak;      /* peak of the normalised waveform  */
volatile float    dbg_scale;     /* codes per unit amplitude          */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);

/* USER CODE BEGIN PFP */
static void  Configure_TIM6_Runtime(void);
static float RaisedCosine(float x);
static void  Compute_IQ(uint32_t n, float *out_i, float *out_q);
static void  Build_Buffers(void);
/* USER CODE END PFP */

/* USER CODE BEGIN 0 */

/**
  * @brief  Set TIM6 PSC/ARR at run time from the real bus clock.
  *
  * Hard-coding ARR = 79 is only correct while the timer clock is 8 MHz.
  * Reading the clock back guarantees the DAC sample rate is exactly FS_HZ,
  * or the code stops in Error_Handler() instead of running at the wrong rate.
  */
static void Configure_TIM6_Runtime(void)
{
    uint32_t pclk1  = HAL_RCC_GetPCLK1Freq();
    uint32_t timclk = pclk1;
    uint32_t div;

    /* On STM32F3 the timer clock equals PCLK1 when the APB1 prescaler is 1,
     * and twice PCLK1 for any larger prescaler.
     */
    if ((RCC->CFGR & RCC_CFGR_PPRE1) != RCC_CFGR_PPRE1_DIV1)
    {
        timclk = pclk1 * 2U;
    }

    dbg_pclk1  = pclk1;
    dbg_timclk = timclk;

    /* The sample rate must come out exactly, with no rounding. */
    if ((timclk % FS_HZ_INT) != 0U)
    {
        Error_Handler();
    }

    div = timclk / FS_HZ_INT;

    if ((div == 0U) || (div > 65536U))
    {
        Error_Handler();
    }

    dbg_arr = div - 1U;

    htim6.Instance->PSC = 0U;
    htim6.Instance->ARR = dbg_arr;
    htim6.Instance->EGR = TIM_EGR_UG;   /* load the new values immediately */
}

/**
  * @brief  Raised-cosine pulse, normalised so p(0) = 1.
  * @param  x  time in symbol periods (t / T)
  */
static float RaisedCosine(float x)
{
    float sinc;
    float denom;
    float z;

    /* Limit at x = 0 */
    if (fabsf(x) < 1e-7f)
    {
        return 1.0f;
    }

    denom = 1.0f - (2.0f * RC_BETA * x) * (2.0f * RC_BETA * x);

    /* Removable singularity at x = +/- 1 / (2*beta).
     * Limit value is (pi/4) * sinc(1 / (2*beta)).
     */
    if (fabsf(denom) < 1e-5f)
    {
        z = 1.0f / (2.0f * RC_BETA);
        return (PI / 4.0f) * (sinf(PI * z) / (PI * z));
    }

    sinc = sinf(PI * x) / (PI * x);

    return sinc * cosf(PI * RC_BETA * x) / denom;
}

/**
  * @brief  Return the I and Q value for sample n.
  *
  * With NRZ, symbol k occupies samples [k*SPS, (k+1)*SPS).
  * With raised cosine, the pulse for symbol k peaks at sample k*SPS, so the
  * best sampling instant in the receiver is t = k / Rs.
  */
static void Compute_IQ(uint32_t n, float *out_i, float *out_q)
{
#if (PULSE_SHAPE == 0)

    uint32_t k = n / SPS;
    uint8_t  d = symbol_pattern[k];

    *out_i = qpsk_i[d];
    *out_q = qpsk_q[d];

#else

    float acc_i = 0.0f;
    float acc_q = 0.0f;
    float t_sym = (float)n / (float)SPS;   /* time measured in symbols */
    float half  = (float)PATTERN_LEN_SYM * 0.5f;
    uint32_t k;

    for (k = 0U; k < PATTERN_LEN_SYM; k++)
    {
        float tau = t_sym - (float)k;
        float p;
        uint8_t d;

        /* Wrap into [-half, +half].  The buffer is played circularly, so
         * the pulse train must be circular too, or the seam shows up as a
         * glitch once per buffer.
         */
        while (tau >  half) { tau -= (float)PATTERN_LEN_SYM; }
        while (tau < -half) { tau += (float)PATTERN_LEN_SYM; }

        if (fabsf(tau) > (float)RC_HALF_SPAN_SYM)
        {
            continue;
        }

        p = RaisedCosine(tau);
        d = symbol_pattern[k];

        acc_i += qpsk_i[d] * p;
        acc_q += qpsk_q[d] * p;
    }

    *out_i = acc_i;
    *out_q = acc_q;

#endif
}

/**
  * @brief  Fill the DAC buffers.
  *
  * Two passes.  The first finds the true peak of the waveform, the second
  * scales it so the peak lands exactly on DAC_MID +/- DAC_AMP.  This means
  * the amplitude is known and identical whatever the mode or roll-off, and
  * clipping cannot happen by accident.
  *
  * In baseband mode I and Q share one scale factor.  Scaling them
  * separately would introduce I/Q gain imbalance.
  */
static void Build_Buffers(void)
{
    uint32_t n;
    float peak = 0.0f;
    float scale;

    /* ---- Pass 1: find the peak ---- */
    for (n = 0U; n < BUF_LEN; n++)
    {
        float I, Q;

        Compute_IQ(n, &I, &Q);

#if (QPSK_MODE == 0)
        {
            float cyc = (float)n * (FC_HZ / FS_HZ);
            float ph;
            float v;

            cyc -= floorf(cyc);              /* keep the trig argument small */
            ph   = 2.0f * PI * cyc;

            v = I * cosf(ph) - Q * sinf(ph);

            if (fabsf(v) > peak) { peak = fabsf(v); }
        }
#else
        if (fabsf(I) > peak) { peak = fabsf(I); }
        if (fabsf(Q) > peak) { peak = fabsf(Q); }
#endif
    }

    if (peak <= 0.0f)
    {
        Error_Handler();
    }

    dbg_peak  = peak;
    scale     = DAC_AMP / peak;
    dbg_scale = scale;

    /* ---- Pass 2: write the codes ---- */
    for (n = 0U; n < BUF_LEN; n++)
    {
        float I, Q;
        float vi, vq;
        float ci, cq;

        Compute_IQ(n, &I, &Q);

#if (QPSK_MODE == 0)
        {
            float cyc = (float)n * (FC_HZ / FS_HZ);
            float ph;

            cyc -= floorf(cyc);
            ph   = 2.0f * PI * cyc;

            vi = (I * cosf(ph) - Q * sinf(ph)) * scale;
            vq = 0.0f;                        /* PA5 held at midpoint */
        }
#else
        vi = I * scale;
        vq = Q * scale;
#endif

        ci = (float)DAC_MID + vi;
        cq = (float)DAC_MID + vq;

        if (ci <    0.0f) { ci =    0.0f; }
        if (ci > 4095.0f) { ci = 4095.0f; }
        if (cq <    0.0f) { cq =    0.0f; }
        if (cq > 4095.0f) { cq = 4095.0f; }

        dac_buf_i[n] = (uint16_t)(ci + 0.5f);
        dac_buf_q[n] = (uint16_t)(cq + 0.5f);
    }

    dbg_buflen = BUF_LEN;
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

    HAL_Init();

    /* USER CODE BEGIN Init */
    /* USER CODE END Init */

    SystemClock_Config();

    /* USER CODE BEGIN SysInit */
    /* USER CODE END SysInit */

    MX_DMA_Init();
    MX_GPIO_Init();
    MX_DAC1_Init();
    MX_TIM6_Init();

    /* USER CODE BEGIN 2 */

    /* Set the sample-rate divider from the measured bus clock. */
    Configure_TIM6_Runtime();

    /* Build the waveform before anything starts clocking out. */
    Build_Buffers();

    /* Start the DMA channels FIRST, then the timer.
     *
     * Order matters.  Both channels are fed by the same TIM6 TRGO event,
     * so if both DMA channels are armed before the first trigger arrives,
     * I and Q stay sample-aligned for as long as the board runs.  Starting
     * the timer first would let one channel begin a sample ahead of the
     * other, which shows up later as a fixed I/Q skew.
     */
    if (HAL_DAC_Start_DMA(&hdac1,
                          DAC_CHANNEL_1,
                          (uint32_t *)dac_buf_i,
                          BUF_LEN,
                          DAC_ALIGN_12B_R) != HAL_OK)
    {
        Error_Handler();
    }

#if (QPSK_MODE == 1)
    if (HAL_DAC_Start_DMA(&hdac1,
                          DAC_CHANNEL_2,
                          (uint32_t *)dac_buf_q,
                          BUF_LEN,
                          DAC_ALIGN_12B_R) != HAL_OK)
    {
        Error_Handler();
    }
#endif

    if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
    {
        Error_Handler();
    }

    /* USER CODE END 2 */

    /* Infinite loop */
    /* USER CODE BEGIN WHILE */
    while (1)
    {
        /* Everything runs from DMA.  The CPU has nothing left to do. */
        /* USER CODE END WHILE */

        /* USER CODE BEGIN 3 */
    }
    /* USER CODE END 3 */
}

/**
  * @brief  System Clock Configuration
  *
  * Direct 8 MHz HSI, PLL off, FLASH_LATENCY_0, all prescalers /1.
  * SYSCLK = HCLK = PCLK1 = PCLK2 = 8 MHz, so the TIM6 clock is 8 MHz.
  *
  * The PLL is deliberately left off.  RCC_PLLSOURCE_HSI is ambiguous on
  * the F303 between HSI/2 and HSI/PREDIV, and getting it wrong produces a
  * system clock high enough to collapse the DAC output.
  *
  * @retval None
  */
void SystemClock_Config(void)
{
    RCC_OscInitTypeDef RCC_OscInitStruct = {0};
    RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

    RCC_OscInitStruct.OscillatorType      = RCC_OSCILLATORTYPE_HSI;
    RCC_OscInitStruct.HSIState            = RCC_HSI_ON;
    RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
    RCC_OscInitStruct.PLL.PLLState        = RCC_PLL_NONE;

    if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
    {
        Error_Handler();
    }

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
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
    /* USER CODE BEGIN Error_Handler_Debug */
    __disable_irq();
    while (1)
    {
    }
    /* USER CODE END Error_Handler_Debug */
}

#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  */
void assert_failed(uint8_t *file, uint32_t line)
{
    /* USER CODE BEGIN 6 */
    /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
