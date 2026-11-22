/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    dma.c
  * @brief   This file provides code for the configuration
  *          of all the requested memory to memory DMA transfers.
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "dma.h"

/* USER CODE BEGIN 0 */

/* ---------------------------------------------------------------------------
 *  NOTE FOR THE THD BUILD
 *
 *  The two memory-to-memory channels that CubeMX previously generated here
 *  (hdma_memtomem_dma1_channel1 and hdma_memtomem_dma2_channel1) have been
 *  removed. Nothing in this project used them. They only enabled the DMA2
 *  clock, which adds a little supply current and switching noise for no
 *  benefit during a distortion measurement.
 *
 *  Remove them in CubeMX as well:
 *      System Core > DMA > DMA1 / DMA2 tab > select the MemToMem entry > Delete
 *  Otherwise the next code generation will put them back.
 *
 *  The DAC channel itself (hdma_dac1_ch1 on DMA1_Channel3) is NOT configured
 *  in this file. CubeMX places it in HAL_DAC_MspInit() inside dac.c. That is
 *  correct and expected. Check there that the settings are:
 *
 *      Instance             = DMA1_Channel3
 *      Direction            = DMA_MEMORY_TO_PERIPH
 *      PeriphInc            = DMA_PINC_DISABLE
 *      MemInc               = DMA_MINC_ENABLE
 *      PeriphDataAlignment  = DMA_PDATAALIGN_HALFWORD
 *      MemDataAlignment     = DMA_MDATAALIGN_HALFWORD
 *      Mode                 = DMA_CIRCULAR
 *      Priority             = DMA_PRIORITY_HIGH
 *
 *  Circular mode and the two half-word alignments are the ones that will
 *  silently ruin the output if they are wrong.
 * ------------------------------------------------------------------------- */

/* USER CODE END 0 */

/*----------------------------------------------------------------------------*/
/* Configure DMA                                                              */
/*----------------------------------------------------------------------------*/

/* USER CODE BEGIN 1 */

/* USER CODE END 1 */

/**
  * Enable DMA controller clock
  */
void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Channel3_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Channel3_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Channel3_IRQn);

  /* The transfer-complete and half-transfer interrupts on this channel are
   * masked in main.c immediately after HAL_DAC_Start_DMA(). They are not
   * disabled here, because HAL_DAC_Start_DMA() re-enables them internally
   * every time it is called. See Start_Output() in main.c.
   *
   * Reason for masking them: the 1000-sample buffer wraps every 10 ms, so
   * those interrupts fire at 100 Hz. Harmonics of 100 Hz fall on 1 kHz,
   * 2 kHz and 3 kHz, which are exactly the fundamental and the first
   * harmonics being measured. Keeping the core asleep removes that
   * contribution from the spectrum.
   *
   * The NVIC entry is left enabled so that DMA transfer-error interrupts
   * are still delivered. */
}

/* USER CODE BEGIN 2 */

/* USER CODE END 2 */
