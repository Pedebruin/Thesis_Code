/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Filters_private.h
 *
 * Code generated for Simulink model 'Filters'.
 *
 * Model version                  : 3.111
 * Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
 * C/C++ source code generated on : Thu Jun 30 13:51:26 2022
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Texas Instruments->C2000
 * Code generation objectives:
 *    1. RAM efficiency
 *    2. Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Filters_private_h_
#define RTW_HEADER_Filters_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "Filters.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

void isr_int1pie13_task_fcn(void);
void isr_int1pie14_task_fcn(void);
void isr_int1pie15_task_fcn(void);
extern void configureGPIOExtInterrupt(void);
void isr_int1pie13_task_fcn(void);

#endif                                 /* RTW_HEADER_Filters_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
