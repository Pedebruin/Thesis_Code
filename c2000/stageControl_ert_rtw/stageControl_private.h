/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: stageControl_private.h
 *
 * Code generated for Simulink model 'stageControl'.
 *
 * Model version                  : 3.112
 * Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
 * C/C++ source code generated on : Tue Sep  6 11:00:21 2022
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Texas Instruments->C2000
 * Code generation objectives:
 *    1. RAM efficiency
 *    2. Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_stageControl_private_h_
#define RTW_HEADER_stageControl_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

void config_QEP_eQEP1(uint32_T pcmaximumvalue, uint32_T pcInitialvalue, uint32_T
                      unittimerperiod, uint32_T comparevalue, uint16_T
                      watchdogtimer, uint16_T qdecctl, uint16_T qepctl, uint16_T
                      qposctl, uint16_T qcapctl, uint16_T qeint);
void InitAdcA (void);
void config_ADCA_SOC0 (void);
void InitAdcD (void);
void config_ADCD_SOC0 (void);
extern uint16_T MW_adcAInitFlag;
extern uint16_T MW_adcDInitFlag;
extern uint32_T MWDSP_EPH_R_R(real32_T evt, uint32_T *sta);

#endif                                 /* RTW_HEADER_stageControl_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
