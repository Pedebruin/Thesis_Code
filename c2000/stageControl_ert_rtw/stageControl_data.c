/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: stageControl_data.c
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

#include "stageControl.h"

/* Block parameters (default storage) */
P_stageControl_T stageControl_P = {
  /* Variable: simulinkSettings
   * Referenced by:
   *   '<Root>/Constant1'
   *   '<Root>/Constant3'
   *   '<Root>/Pulse Generator2'
   *   '<Root>/Gain'
   *   '<Root>/Sine Wave1'
   *   '<S1>/Step'
   *   '<S1>/Step1'
   *   '<S3>/Constant'
   *   '<S3>/Resolution'
   */
  {
    0.0001,
    0.0001,
    3.95e-8,
    30000.0,
    -1.3636e-6,
    -83.8084733536745,
    269.15348039269168,
    100.0,
    100.0,
    300.0,
    60.0,
    1000.0,
    0.0001,
    0.0003,
    10.0,
    10.0,
    10.0,
    0.05
  },

  /* Mask Parameter: Counter1_InitialCount
   * Referenced by: '<Root>/Counter1'
   */
  0U,

  /* Expression: 100
   * Referenced by: '<S1>/Saturation'
   */
  100.0,

  /* Expression: -100
   * Referenced by: '<S1>/Saturation'
   */
  -100.0,

  /* Expression: [0]
   * Referenced by: '<S1>/Output'
   */
  0.0,

  /* Computed Parameter: DiscreteStateSpace_A
   * Referenced by: '<S1>/Discrete State-Space'
   */
  { 1.1518358019806489, -0.30367160396129789, 0.5 },

  /* Computed Parameter: DiscreteStateSpace_B
   * Referenced by: '<S1>/Discrete State-Space'
   */
  8192.0,

  /* Computed Parameter: DiscreteStateSpace_C
   * Referenced by: '<S1>/Discrete State-Space'
   */
  { -3359.4942193134066, 6719.0939892419046 },

  /* Computed Parameter: DiscreteStateSpace_D
   * Referenced by: '<S1>/Discrete State-Space'
   */
  3.6506771322837226E+7,

  /* Expression: 0
   * Referenced by: '<S1>/Discrete State-Space'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/Step'
   */
  0.0,

  /* Expression: 100
   * Referenced by: '<S1>/Step'
   */
  100.0,

  /* Expression: 100
   * Referenced by: '<S1>/Step1'
   */
  100.0,

  /* Expression: 0
   * Referenced by: '<S1>/Step1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S1>/Pulse Generator1'
   */
  1.0,

  /* Computed Parameter: PulseGenerator1_Period
   * Referenced by: '<S1>/Pulse Generator1'
   */
  2.0,

  /* Computed Parameter: PulseGenerator1_Duty
   * Referenced by: '<S1>/Pulse Generator1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S1>/Pulse Generator1'
   */
  0.0,

  /* Expression: (2^12-1)/2
   * Referenced by: '<Root>/Constant'
   */
  2047.5,

  /* Expression: 3.3/(2^12-1)
   * Referenced by: '<Root>/Gain1'
   */
  0.00080586080586080586,

  /* Expression: 0
   * Referenced by: '<Root>/Pulse Generator2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave1'
   */
  0.0,

  /* Expression: (2^12-2)/100
   * Referenced by: '<S2>/Gain'
   */
  40.94,

  /* Expression: 4094
   * Referenced by: '<S2>/Saturation'
   */
  4094.0,

  /* Expression: 0
   * Referenced by: '<S2>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S2>/Saturation1'
   */
  0.0,

  /* Expression: -4095
   * Referenced by: '<S2>/Saturation1'
   */
  -4095.0,

  /* Computed Parameter: ManualSwitch_CurrentSetting
   * Referenced by: '<S1>/Manual Switch'
   */
  1U,

  /* Computed Parameter: ManualSwitch1_CurrentSetting
   * Referenced by: '<Root>/Manual Switch1'
   */
  0U,

  /* Computed Parameter: ManualSwitch_CurrentSetting_b
   * Referenced by: '<Root>/Manual Switch'
   */
  1U
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
