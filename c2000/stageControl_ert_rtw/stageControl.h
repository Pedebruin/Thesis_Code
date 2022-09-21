/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: stageControl.h
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

#ifndef RTW_HEADER_stageControl_h_
#define RTW_HEADER_stageControl_h_
#ifndef stageControl_COMMON_INCLUDES_
#define stageControl_COMMON_INCLUDES_
#include "rtwtypes.h"
#define __TMWTYPES__                                             /* Inferred types compatibility mode */
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "ext_mode.h"
#include "c2000BoardSupport.h"
#include "MW_f2837xD_includes.h"
#include "MW_c28xIPC.h"
#include "MW_c2000DAC.h"
#endif                                 /* stageControl_COMMON_INCLUDES_ */

#include "stageControl_types.h"
#include <string.h>
#include "MW_target_hardware_resources.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWExtModeInfo
#define rtmGetRTWExtModeInfo(rtm)      ((rtm)->extModeInfo)
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   ((rtm)->Timing.taskTime0)
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                (&(rtm)->Timing.taskTime0)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Strain;                       /* '<Root>/Gain' */
  real_T Position;                     /* '<S3>/Resolution' */
  real_T Reference;                    /* '<Root>/Manual Switch1' */
  real_T ManualSwitch;                 /* '<S1>/Manual Switch' */
  real_T PulseGenerator1;              /* '<S1>/Pulse Generator1' */
  real32_T CastToSingle2;              /* '<Root>/Cast To Single2' */
  real32_T eQEP;                       /* '<S3>/eQEP' */
  real32_T CastToSingle;               /* '<Root>/Cast To Single' */
  real32_T GPIO97Button11;             /* '<Root>/GPIO97 (Button 1)1' */
  real32_T CastToSingle1;              /* '<Root>/Cast To Single1' */
  uint16_T ADC1;                       /* '<Root>/ADC1' */
  uint16_T Acc1;                       /* '<Root>/ADC2' */
} B_stageControl_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteStateSpace_DSTATE[2]; /* '<S1>/Discrete State-Space' */
  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<Root>/Scope' */

  int32_T clockTickCounter;            /* '<Root>/Pulse Generator2' */
  int32_T counter;                     /* '<Root>/Sine Wave1' */
  int32_T clockTickCounter_n;          /* '<S1>/Pulse Generator1' */
  uint32_T Counter1_ClkEphState;       /* '<Root>/Counter1' */
  int16_T PIDcontrol_SubsysRanBC;      /* '<Root>/PID control' */
  uint16_T Counter1_Count;             /* '<Root>/Counter1' */
  boolean_T PIDcontrol_MODE;           /* '<Root>/PID control' */
} DW_stageControl_T;

/* Parameters (default storage) */
struct P_stageControl_T_ {
  struct_PLgAzJHXVR5VPvplI66pVB simulinkSettings;/* Variable: simulinkSettings
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
  uint16_T Counter1_InitialCount;      /* Mask Parameter: Counter1_InitialCount
                                        * Referenced by: '<Root>/Counter1'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 100
                                        * Referenced by: '<S1>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -100
                                        * Referenced by: '<S1>/Saturation'
                                        */
  real_T Output_Y0;                    /* Expression: [0]
                                        * Referenced by: '<S1>/Output'
                                        */
  real_T DiscreteStateSpace_A[3];    /* Computed Parameter: DiscreteStateSpace_A
                                      * Referenced by: '<S1>/Discrete State-Space'
                                      */
  real_T DiscreteStateSpace_B;       /* Computed Parameter: DiscreteStateSpace_B
                                      * Referenced by: '<S1>/Discrete State-Space'
                                      */
  real_T DiscreteStateSpace_C[2];    /* Computed Parameter: DiscreteStateSpace_C
                                      * Referenced by: '<S1>/Discrete State-Space'
                                      */
  real_T DiscreteStateSpace_D;       /* Computed Parameter: DiscreteStateSpace_D
                                      * Referenced by: '<S1>/Discrete State-Space'
                                      */
  real_T DiscreteStateSpace_InitialCondi;/* Expression: 0
                                          * Referenced by: '<S1>/Discrete State-Space'
                                          */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S1>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 100
                                        * Referenced by: '<S1>/Step'
                                        */
  real_T Step1_Y0;                     /* Expression: 100
                                        * Referenced by: '<S1>/Step1'
                                        */
  real_T Step1_YFinal;                 /* Expression: 0
                                        * Referenced by: '<S1>/Step1'
                                        */
  real_T PulseGenerator1_Amp;          /* Expression: 1
                                        * Referenced by: '<S1>/Pulse Generator1'
                                        */
  real_T PulseGenerator1_Period;   /* Computed Parameter: PulseGenerator1_Period
                                    * Referenced by: '<S1>/Pulse Generator1'
                                    */
  real_T PulseGenerator1_Duty;       /* Computed Parameter: PulseGenerator1_Duty
                                      * Referenced by: '<S1>/Pulse Generator1'
                                      */
  real_T PulseGenerator1_PhaseDelay;   /* Expression: 0
                                        * Referenced by: '<S1>/Pulse Generator1'
                                        */
  real_T Constant_Value;               /* Expression: (2^12-1)/2
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain1_Gain;                   /* Expression: 3.3/(2^12-1)
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T PulseGenerator2_PhaseDelay;   /* Expression: 0
                                        * Referenced by: '<Root>/Pulse Generator2'
                                        */
  real_T SineWave1_Bias;               /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave1'
                                        */
  real_T SineWave1_Offset;             /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave1'
                                        */
  real_T Gain_Gain;                    /* Expression: (2^12-2)/100
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T Saturation_UpperSat_d;        /* Expression: 4094
                                        * Referenced by: '<S2>/Saturation'
                                        */
  real_T Saturation_LowerSat_p;        /* Expression: 0
                                        * Referenced by: '<S2>/Saturation'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 0
                                        * Referenced by: '<S2>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: -4095
                                        * Referenced by: '<S2>/Saturation1'
                                        */
  uint16_T ManualSwitch_CurrentSetting;
                              /* Computed Parameter: ManualSwitch_CurrentSetting
                               * Referenced by: '<S1>/Manual Switch'
                               */
  uint16_T ManualSwitch1_CurrentSetting;
                             /* Computed Parameter: ManualSwitch1_CurrentSetting
                              * Referenced by: '<Root>/Manual Switch1'
                              */
  uint16_T ManualSwitch_CurrentSetting_b;
                            /* Computed Parameter: ManualSwitch_CurrentSetting_b
                             * Referenced by: '<Root>/Manual Switch'
                             */
};

/* Real-time Model Data Structure */
struct tag_RTM_stageControl_T {
  const char_T *errorStatus;
  RTWExtModeInfo *extModeInfo;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTick2;
    struct {
      uint16_T TID[3];
    } TaskCounters;

    time_T tFinal;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (default storage) */
extern P_stageControl_T stageControl_P;

/* Block signals (default storage) */
extern B_stageControl_T stageControl_B;

/* Block states (default storage) */
extern DW_stageControl_T stageControl_DW;

/* Model entry point functions */
extern void stageControl_initialize(void);
extern void stageControl_step(void);
extern void stageControl_terminate(void);

/* Real-time Model object */
extern RT_MODEL_stageControl_T *const stageControl_M;
extern volatile boolean_T stopRequested;
extern volatile boolean_T runModel;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S1>/Rate Transition' : Eliminated since input and output rates are identical
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'stageControl'
 * '<S1>'   : 'stageControl/PID control'
 * '<S2>'   : 'stageControl/Stage input'
 * '<S3>'   : 'stageControl/eQEP1'
 */
#endif                                 /* RTW_HEADER_stageControl_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
