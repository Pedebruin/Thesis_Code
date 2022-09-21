/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Filters.h
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

#ifndef RTW_HEADER_Filters_h_
#define RTW_HEADER_Filters_h_
#ifndef Filters_COMMON_INCLUDES_
#define Filters_COMMON_INCLUDES_
#include "rtwtypes.h"
#define __TMWTYPES__                                             /* Inferred types compatibility mode */
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "MW_c28xIPC.h"
#include "c2000BoardSupport.h"
#include "MW_f2837xD_includes.h"
#endif                                 /* Filters_COMMON_INCLUDES_ */

#include "Filters_types.h"
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

extern void configureCPU2Peripherals(uint32_T gpioNumber, uint32_T gpGRegValA,
  uint32_T gpRegValA);

/* Block signals (default storage) */
typedef struct {
  real_T PulseGenerator1;              /* '<Root>/Pulse Generator1' */
  real32_T IPCReceive_o1;              /* '<S3>/IPC Receive' */
  real32_T IPCReceive_o1_l;            /* '<S1>/IPC Receive' */
  uint16_T IPCReceive_o2;              /* '<S3>/IPC Receive' */
  uint16_T IPCReceive_o1_j;            /* '<S2>/IPC Receive' */
  uint16_T IPCReceive_o2_l;            /* '<S2>/IPC Receive' */
  uint16_T IPCReceive_o2_d;            /* '<S1>/IPC Receive' */
} B_Filters_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  int32_T clockTickCounter;            /* '<Root>/Pulse Generator1' */
  int16_T FunctionCallSubsystem2_SubsysRa;/* '<Root>/Function-Call Subsystem2' */
  int16_T FunctionCallSubsystem1_SubsysRa;/* '<Root>/Function-Call Subsystem1' */
  int16_T FunctionCallSubsystem_SubsysRan;/* '<Root>/Function-Call Subsystem' */
} DW_Filters_T;

/* Parameters (default storage) */
struct P_Filters_T_ {
  real_T PulseGenerator1_Amp;          /* Expression: 1
                                        * Referenced by: '<Root>/Pulse Generator1'
                                        */
  real_T PulseGenerator1_Period;   /* Computed Parameter: PulseGenerator1_Period
                                    * Referenced by: '<Root>/Pulse Generator1'
                                    */
  real_T PulseGenerator1_Duty;       /* Computed Parameter: PulseGenerator1_Duty
                                      * Referenced by: '<Root>/Pulse Generator1'
                                      */
  real_T PulseGenerator1_PhaseDelay;   /* Expression: 0
                                        * Referenced by: '<Root>/Pulse Generator1'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Filters_T {
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
    uint8_T rtmL2HLastBufWr;
    uint32_T rtmL2HDbBufClockTick[2];
    uint32_T clockTick1;
    uint32_T clockTick2;
    uint32_T clockTick3;
    time_T tFinal;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (default storage) */
extern P_Filters_T Filters_P;

/* Block signals (default storage) */
extern B_Filters_T Filters_B;

/* Block states (default storage) */
extern DW_Filters_T Filters_DW;

/* Model entry point functions */
extern void Filters_initialize(void);
extern void Filters_step(void);
extern void Filters_terminate(void);
extern volatile boolean_T runModel;

/* Real-time Model object */
extern RT_MODEL_Filters_T *const Filters_M;
extern volatile boolean_T stopRequested;
extern volatile boolean_T runModel;

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
 * '<Root>' : 'Filters'
 * '<S1>'   : 'Filters/Function-Call Subsystem'
 * '<S2>'   : 'Filters/Function-Call Subsystem1'
 * '<S3>'   : 'Filters/Function-Call Subsystem2'
 */
#endif                                 /* RTW_HEADER_Filters_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
