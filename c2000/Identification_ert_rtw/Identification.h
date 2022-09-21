/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Identification.h
 *
 * Code generated for Simulink model 'Identification'.
 *
 * Model version                  : 4.0
 * Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
 * C/C++ source code generated on : Mon May 23 15:26:25 2022
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Texas Instruments->C2000
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Identification_h_
#define RTW_HEADER_Identification_h_
#ifndef Identification_COMMON_INCLUDES_
#define Identification_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "ext_mode.h"
#include "c2000BoardSupport.h"
#include "MW_f2837xD_includes.h"
#include "MW_c2000DAC.h"
#endif                                 /* Identification_COMMON_INCLUDES_ */

#include "Identification_types.h"
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
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Position;                     /* '<S3>/Resolution' */
  real32_T Input;                      /* '<Root>/Gain' */
  real32_T Gain;                       /* '<S1>/Gain' */
  real32_T inputPLUS;                  /* '<S1>/Saturation' */
  real32_T inputMIN;                   /* '<S1>/Saturation1' */
  real32_T eQEP_o1;                    /* '<S3>/eQEP' */
  real32_T eQEP_o2;                    /* '<S3>/eQEP' */
  real32_T eQEP_o3;                    /* '<S3>/eQEP' */
  real32_T SampleChirp;                /* '<S2>/SampleChirp' */
} BlockIO_Identification;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real32_T SampleChirp_CURRENT_STEP;   /* '<S2>/SampleChirp' */
  real32_T SampleChirp_ACC_PHASE;      /* '<S2>/SampleChirp' */
  real32_T SampleChirp_BETA;           /* '<S2>/SampleChirp' */
  real32_T SampleChirp_MIN_FREQ;       /* '<S2>/SampleChirp' */
  real32_T SampleChirp_PERIOD_THETA;   /* '<S2>/SampleChirp' */
  int16_T Subsystem_SubsysRanBC;       /* '<Root>/Subsystem' */
  boolean_T SampleChirp_SWEEP_DIRECTION;/* '<S2>/SampleChirp' */
  boolean_T Subsystem_MODE;            /* '<Root>/Subsystem' */
} D_Work_Identification;

/* Real-time Model Data Structure */
struct tag_RTM_Identification {
  const char_T *errorStatus;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;

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
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
extern BlockIO_Identification Identification_B;

/* Block states (default storage) */
extern D_Work_Identification Identification_DWork;

/* Model entry point functions */
extern void Identification_initialize(void);
extern void Identification_step(void);
extern void Identification_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Identification *const Identification_M;
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
 * '<Root>' : 'Identification'
 * '<S1>'   : 'Identification/Stage input'
 * '<S2>'   : 'Identification/Subsystem'
 * '<S3>'   : 'Identification/eQEP'
 */
#endif                                 /* RTW_HEADER_Identification_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
