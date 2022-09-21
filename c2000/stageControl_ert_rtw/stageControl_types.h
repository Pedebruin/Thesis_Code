/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: stageControl_types.h
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

#ifndef RTW_HEADER_stageControl_types_h_
#define RTW_HEADER_stageControl_types_h_
#include "rtwtypes.h"

/* Model Code Variants */
#ifndef DEFINED_TYPEDEF_FOR_struct_PLgAzJHXVR5VPvplI66pVB_
#define DEFINED_TYPEDEF_FOR_struct_PLgAzJHXVR5VPvplI66pVB_

typedef struct {
  real_T measurementSampleTime;
  real_T PIDSampleTime;
  real_T laserResolution;
  real_T maxEncoderPosition;
  real_T strainGain;
  real_T accSensitivity;
  real_T matchingGain;
  real_T inputGain;
  real_T chirpInitFreq;
  real_T chirpTargetFreq;
  real_T chirpTargetTime;
  real_T controlBandwidth;
  real_T referenceSampleTime;
  real_T referenceAmplitude;
  real_T referenceFrequency;
  real_T impulseTime;
  real_T impulseSamples;
  real_T inputForce;
} struct_PLgAzJHXVR5VPvplI66pVB;

#endif

#ifndef struct_tag_mJwWZQsvJAXvgRMjJItCbG
#define struct_tag_mJwWZQsvJAXvgRMjJItCbG

struct tag_mJwWZQsvJAXvgRMjJItCbG
{
  int32_T isInitialized;
};

#endif                                 /* struct_tag_mJwWZQsvJAXvgRMjJItCbG */

#ifndef typedef_codertarget_tic2000_blocks_DA_T
#define typedef_codertarget_tic2000_blocks_DA_T

typedef struct tag_mJwWZQsvJAXvgRMjJItCbG codertarget_tic2000_blocks_DA_T;

#endif                             /* typedef_codertarget_tic2000_blocks_DA_T */

/* Parameters (default storage) */
typedef struct P_stageControl_T_ P_stageControl_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_stageControl_T RT_MODEL_stageControl_T;

#endif                                 /* RTW_HEADER_stageControl_types_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
