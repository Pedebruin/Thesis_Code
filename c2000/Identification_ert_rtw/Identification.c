/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Identification.c
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

#include "Identification.h"
#include <math.h>
#include "rtwtypes.h"
#include "Identification_private.h"

/* Block signals (default storage) */
BlockIO_Identification Identification_B;

/* Block states (default storage) */
D_Work_Identification Identification_DWork;

/* Real-time model */
static RT_MODEL_Identification Identification_M_;
RT_MODEL_Identification *const Identification_M = &Identification_M_;

/* Model step function */
void Identification_step(void)
{
  real_T tmp_0;
  int16_T tmp;

  /* Reset subsysRan breadcrumbs */
  srClearBC(Identification_DWork.Subsystem_SubsysRanBC);

  /* Step: '<Root>/Step1' incorporates:
   *  Step: '<Root>/Step'
   */
  tmp_0 = Identification_M->Timing.t[0];
  if (tmp_0 < 80.0) {
    tmp = 0;
  } else {
    tmp = -1;
  }

  /* End of Step: '<Root>/Step1' */

  /* Outputs for Enabled SubSystem: '<Root>/Subsystem' incorporates:
   *  EnablePort: '<S2>/Enable'
   */
  /* Sum: '<Root>/Sum' incorporates:
   *  Step: '<Root>/Step'
   */
  if (!(tmp_0 < 20.0) + tmp > 0) {
    if (!Identification_DWork.Subsystem_MODE) {
      /* InitializeConditions for S-Function (sdspchirp): '<S2>/SampleChirp' */

      /* DSP System Toolbox Chirp (sdspchirp) - '<S2>/SampleChirp' */
      /* Unidirectional Logarithmic  */
      if ((300.0F > 100.0F) && (100.0F > 0.0F) ) {
        Identification_DWork.SampleChirp_MIN_FREQ = 300.0F / 100.0F;
        Identification_DWork.SampleChirp_BETA = 60.0F / log
          (Identification_DWork.SampleChirp_MIN_FREQ);
      }

      Identification_DWork.SampleChirp_PERIOD_THETA =
        Identification_DWork.SampleChirp_BETA * 100.0F * (powf
        (Identification_DWork.SampleChirp_MIN_FREQ , 60.0F / 60.0F) - 1.0F);
      Identification_DWork.SampleChirp_SWEEP_DIRECTION = (300.0F > 100.0F) ? 0 :
        1;
      Identification_DWork.SampleChirp_ACC_PHASE = 0.0F;
      Identification_DWork.SampleChirp_CURRENT_STEP = 0.0F;
      Identification_DWork.Subsystem_MODE = true;
    }

    /* S-Function (sdspchirp): '<S2>/SampleChirp' */
    /* DSP System Toolbox Chirp (sdspchirp) - '<S2>/SampleChirp' */
    /* Unidirectional Logarithmic  */
    {
      real32_T *y = &Identification_B.SampleChirp;
      real32_T instantPhase = 0.0;
      instantPhase = (Identification_DWork.SampleChirp_SWEEP_DIRECTION == 0) ?
        Identification_DWork.SampleChirp_BETA*100.0F*(powf
        (Identification_DWork.SampleChirp_MIN_FREQ,
         Identification_DWork.SampleChirp_CURRENT_STEP / 60.0F) - 1.0F) :
        Identification_DWork.SampleChirp_PERIOD_THETA -
          Identification_DWork.SampleChirp_BETA*100.0F*(powf
          (Identification_DWork.SampleChirp_MIN_FREQ,(60.0F-
          Identification_DWork.SampleChirp_CURRENT_STEP) / 60.0F) - 1.0F);
      instantPhase -= (int_T)instantPhase;
      *y = cosf((real32_T)DSP_TWO_PI * (instantPhase +
                 Identification_DWork.SampleChirp_ACC_PHASE) + 0.0F);
      Identification_DWork.SampleChirp_CURRENT_STEP += 0.0001F;/* Go to next time step */
      if (Identification_DWork.SampleChirp_CURRENT_STEP > (60.0F + FLT_EPSILON))
      {
        Identification_DWork.SampleChirp_CURRENT_STEP =
          Identification_DWork.SampleChirp_CURRENT_STEP - 60.0F;
        Identification_DWork.SampleChirp_ACC_PHASE += instantPhase;
        if (Identification_DWork.SampleChirp_ACC_PHASE > 1.0F) {
          Identification_DWork.SampleChirp_ACC_PHASE -= floor
            (Identification_DWork.SampleChirp_ACC_PHASE);
        }
      }
    }

    srUpdateBC(Identification_DWork.Subsystem_SubsysRanBC);
  } else if (Identification_DWork.Subsystem_MODE) {
    /* Disable for S-Function (sdspchirp): '<S2>/SampleChirp' incorporates:
     *  Outport: '<S2>/Out1'
     */
    Identification_B.SampleChirp = 0.0F;
    Identification_DWork.Subsystem_MODE = false;
  }

  /* End of Sum: '<Root>/Sum' */
  /* End of Outputs for SubSystem: '<Root>/Subsystem' */

  /* Gain: '<Root>/Gain' */
  Identification_B.Input = 100.0F * Identification_B.SampleChirp;

  /* Gain: '<S1>/Gain' */
  Identification_B.Gain = 40.94F * Identification_B.Input;

  /* Saturate: '<S1>/Saturation' */
  if (Identification_B.Gain > 4094.0F) {
    /* Saturate: '<S1>/Saturation' */
    Identification_B.inputPLUS = 4094.0F;
  } else if (Identification_B.Gain < 0.0F) {
    /* Saturate: '<S1>/Saturation' */
    Identification_B.inputPLUS = 0.0F;
  } else {
    /* Saturate: '<S1>/Saturation' */
    Identification_B.inputPLUS = Identification_B.Gain;
  }

  /* End of Saturate: '<S1>/Saturation' */

  /* MATLABSystem: '<S1>/DAC' */
  MW_C2000DACSat(0U, Identification_B.inputPLUS);

  /* Saturate: '<S1>/Saturation1' */
  if (Identification_B.Gain > 0.0F) {
    /* Saturate: '<S1>/Saturation1' */
    Identification_B.inputMIN = 0.0F;
  } else if (Identification_B.Gain < -4095.0F) {
    /* Saturate: '<S1>/Saturation1' */
    Identification_B.inputMIN = -4095.0F;
  } else {
    /* Saturate: '<S1>/Saturation1' */
    Identification_B.inputMIN = Identification_B.Gain;
  }

  /* End of Saturate: '<S1>/Saturation1' */

  /* MATLABSystem: '<S1>/DAC1' incorporates:
   *  Abs: '<S1>/Abs'
   */
  MW_C2000DACSat(1U, fabsf(Identification_B.inputMIN));

  /* S-Function (c280xqep): '<S3>/eQEP' */
  {
    Identification_B.eQEP_o1 = EQep1Regs.QPOSCNT;/*eQEP Position Counter*/
    Identification_B.eQEP_o2 = EQep1Regs.QEPSTS.bit.PCEF;
    /* Position counter error flag. This bit is not sticky and it is updated for every index event*/

    /* V1.1 Added UPEVNT (bit 7) This reflects changes made as of F280x Rev A devices==>Currently, our target board "TMS320F2808eZdsp" is Rev 0.
     *         if(EQep1Regs.QEPSTS.bit.UPEVNT==1){
     */
    if (EQep1Regs.QEPSTS.bit.COEF ==0 && EQep1Regs.QEPSTS.bit.CDEF ==0)
      Identification_B.eQEP_o3 = EQep1Regs.QCPRD;
                   /* eQEP Capture Period (QCPRD) Register : No Capture overflow
                      else
                      Identification_B.eQEP_o3 = 65535;      eQEP Capture Period (QCPRD) Register : Capture overflow, saturate the result
                      EQep1Regs.QEPSTS.bit.UPEVNT==1;
                      }*/
    if (EQep1Regs.QEPSTS.bit.COEF ==1)
      EQep1Regs.QEPSTS.bit.COEF = 1;
    if (EQep1Regs.QEPSTS.bit.CDEF ==1)
      EQep1Regs.QEPSTS.bit.CDEF = 1;
  }

  /* Gain: '<S3>/Resolution' incorporates:
   *  Constant: '<S3>/Constant'
   *  Sum: '<S3>/Sum'
   */
  Identification_B.Position = (Identification_B.eQEP_o1 - 10000.0) * 1.0E-8;

  {                                    /* Sample time: [0.0s, 0.0s] */
    extmodeErrorCode_T errorCode = EXTMODE_SUCCESS;
    extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T)
      ((Identification_M->Timing.clockTick0 * 1) + 0)
      ;

    /* Trigger External Mode event */
    errorCode = extmodeEvent(0,currentTime);
    if (errorCode != EXTMODE_SUCCESS) {
      /* Code to handle External Mode event errors
         may be added here */
    }
  }

  {                                    /* Sample time: [0.0001s, 0.0s] */
    extmodeErrorCode_T errorCode = EXTMODE_SUCCESS;
    extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T)
      ((Identification_M->Timing.clockTick1 * 1) + 0)
      ;

    /* Trigger External Mode event */
    errorCode = extmodeEvent(1,currentTime);
    if (errorCode != EXTMODE_SUCCESS) {
      /* Code to handle External Mode event errors
         may be added here */
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   */
  Identification_M->Timing.t[0] =
    ((time_T)(++Identification_M->Timing.clockTick0)) *
    Identification_M->Timing.stepSize0;

  {
    /* Update absolute timer for sample time: [0.0001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.0001, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     */
    Identification_M->Timing.clockTick1++;
  }
}

/* Model initialize function */
void Identification_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)Identification_M, 0,
                sizeof(RT_MODEL_Identification));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Identification_M->solverInfo,
                          &Identification_M->Timing.simTimeStep);
    rtsiSetTPtr(&Identification_M->solverInfo, &rtmGetTPtr(Identification_M));
    rtsiSetStepSizePtr(&Identification_M->solverInfo,
                       &Identification_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&Identification_M->solverInfo, (&rtmGetErrorStatus
      (Identification_M)));
    rtsiSetRTModelPtr(&Identification_M->solverInfo, Identification_M);
  }

  rtsiSetSimTimeStep(&Identification_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&Identification_M->solverInfo,"FixedStepDiscrete");
  rtmSetTPtr(Identification_M, &Identification_M->Timing.tArray[0]);
  rtmSetTFinal(Identification_M, -1);
  Identification_M->Timing.stepSize0 = 0.0001;

  /* External mode info */
  Identification_M->Sizes.checksums[0] = (3789599025U);
  Identification_M->Sizes.checksums[1] = (1398328255U);
  Identification_M->Sizes.checksums[2] = (719049494U);
  Identification_M->Sizes.checksums[3] = (2920486555U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[4];
    Identification_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    systemRan[3] = (sysRanDType *)&Identification_DWork.Subsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(Identification_M->extModeInfo,
      &Identification_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(Identification_M->extModeInfo,
                        Identification_M->Sizes.checksums);
    rteiSetTPtr(Identification_M->extModeInfo, rtmGetTPtr(Identification_M));
  }

  /* block I/O */
  (void) memset(((void *) &Identification_B), 0,
                sizeof(BlockIO_Identification));

  /* states (dwork) */
  (void) memset((void *)&Identification_DWork, 0,
                sizeof(D_Work_Identification));

  /* Start for MATLABSystem: '<S1>/DAC' */
  MW_ConfigureDACA();

  /* Start for MATLABSystem: '<S1>/DAC1' */
  MW_ConfigureDACB();

  /* Start for S-Function (c280xqep): '<S3>/eQEP' */
  config_QEP_eQEP1(20000U, 10000U, 0, 0, 0, 0, 136, 32768, 32887, 0);

  /* SystemInitialize for Enabled SubSystem: '<Root>/Subsystem' */

  /* InitializeConditions for S-Function (sdspchirp): '<S2>/SampleChirp' */

  /* DSP System Toolbox Chirp (sdspchirp) - '<S2>/SampleChirp' */
  /* Unidirectional Logarithmic  */
  if ((300.0F > 100.0F) && (100.0F > 0.0F) ) {
    Identification_DWork.SampleChirp_MIN_FREQ = 300.0F / 100.0F;
    Identification_DWork.SampleChirp_BETA = 60.0F / log
      (Identification_DWork.SampleChirp_MIN_FREQ);
  }

  Identification_DWork.SampleChirp_PERIOD_THETA =
    Identification_DWork.SampleChirp_BETA * 100.0F * (powf
    (Identification_DWork.SampleChirp_MIN_FREQ , 60.0F / 60.0F) - 1.0F);
  Identification_DWork.SampleChirp_SWEEP_DIRECTION = (300.0F > 100.0F) ? 0 : 1;
  Identification_DWork.SampleChirp_ACC_PHASE = 0.0F;
  Identification_DWork.SampleChirp_CURRENT_STEP = 0.0F;

  /* End of SystemInitialize for SubSystem: '<Root>/Subsystem' */
}

/* Model terminate function */
void Identification_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
