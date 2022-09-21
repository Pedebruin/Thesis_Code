/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: stageControl.c
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
#include "rtwtypes.h"
#include "stageControl_private.h"
#include <math.h>

/* Block signals (default storage) */
B_stageControl_T stageControl_B;

/* Block states (default storage) */
DW_stageControl_T stageControl_DW;

/* Real-time model */
static RT_MODEL_stageControl_T stageControl_M_;
RT_MODEL_stageControl_T *const stageControl_M = &stageControl_M_;
static void rate_scheduler(void);

#ifndef __TMS320C28XX_CLA__

uint16_T MW_adcAInitFlag = 0;

#endif

#ifndef __TMS320C28XX_CLA__

uint16_T MW_adcDInitFlag = 0;

#endif

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (stageControl_M->Timing.TaskCounters.TID[1])++;
  if ((stageControl_M->Timing.TaskCounters.TID[1]) > 99) {/* Sample time: [0.01s, 0.0s] */
    stageControl_M->Timing.TaskCounters.TID[1] = 0;
  }

  (stageControl_M->Timing.TaskCounters.TID[2])++;
  if ((stageControl_M->Timing.TaskCounters.TID[2]) > 4999) {/* Sample time: [0.5s, 0.0s] */
    stageControl_M->Timing.TaskCounters.TID[2] = 0;
  }
}

uint32_T MWDSP_EPH_R_R(real32_T evt, uint32_T *sta)
{
  uint32_T previousState;
  uint32_T retVal;
  int16_T curState;
  int16_T lastzcevent;
  int16_T newState;
  int16_T newStateR;

  /* S-Function (sdspcount2): '<Root>/Counter1' */
  /* Detect rising edge events */
  previousState = *sta;
  retVal = 0UL;
  lastzcevent = 0;
  newState = 5;
  newStateR = 5;
  if (evt > 0.0F) {
    curState = 2;
  } else {
    curState = !(evt < 0.0F);
  }

  if (*sta == 5UL) {
    newStateR = curState;
  } else if ((uint32_T)curState != *sta) {
    if (*sta == 3UL) {
      if ((uint16_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 2;
        previousState = 1UL;
      }
    }

    if (previousState == 4UL) {
      if ((uint16_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 3;
        previousState = 1UL;
      }
    }

    if ((previousState == 1UL) && ((uint16_T)curState == 2U)) {
      retVal = 2UL;
    }

    if (previousState == 0UL) {
      retVal = 2UL;
    }

    if ((uint16_T)retVal == (uint16_T)lastzcevent) {
      retVal = 0UL;
    }

    if (((uint16_T)curState == 1U) && ((uint16_T)retVal == 2U)) {
      newState = 3;
    } else {
      newState = curState;
    }
  }

  if ((uint16_T)newStateR != 5U) {
    *sta = (uint32_T)newStateR;
    retVal = 0UL;
  }

  if ((uint16_T)newState != 5U) {
    *sta = (uint32_T)newState;
  }

  /* End of S-Function (sdspcount2): '<Root>/Counter1' */
  return retVal;
}

/* Model step function */
void stageControl_step(void)
{
  /* local block i/o variables */
  real_T rtb_input1;
  real_T rtb_DiscreteStateSpace;
  real_T rtb_PulseGenerator2;
  real_T rtb_inputMIN;
  boolean_T rtb_Enable;

  /* Reset subsysRan breadcrumbs */
  srClearBC(stageControl_DW.PIDcontrol_SubsysRanBC);

  /* S-Function (c2802xadc): '<Root>/ADC1' */
  {
    /*  Internal Reference Voltage : Fixed scale 0 to 3.3 V range.  */
    /*  External Reference Voltage : Allowable ranges of VREFHI(ADCINA0) = 3.3 and VREFLO(tied to ground) = 0  */
    AdcaRegs.ADCSOCFRC1.bit.SOC0 = 1;

    /* Wait for the period of Sampling window and EOC result to be latched after trigger */
#ifndef __TMS320C28XX_CLA__

    asm(" RPT #75|| NOP");

#endif

#ifdef __TMS320C28XX_CLA__

    float wait_index;
    for (wait_index= 12; wait_index > 0; wait_index--)
      __mnop();

#endif

    stageControl_B.ADC1 = (AdcaResultRegs.ADCRESULT0);
  }

  /* Gain: '<Root>/Gain' incorporates:
   *  Constant: '<Root>/Constant'
   *  Gain: '<Root>/Gain1'
   *  Sum: '<Root>/Sum'
   */
  stageControl_B.Strain = ((real_T)stageControl_B.ADC1 -
    stageControl_P.Constant_Value) * stageControl_P.Gain1_Gain *
    stageControl_P.simulinkSettings.strainGain;

  /* DataTypeConversion: '<Root>/Cast To Single2' */
  stageControl_B.CastToSingle2 = (real32_T)stageControl_B.Strain;

  /* S-Function (c28xipc_tx): '<Root>/IPC Transmit1' */
  MW_IPC_Transmit(CHANNEL2, (uint32_t *)&stageControl_B.CastToSingle2, 1, 8, 0,
                  0);

  /* S-Function (c2802xadc): '<Root>/ADC2' */
  {
    /*  Internal Reference Voltage : Fixed scale 0 to 3.3 V range.  */
    /*  External Reference Voltage : Allowable ranges of VREFHI(ADCINA0) = 3.3 and VREFLO(tied to ground) = 0  */
    AdcdRegs.ADCSOCFRC1.bit.SOC0 = 1;

    /* Wait for the period of Sampling window and EOC result to be latched after trigger */
#ifndef __TMS320C28XX_CLA__

    asm(" RPT #170|| NOP");

#endif

#ifdef __TMS320C28XX_CLA__

    float wait_index;
    for (wait_index= 28; wait_index > 0; wait_index--)
      __mnop();

#endif

    stageControl_B.Acc1 = (AdcdResultRegs.ADCRESULT0);
  }

  /* S-Function (c280xqep): '<S3>/eQEP' */
  {
    stageControl_B.eQEP = EQep1Regs.QPOSCNT;/*eQEP Position Counter*/
  }

  /* Gain: '<S3>/Resolution' incorporates:
   *  Constant: '<S3>/Constant'
   *  Sum: '<S3>/Sum'
   */
  stageControl_B.Position = (stageControl_B.eQEP -
    stageControl_P.simulinkSettings.maxEncoderPosition / 2.0) *
    stageControl_P.simulinkSettings.laserResolution;

  /* DataTypeConversion: '<Root>/Cast To Single' */
  stageControl_B.CastToSingle = (real32_T)stageControl_B.Position;

  /* S-Function (c28xipc_tx): '<Root>/IPC Transmit' */
  MW_IPC_Transmit(CHANNEL0, (uint32_t *)&stageControl_B.CastToSingle, 1, 8, 0, 0);
  if (stageControl_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (c280xgpio_di): '<Root>/GPIO97 (Button 1)1' */
    {
      stageControl_B.GPIO97Button11 = GpioDataRegs.GPBDAT.bit.GPIO56;
    }

    /* S-Function (sdspcount2): '<Root>/Counter1' */
    if (MWDSP_EPH_R_R(stageControl_B.GPIO97Button11,
                      &stageControl_DW.Counter1_ClkEphState) != 0UL) {
      if (stageControl_DW.Counter1_Count < 1U) {
        stageControl_DW.Counter1_Count = 1U;
      } else {
        stageControl_DW.Counter1_Count = 0U;
      }
    }

    /* DataTypeConversion: '<Root>/Cast To Boolean1' incorporates:
     *  S-Function (sdspcount2): '<Root>/Counter1'
     */
    rtb_Enable = (stageControl_DW.Counter1_Count != 0U);
  }

  /* DiscretePulseGenerator: '<Root>/Pulse Generator2' incorporates:
   *  Sin: '<Root>/Sine Wave1'
   */
  rtb_PulseGenerator2 = (stageControl_DW.clockTickCounter < 1.0 /
    stageControl_P.simulinkSettings.referenceSampleTime * 2.5) &&
    (stageControl_DW.clockTickCounter >= 0L) ?
    stageControl_P.simulinkSettings.referenceAmplitude : 0.0;
  rtb_inputMIN = 1.0 / stageControl_P.simulinkSettings.referenceSampleTime;
  if (stageControl_DW.clockTickCounter >= rtb_inputMIN * 5.0 - 1.0) {
    stageControl_DW.clockTickCounter = 0L;
  } else {
    stageControl_DW.clockTickCounter++;
  }

  /* End of DiscretePulseGenerator: '<Root>/Pulse Generator2' */

  /* ManualSwitch: '<Root>/Manual Switch1' */
  if (stageControl_P.ManualSwitch1_CurrentSetting == 1U) {
    /* ManualSwitch: '<Root>/Manual Switch1' incorporates:
     *  Constant: '<Root>/Constant3'
     *  Sum: '<Root>/Sum3'
     */
    stageControl_B.Reference = rtb_PulseGenerator2 -
      stageControl_P.simulinkSettings.referenceAmplitude / 2.0;
  } else {
    /* ManualSwitch: '<Root>/Manual Switch1' incorporates:
     *  Sin: '<Root>/Sine Wave1'
     */
    stageControl_B.Reference = sin(((real_T)stageControl_DW.counter +
      stageControl_P.SineWave1_Offset) * 2.0 * 3.1415926535897931 /
      (rtb_inputMIN / stageControl_P.simulinkSettings.referenceFrequency)) *
      stageControl_P.simulinkSettings.referenceAmplitude +
      stageControl_P.SineWave1_Bias;
  }

  /* End of ManualSwitch: '<Root>/Manual Switch1' */

  /* Outputs for Enabled SubSystem: '<Root>/PID control' incorporates:
   *  EnablePort: '<S1>/Enable'
   */
  if (stageControl_M->Timing.TaskCounters.TID[1] == 0) {
    if (rtb_Enable) {
      stageControl_DW.PIDcontrol_MODE = true;
    } else if (stageControl_DW.PIDcontrol_MODE) {
      /* Disable for ManualSwitch: '<S1>/Manual Switch' incorporates:
       *  Outport: '<S1>/Output'
       */
      stageControl_B.ManualSwitch = stageControl_P.Output_Y0;
      stageControl_DW.PIDcontrol_MODE = false;
    }
  }

  if (stageControl_DW.PIDcontrol_MODE) {
    /* Sum: '<S1>/Sum1' */
    rtb_input1 = stageControl_B.Position - stageControl_B.Reference;

    /* DiscreteStateSpace: '<S1>/Discrete State-Space' */
    {
      rtb_DiscreteStateSpace = (stageControl_P.DiscreteStateSpace_C[0])*
        stageControl_DW.DiscreteStateSpace_DSTATE[0]
        + (stageControl_P.DiscreteStateSpace_C[1])*
        stageControl_DW.DiscreteStateSpace_DSTATE[1];
      rtb_DiscreteStateSpace += stageControl_P.DiscreteStateSpace_D*rtb_input1;
    }

    /* ManualSwitch: '<S1>/Manual Switch' incorporates:
     *  Saturate: '<S1>/Saturation'
     */
    if (stageControl_P.ManualSwitch_CurrentSetting == 1U) {
      /* Step: '<S1>/Step' incorporates:
       *  Step: '<S1>/Step1'
       */
      rtb_inputMIN = stageControl_M->Timing.taskTime0;
      if (rtb_inputMIN < stageControl_P.simulinkSettings.impulseTime) {
        rtb_PulseGenerator2 = stageControl_P.Step_Y0;
      } else {
        rtb_PulseGenerator2 = stageControl_P.Step_YFinal;
      }

      /* End of Step: '<S1>/Step' */

      /* Step: '<S1>/Step1' */
      if (rtb_inputMIN < stageControl_P.simulinkSettings.impulseSamples *
          stageControl_P.simulinkSettings.measurementSampleTime +
          stageControl_P.simulinkSettings.impulseTime) {
        rtb_inputMIN = stageControl_P.Step1_Y0;
      } else {
        rtb_inputMIN = stageControl_P.Step1_YFinal;
      }

      /* ManualSwitch: '<S1>/Manual Switch' incorporates:
       *  Product: '<S1>/Product'
       */
      stageControl_B.ManualSwitch = rtb_PulseGenerator2 * rtb_inputMIN;
    } else if (rtb_DiscreteStateSpace > stageControl_P.Saturation_UpperSat) {
      /* Saturate: '<S1>/Saturation' incorporates:
       *  ManualSwitch: '<S1>/Manual Switch'
       */
      stageControl_B.ManualSwitch = stageControl_P.Saturation_UpperSat;
    } else if (rtb_DiscreteStateSpace < stageControl_P.Saturation_LowerSat) {
      /* Saturate: '<S1>/Saturation' incorporates:
       *  ManualSwitch: '<S1>/Manual Switch'
       */
      stageControl_B.ManualSwitch = stageControl_P.Saturation_LowerSat;
    } else {
      /* ManualSwitch: '<S1>/Manual Switch' incorporates:
       *  Saturate: '<S1>/Saturation'
       */
      stageControl_B.ManualSwitch = rtb_DiscreteStateSpace;
    }

    /* End of ManualSwitch: '<S1>/Manual Switch' */
    if (stageControl_M->Timing.TaskCounters.TID[2] == 0) {
      /* DiscretePulseGenerator: '<S1>/Pulse Generator1' */
      stageControl_B.PulseGenerator1 = (stageControl_DW.clockTickCounter_n <
        stageControl_P.PulseGenerator1_Duty) &&
        (stageControl_DW.clockTickCounter_n >= 0L) ?
        stageControl_P.PulseGenerator1_Amp : 0.0;

      /* DiscretePulseGenerator: '<S1>/Pulse Generator1' */
      if (stageControl_DW.clockTickCounter_n >=
          stageControl_P.PulseGenerator1_Period - 1.0) {
        stageControl_DW.clockTickCounter_n = 0L;
      } else {
        stageControl_DW.clockTickCounter_n++;
      }

      /* S-Function (c280xgpio_do): '<S1>/RED LED' */
      {
        if (stageControl_B.PulseGenerator1)
          GpioDataRegs.GPBSET.bit.GPIO34 = 1;
        else
          GpioDataRegs.GPBCLEAR.bit.GPIO34 = 1;
      }
    }

    /* Update for DiscreteStateSpace: '<S1>/Discrete State-Space' */
    {
      real_T xnew[2];
      xnew[0] = (stageControl_P.DiscreteStateSpace_A[0])*
        stageControl_DW.DiscreteStateSpace_DSTATE[0]
        + (stageControl_P.DiscreteStateSpace_A[1])*
        stageControl_DW.DiscreteStateSpace_DSTATE[1];
      xnew[0] += stageControl_P.DiscreteStateSpace_B*rtb_input1;
      xnew[1] = (stageControl_P.DiscreteStateSpace_A[2])*
        stageControl_DW.DiscreteStateSpace_DSTATE[0];
      (void) memcpy(&stageControl_DW.DiscreteStateSpace_DSTATE[0], xnew,
                    sizeof(real_T)*2);
    }

    srUpdateBC(stageControl_DW.PIDcontrol_SubsysRanBC);
  }

  /* End of Outputs for SubSystem: '<Root>/PID control' */

  /* DataTypeConversion: '<Root>/Cast To Single1' */
  stageControl_B.CastToSingle1 = (real32_T)stageControl_B.ManualSwitch;

  /* S-Function (c28xipc_tx): '<Root>/IPC Transmit2' */
  MW_IPC_Transmit(CHANNEL1, (uint32_t *)&stageControl_B.CastToSingle1, 1, 8, 0,
                  0);

  /* ManualSwitch: '<Root>/Manual Switch' incorporates:
   *  Constant: '<Root>/Constant1'
   */
  if (stageControl_P.ManualSwitch_CurrentSetting_b == 1U) {
    rtb_inputMIN = stageControl_B.ManualSwitch;
  } else {
    rtb_inputMIN = stageControl_P.simulinkSettings.inputForce *
      stageControl_P.simulinkSettings.matchingGain;
  }

  /* End of ManualSwitch: '<Root>/Manual Switch' */

  /* Gain: '<S2>/Gain' */
  rtb_inputMIN *= stageControl_P.Gain_Gain;

  /* Saturate: '<S2>/Saturation' */
  if (rtb_inputMIN > stageControl_P.Saturation_UpperSat_d) {
    rtb_PulseGenerator2 = stageControl_P.Saturation_UpperSat_d;
  } else if (rtb_inputMIN < stageControl_P.Saturation_LowerSat_p) {
    rtb_PulseGenerator2 = stageControl_P.Saturation_LowerSat_p;
  } else {
    rtb_PulseGenerator2 = rtb_inputMIN;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* MATLABSystem: '<S2>/DAC' */
  MW_C2000DACSat(0U, rtb_PulseGenerator2);

  /* Saturate: '<S2>/Saturation1' */
  if (rtb_inputMIN > stageControl_P.Saturation1_UpperSat) {
    rtb_inputMIN = stageControl_P.Saturation1_UpperSat;
  } else if (rtb_inputMIN < stageControl_P.Saturation1_LowerSat) {
    rtb_inputMIN = stageControl_P.Saturation1_LowerSat;
  }

  /* End of Saturate: '<S2>/Saturation1' */

  /* MATLABSystem: '<S2>/DAC1' incorporates:
   *  Abs: '<S2>/Abs'
   */
  MW_C2000DACSat(1U, fabs(rtb_inputMIN));

  /* Update for Sin: '<Root>/Sine Wave1' */
  stageControl_DW.counter++;
  if (1.0 / stageControl_P.simulinkSettings.referenceSampleTime /
      stageControl_P.simulinkSettings.referenceFrequency ==
      stageControl_DW.counter) {
    stageControl_DW.counter = 0L;
  }

  /* End of Update for Sin: '<Root>/Sine Wave1' */
  {                                    /* Sample time: [0.0001s, 0.0s] */
    extmodeErrorCode_T errorCode = EXTMODE_SUCCESS;
    extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T)
      ((stageControl_M->Timing.clockTick0 * 1) + 0)
      ;

    /* Trigger External Mode event */
    errorCode = extmodeEvent(0,currentTime);
    if (errorCode != EXTMODE_SUCCESS) {
      /* Code to handle External Mode event errors
         may be added here */
    }
  }

  if (stageControl_M->Timing.TaskCounters.TID[1] == 0) {/* Sample time: [0.01s, 0.0s] */
    extmodeErrorCode_T errorCode = EXTMODE_SUCCESS;
    extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T)
      ((stageControl_M->Timing.clockTick1 * 100) + 0)
      ;

    /* Trigger External Mode event */
    errorCode = extmodeEvent(1,currentTime);
    if (errorCode != EXTMODE_SUCCESS) {
      /* Code to handle External Mode event errors
         may be added here */
    }
  }

  if (stageControl_M->Timing.TaskCounters.TID[2] == 0) {/* Sample time: [0.5s, 0.0s] */
    extmodeErrorCode_T errorCode = EXTMODE_SUCCESS;
    extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T)
      ((stageControl_M->Timing.clockTick2 * 5000) + 0)
      ;

    /* Trigger External Mode event */
    errorCode = extmodeEvent(2,currentTime);
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
  stageControl_M->Timing.taskTime0 =
    ((time_T)(++stageControl_M->Timing.clockTick0)) *
    stageControl_M->Timing.stepSize0;
  if (stageControl_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update absolute timer for sample time: [0.01s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.01, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     */
    stageControl_M->Timing.clockTick1++;
  }

  if (stageControl_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update absolute timer for sample time: [0.5s, 0.0s] */
    /* The "clockTick2" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.5, which is the step size
     * of the task. Size of "clockTick2" ensures timer will not overflow during the
     * application lifespan selected.
     */
    stageControl_M->Timing.clockTick2++;
  }

  rate_scheduler();
}

/* Model initialize function */
void stageControl_initialize(void)
{
  /* Registration code */
  rtmSetTFinal(stageControl_M, -1);
  stageControl_M->Timing.stepSize0 = 0.0001;

  /* External mode info */
  stageControl_M->Sizes.checksums[0] = (793540281U);
  stageControl_M->Sizes.checksums[1] = (925772896U);
  stageControl_M->Sizes.checksums[2] = (2286743409U);
  stageControl_M->Sizes.checksums[3] = (3201996211U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[8];
    stageControl_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)&stageControl_DW.PIDcontrol_SubsysRanBC;
    systemRan[2] = (sysRanDType *)&stageControl_DW.PIDcontrol_SubsysRanBC;
    systemRan[3] = (sysRanDType *)&stageControl_DW.PIDcontrol_SubsysRanBC;
    systemRan[4] = &rtAlwaysEnabled;
    systemRan[5] = &rtAlwaysEnabled;
    systemRan[6] = &rtAlwaysEnabled;
    systemRan[7] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(stageControl_M->extModeInfo,
      &stageControl_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(stageControl_M->extModeInfo,
                        stageControl_M->Sizes.checksums);
    rteiSetTPtr(stageControl_M->extModeInfo, rtmGetTPtr(stageControl_M));
  }

  /* Start for S-Function (c2802xadc): '<Root>/ADC1' */
  if (MW_adcAInitFlag == 0) {
    InitAdcA();
    MW_adcAInitFlag = 1;
  }

  config_ADCA_SOC0 ();

  /* Start for S-Function (c28xipc_tx): '<Root>/IPC Transmit1' */
  IPCInit(CHANNEL2, 1, 0);

  /* Start for S-Function (c2802xadc): '<Root>/ADC2' */
  if (MW_adcDInitFlag == 0) {
    InitAdcD();
    MW_adcDInitFlag = 1;
  }

  config_ADCD_SOC0 ();

  /* Start for S-Function (c280xqep): '<S3>/eQEP' */
  config_QEP_eQEP1(30000U, 15000U, 0, 0, 0, 0, 136, 32768, 32887, 0);

  /* Start for S-Function (c28xipc_tx): '<Root>/IPC Transmit' */
  IPCInit(CHANNEL0, 1, 0);

  /* Start for S-Function (c280xgpio_di): '<Root>/GPIO97 (Button 1)1' */
  EALLOW;
  GpioCtrlRegs.GPBMUX2.all &= 0xFFFCFFFF;
  GpioCtrlRegs.GPBDIR.all &= 0xFEFFFFFF;
  EDIS;

  /* Start for Enabled SubSystem: '<Root>/PID control' */

  /* Start for S-Function (c280xgpio_do): '<S1>/RED LED' */
  EALLOW;
  GpioCtrlRegs.GPBMUX1.all &= 0xFFFFFFCF;
  GpioCtrlRegs.GPBDIR.all |= 0x4;
  EDIS;

  /* End of Start for SubSystem: '<Root>/PID control' */

  /* Start for S-Function (c28xipc_tx): '<Root>/IPC Transmit2' */
  IPCInit(CHANNEL1, 1, 0);

  /* Start for MATLABSystem: '<S2>/DAC' */
  MW_ConfigureDACA();

  /* Start for MATLABSystem: '<S2>/DAC1' */
  MW_ConfigureDACB();

  /* InitializeConditions for S-Function (sdspcount2): '<Root>/Counter1' */
  stageControl_DW.Counter1_ClkEphState = 5UL;
  stageControl_DW.Counter1_Count = stageControl_P.Counter1_InitialCount;

  /* SystemInitialize for Enabled SubSystem: '<Root>/PID control' */
  /* InitializeConditions for DiscreteStateSpace: '<S1>/Discrete State-Space' */
  stageControl_DW.DiscreteStateSpace_DSTATE[0] =
    stageControl_P.DiscreteStateSpace_InitialCondi;
  stageControl_DW.DiscreteStateSpace_DSTATE[1] =
    stageControl_P.DiscreteStateSpace_InitialCondi;

  /* SystemInitialize for ManualSwitch: '<S1>/Manual Switch' incorporates:
   *  Outport: '<S1>/Output'
   */
  stageControl_B.ManualSwitch = stageControl_P.Output_Y0;

  /* End of SystemInitialize for SubSystem: '<Root>/PID control' */
}

/* Model terminate function */
void stageControl_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
