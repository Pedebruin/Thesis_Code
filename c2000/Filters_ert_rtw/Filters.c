/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Filters.c
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

#include "Filters.h"
#include "Filters_private.h"
#include "rtwtypes.h"
#include "xcp.h"
#include "ext_mode.h"

extmodeSimulationTime_T currentTime = (extmodeSimulationTime_T) 0;

/* Block signals (default storage) */
B_Filters_T Filters_B;

/* Block states (default storage) */
DW_Filters_T Filters_DW;

/* Real-time model */
static RT_MODEL_Filters_T Filters_M_;
RT_MODEL_Filters_T *const Filters_M = &Filters_M_;

/* Hardware Interrupt Block: '<Root>/C28x Hardware Interrupt' */
void isr_int1pie13_task_fcn(void)
{
  if (1 == runModel) {
    /* Call the system: <Root>/Function-Call Subsystem */
    {
      /* Reset subsysRan breadcrumbs */
      srClearBC(Filters_DW.FunctionCallSubsystem_SubsysRan);

      /* S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */

      /* Output and update for function-call system: '<Root>/Function-Call Subsystem' */

      /* Asynchronous task reads absolute time. Data (absolute time)
         transfers from low priority task (base rate) to high priority
         task (asynchronous rate). Double buffers are used to ensure
         data integrity.
         -- rtmL2HLastBufWr is the buffer index that is written last.
       */
      Filters_M->Timing.clockTick1 = Filters_M->
        Timing.rtmL2HDbBufClockTick[Filters_M->Timing.rtmL2HLastBufWr];

      /* S-Function (c28xipc_rx): '<S1>/IPC Receive' */
      MW_IPC_Receive(CHANNEL0, (uint32_t *)&Filters_B.IPCReceive_o1_l,
                     &Filters_B.IPCReceive_o2_d, 1, 8, 0, 0);
      Filters_DW.FunctionCallSubsystem_SubsysRan = 4;

      /* End of Outputs for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */
    }

    currentTime = (extmodeSimulationTime_T) Filters_M->Timing.clockTick0;
    extmodeEvent(1,currentTime);
  }
}

/* Hardware Interrupt Block: '<Root>/C28x Hardware Interrupt' */
void isr_int1pie14_task_fcn(void)
{
  if (1 == runModel) {
    /* Call the system: <Root>/Function-Call Subsystem1 */
    {
      /* Reset subsysRan breadcrumbs */
      srClearBC(Filters_DW.FunctionCallSubsystem1_SubsysRa);

      /* S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */

      /* Output and update for function-call system: '<Root>/Function-Call Subsystem1' */

      /* Asynchronous task reads absolute time. Data (absolute time)
         transfers from low priority task (base rate) to high priority
         task (asynchronous rate). Double buffers are used to ensure
         data integrity.
         -- rtmL2HLastBufWr is the buffer index that is written last.
       */
      Filters_M->Timing.clockTick2 = Filters_M->
        Timing.rtmL2HDbBufClockTick[Filters_M->Timing.rtmL2HLastBufWr];

      /* S-Function (c28xipc_rx): '<S2>/IPC Receive' */
      MW_IPC_Receive(CHANNEL1, (uint32_t *)&Filters_B.IPCReceive_o1_j,
                     &Filters_B.IPCReceive_o2_l, 1, 4, 0, 0);
      Filters_DW.FunctionCallSubsystem1_SubsysRa = 4;

      /* End of Outputs for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */
    }

    currentTime = (extmodeSimulationTime_T) Filters_M->Timing.clockTick0;
    extmodeEvent(2,currentTime);
  }
}

/* Hardware Interrupt Block: '<Root>/C28x Hardware Interrupt' */
void isr_int1pie15_task_fcn(void)
{
  if (1 == runModel) {
    /* Call the system: <Root>/Function-Call Subsystem2 */
    {
      /* Reset subsysRan breadcrumbs */
      srClearBC(Filters_DW.FunctionCallSubsystem2_SubsysRa);

      /* S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */

      /* Output and update for function-call system: '<Root>/Function-Call Subsystem2' */

      /* Asynchronous task reads absolute time. Data (absolute time)
         transfers from low priority task (base rate) to high priority
         task (asynchronous rate). Double buffers are used to ensure
         data integrity.
         -- rtmL2HLastBufWr is the buffer index that is written last.
       */
      Filters_M->Timing.clockTick3 = Filters_M->
        Timing.rtmL2HDbBufClockTick[Filters_M->Timing.rtmL2HLastBufWr];

      /* S-Function (c28xipc_rx): '<S3>/IPC Receive' */
      MW_IPC_Receive(CHANNEL2, (uint32_t *)&Filters_B.IPCReceive_o1,
                     &Filters_B.IPCReceive_o2, 1, 8, 0, 0);
      Filters_DW.FunctionCallSubsystem2_SubsysRa = 4;

      /* End of Outputs for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */
    }

    currentTime = (extmodeSimulationTime_T) Filters_M->Timing.clockTick0;
    extmodeEvent(3,currentTime);
  }
}

/* Model step function */
void Filters_step(void)
{
  /* DiscretePulseGenerator: '<Root>/Pulse Generator1' */
  Filters_B.PulseGenerator1 = (Filters_DW.clockTickCounter <
    Filters_P.PulseGenerator1_Duty) && (Filters_DW.clockTickCounter >= 0L) ?
    Filters_P.PulseGenerator1_Amp : 0.0;

  /* DiscretePulseGenerator: '<Root>/Pulse Generator1' */
  if (Filters_DW.clockTickCounter >= Filters_P.PulseGenerator1_Period - 1.0) {
    Filters_DW.clockTickCounter = 0L;
  } else {
    Filters_DW.clockTickCounter++;
  }

  /* S-Function (c280xgpio_do): '<Root>/BLUE LED' */
  {
    if (Filters_B.PulseGenerator1)
      GpioDataRegs.GPASET.bit.GPIO31 = 1;
    else
      GpioDataRegs.GPACLEAR.bit.GPIO31 = 1;
  }

  {                                    /* Sample time: [0.5s, 0.0s] */
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   */
  Filters_M->Timing.taskTime0 =
    ((time_T)(++Filters_M->Timing.clockTick0)) * Filters_M->Timing.stepSize0;

  {
    /* Base rate updates double buffers of absolute time for
       asynchronous task. Double buffers are used to ensure
       data integrity when asynchronous task reads absolute
       time.
       -- rtmL2HLastBufWr is the buffer index that is written last.
     */
    boolean_T bufIdx = !(Filters_M->Timing.rtmL2HLastBufWr != 0U);
    Filters_M->Timing.rtmL2HDbBufClockTick[bufIdx] =
      Filters_M->Timing.clockTick0;
    Filters_M->Timing.rtmL2HLastBufWr = bufIdx ? 1U : 0U;
  }
}

/* Model initialize function */
void Filters_initialize(void)
{
  /* Registration code */
  rtmSetTFinal(Filters_M, -1);
  Filters_M->Timing.stepSize0 = 0.5;

  /* External mode info */
  Filters_M->Sizes.checksums[0] = (2686384147U);
  Filters_M->Sizes.checksums[1] = (382012995U);
  Filters_M->Sizes.checksums[2] = (3562166748U);
  Filters_M->Sizes.checksums[3] = (1625091289U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[4];
    Filters_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)&Filters_DW.FunctionCallSubsystem_SubsysRan;
    systemRan[2] = (sysRanDType *)&Filters_DW.FunctionCallSubsystem1_SubsysRa;
    systemRan[3] = (sysRanDType *)&Filters_DW.FunctionCallSubsystem2_SubsysRa;
    rteiSetModelMappingInfoPtr(Filters_M->extModeInfo,
      &Filters_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(Filters_M->extModeInfo, Filters_M->Sizes.checksums);
    rteiSetTPtr(Filters_M->extModeInfo, rtmGetTPtr(Filters_M));
  }

  /* Start for S-Function (c280xgpio_do): '<Root>/BLUE LED' */
  {
    Uint32 *pulMsgRam = (void *)CPU01_TO_CPU02_PASSMSG;
    Uint32 ulRWord32 = 0;
    Uint32 gpioData = 0x00010000;
    gpioData = gpioData|31;

#ifndef GPIO31UsedByCLA

    {
      IPCLiteLtoRFunctionCall(IPC_FLAG0, pulMsgRam[0], gpioData, IPC_FLAG31);
      while (IPCLiteLtoRGetResult(&ulRWord32,IPC_LENGTH_32_BITS,
              IPC_FLAG31) != STATUS_PASS) {
      }
    }

#endif

  }

  /* Start for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' incorporates:
   *  SubSystem: '<Root>/Function-Call Subsystem'
   */
  /* Start for function-call system: '<Root>/Function-Call Subsystem' */

  /* Start for S-Function (c28xipc_rx): '<S1>/IPC Receive' */
  IPCInit(CHANNEL0, 1, 0);
  ;

  /* Start for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' incorporates:
   *  SubSystem: '<Root>/Function-Call Subsystem1'
   */
  /* Start for function-call system: '<Root>/Function-Call Subsystem1' */

  /* Start for S-Function (c28xipc_rx): '<S2>/IPC Receive' */
  IPCInit(CHANNEL1, 1, 0);
  ;

  /* Start for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' incorporates:
   *  SubSystem: '<Root>/Function-Call Subsystem2'
   */

  /* Start for function-call system: '<Root>/Function-Call Subsystem2' */

  /* Start for S-Function (c28xipc_rx): '<S3>/IPC Receive' */
  IPCInit(CHANNEL2, 1, 0);
  ;

  /* End of Start for S-Function (c28xisr_c2000): '<Root>/C28x Hardware Interrupt' */
}

/* Model terminate function */
void Filters_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
