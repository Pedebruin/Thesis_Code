#include "c2000BoardSupport.h"
#include "MW_f2837xD_includes.h"
#include "rtwtypes.h"
#include "Identification.h"
#include "Identification_private.h"

void init_board ()
{
  DisableDog();
  EALLOW;
  CpuSysRegs.PCLKCR0.bit.DMA = 1;
  CpuSysRegs.PCLKCR6.bit.SD1 = 1;
  CpuSysRegs.PCLKCR6.bit.SD2 = 1;
  EDIS;

#ifdef CPU1

  EALLOW;

  //enable pull-ups on unbonded IOs as soon as possible to reduce power consumption.
  GPIO_EnableUnbondedIOPullups();
  CpuSysRegs.PCLKCR13.bit.ADC_A = 1;
  CpuSysRegs.PCLKCR13.bit.ADC_B = 1;
  CpuSysRegs.PCLKCR13.bit.ADC_C = 1;
  CpuSysRegs.PCLKCR13.bit.ADC_D = 1;

  //check if device is trimmed
  if (*((Uint16 *)0x5D1B6) == 0x0000) {
    //device is not trimmed, apply static calibration values
    AnalogSubsysRegs.ANAREFTRIMA.all = 31709;
    AnalogSubsysRegs.ANAREFTRIMB.all = 31709;
    AnalogSubsysRegs.ANAREFTRIMC.all = 31709;
    AnalogSubsysRegs.ANAREFTRIMD.all = 31709;
  }

  CpuSysRegs.PCLKCR13.bit.ADC_A = 0;
  CpuSysRegs.PCLKCR13.bit.ADC_B = 0;
  CpuSysRegs.PCLKCR13.bit.ADC_C = 0;
  CpuSysRegs.PCLKCR13.bit.ADC_D = 0;
  EDIS;
  InitSysPll(XTAL_OSC,40,0,1);

  //Turn on all peripherals
  //InitPeripheralClocks();
  EALLOW;
  CpuSysRegs.PCLKCR0.bit.CPUTIMER0 = 1;
  CpuSysRegs.PCLKCR0.bit.CPUTIMER1 = 1;
  CpuSysRegs.PCLKCR0.bit.CPUTIMER2 = 1;
  CpuSysRegs.PCLKCR0.bit.HRPWM = 1;
  CpuSysRegs.PCLKCR1.bit.EMIF1 = 1;
  CpuSysRegs.PCLKCR1.bit.EMIF2 = 1;

  /* Assign all Peripherals to CPU2 */
  DevCfgRegs.CPUSEL11.all = 0x0000000F;
  DevCfgRegs.CPUSEL0.all = 0x00000FFF;
  DevCfgRegs.CPUSEL1.all = 0x0000003F;
  DevCfgRegs.CPUSEL2.all = 0x00000007;
  DevCfgRegs.CPUSEL5.all = 0x0000000F;
  DevCfgRegs.CPUSEL6.all = 0x0000000F;
  DevCfgRegs.CPUSEL8.all = 0x00000003;
  DevCfgRegs.CPUSEL14.all = 0x00070000;
  DevCfgRegs.CPUSEL7.all = 0x00000003;
  DevCfgRegs.CPUSEL12.all = 0x000000FF;
  DevCfgRegs.CPUSEL4.all = 0x00000003;

#ifdef MW_DAC_CHANNEL_A

  DevCfgRegs.CPUSEL14.bit.DAC_A = 0;

#endif

#ifdef MW_DAC_CHANNEL_B

  DevCfgRegs.CPUSEL14.bit.DAC_B = 0;

#endif

#ifdef MW_DAC_CHANNEL_C

  DevCfgRegs.CPUSEL14.bit.DAC_C = 0;

#endif

  /* Assign SDFM modules to CPU1 */
#ifdef MW_SDFM_1

  DevCfgRegs.CPUSEL4.bit.SD1 = 0;

#endif

#ifdef MW_SDFM_2

  DevCfgRegs.CPUSEL4.bit.SD2 = 0;

#endif

  /* Assign used eQEP modules to CPU1 */
  DevCfgRegs.CPUSEL2.bit.EQEP1 = 0;

#if MW_EXTMODE_SCIMODULE == 0

  DevCfgRegs.CPUSEL5.bit.SCI_A = 0;

#elif MW_EXTMODE_SCIMODULE == 1

  DevCfgRegs.CPUSEL5.bit.SCI_B = 0;

#elif MW_EXTMODE_SCIMODULE == 2

  DevCfgRegs.CPUSEL5.bit.SCI_C = 0;

#elif MW_EXTMODE_SCIMODULE == 3

  DevCfgRegs.CPUSEL5.bit.SCI_D = 0;

#endif

  /* Assign used SPI modules to CPU1 */
#ifdef MW_SPI_A

  DevCfgRegs.CPUSEL6.bit.SPI_A = 0;

#endif

#ifdef MW_SPI_B

  DevCfgRegs.CPUSEL6.bit.SPI_B = 0;

#endif

#ifdef MW_SPI_C

  DevCfgRegs.CPUSEL6.bit.SPI_C = 0;

#endif

#ifdef MW_SPI_D

  DevCfgRegs.CPUSEL6.bit.SPI_D = 0;

#endif

#if defined MW_CMPSS1_COMPH || defined MW_CMPSS1_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS1 = 0;

#endif

#if defined MW_CMPSS2_COMPH || defined MW_CMPSS2_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS2 = 0;

#endif

#if defined MW_CMPSS3_COMPH || defined MW_CMPSS3_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS3 = 0;

#endif

#if defined MW_CMPSS4_COMPH || defined MW_CMPSS4_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS4 = 0;

#endif

#if defined MW_CMPSS5_COMPH || defined MW_CMPSS5_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS5 = 0;

#endif

#if defined MW_CMPSS6_COMPH || defined MW_CMPSS6_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS6 = 0;

#endif

#if defined MW_CMPSS7_COMPH || defined MW_CMPSS7_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS7 = 0;

#endif

#if defined MW_CMPSS8_COMPH || defined MW_CMPSS8_COMPL

  DevCfgRegs.CPUSEL12.bit.CMPSS8 = 0;

#endif

  EDIS;

#endif

  EALLOW;

  /* Configure low speed peripheral clocks */
  ClkCfgRegs.LOSPCP.bit.LSPCLKDIV = 0U;
  EDIS;

  /* Disable and clear all CPU interrupts */
  DINT;
  IER = 0x0000;
  IFR = 0x0000;
  InitPieCtrl();
  InitPieVectTable();
  initSetGPIOIPC();
  InitCpuTimers();

#ifdef CPU1

  /* initial GPIO qualification settings.... */
  EALLOW;
  GpioCtrlRegs.GPAQSEL1.all = 0x0;
  GpioCtrlRegs.GPAQSEL2.all = 0x0;
  GpioCtrlRegs.GPBQSEL1.all = 0x0;
  GpioCtrlRegs.GPBQSEL2.all = 0x0;
  GpioCtrlRegs.GPCQSEL1.all = 0x0;
  GpioCtrlRegs.GPCQSEL2.all = 0x0;
  GpioCtrlRegs.GPDQSEL1.all = 0x0;
  GpioCtrlRegs.GPDQSEL2.all = 0x0;
  GpioCtrlRegs.GPEQSEL1.all = 0x0;
  GpioCtrlRegs.GPEQSEL2.all = 0x0;
  GpioCtrlRegs.GPFQSEL1.all = 0x0;
  EDIS;

#endif

}
