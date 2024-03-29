#include "c2000BoardSupport.h"
#include "MW_f2837xD_includes.h"
#include "rtwtypes.h"
#include "stageControl.h"
#include "stageControl_private.h"

void config_ADCA_SOC0()
{
  EALLOW;
  AdcaRegs.ADCSOC0CTL.bit.CHSEL = 15;  /* Set SOC0 channel select to ADCIN15*/
  AdcaRegs.ADCSOC0CTL.bit.TRIGSEL = 0;
  AdcaRegs.ADCSOC0CTL.bit.ACQPS = 14.0;
                               /* Set SOC0 S/H Window to 15.0 ADC Clock Cycles*/
  AdcaRegs.ADCINTSOCSEL1.bit.SOC0 = 0;
                                   /* SOCx No ADCINT Interrupt Trigger Select.*/
  AdcaRegs.ADCOFFTRIM.bit.OFFTRIM = AdcaRegs.ADCOFFTRIM.bit.OFFTRIM;/* Set Offset Error Correctino Value*/
  AdcaRegs.ADCCTL1.bit.INTPULSEPOS = 1;
                                /* Late interrupt pulse trips AdcResults latch*/
  AdcaRegs.ADCSOCPRICTL.bit.SOCPRIORITY = 0;/* All in round robin mode SOC Priority*/
  EDIS;
}

void config_ADCD_SOC0()
{
  EALLOW;
  AdcdRegs.ADCSOC0CTL.bit.CHSEL = 0;
                               /* Set SOC0 channel select to ADCIN0 and ADCIN1*/
  AdcdRegs.ADCSOC0CTL.bit.TRIGSEL = 0;
  AdcdRegs.ADCSOC0CTL.bit.ACQPS = 14.0;
                               /* Set SOC0 S/H Window to 15.0 ADC Clock Cycles*/
  AdcdRegs.ADCINTSOCSEL1.bit.SOC0 = 0;
                                   /* SOCx No ADCINT Interrupt Trigger Select.*/
  AdcdRegs.ADCOFFTRIM.bit.OFFTRIM = AdcdRegs.ADCOFFTRIM.bit.OFFTRIM;/* Set Offset Error Correctino Value*/
  AdcdRegs.ADCCTL1.bit.INTPULSEPOS = 1;
                                /* Late interrupt pulse trips AdcResults latch*/
  AdcdRegs.ADCSOCPRICTL.bit.SOCPRIORITY = 0;/* All in round robin mode SOC Priority*/
  EDIS;
}

void InitAdcA()
{
  EALLOW;
  CpuSysRegs.PCLKCR13.bit.ADC_A = 1;
  AdcaRegs.ADCCTL2.bit.PRESCALE = 8;
  AdcSetMode(ADC_ADCA, ADC_RESOLUTION_12BIT, ADC_SIGNALMODE_SINGLE);

  //power up the ADC
  AdcaRegs.ADCCTL1.bit.ADCPWDNZ = 1;

  //delay for 1ms to allow ADC time to power up
  DELAY_US(1000);
  EDIS;
}

void InitAdcD()
{
  EALLOW;
  CpuSysRegs.PCLKCR13.bit.ADC_D = 1;
  AdcdRegs.ADCCTL2.bit.PRESCALE = 8;
  AdcSetMode(ADC_ADCD, ADC_RESOLUTION_16BIT, ADC_SIGNALMODE_DIFFERENTIAL);

  //power up the ADC
  AdcdRegs.ADCCTL1.bit.ADCPWDNZ = 1;

  //delay for 1ms to allow ADC time to power up
  DELAY_US(1000);
  EDIS;
}
