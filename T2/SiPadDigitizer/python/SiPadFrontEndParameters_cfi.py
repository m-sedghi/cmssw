import FWCore.ParameterSet.Config as cms


fftSimParam = cms.PSet(
	NumOfFFT_Points = cms.int32(2048), # Length of signal, This should be an integer number with power of 2
	SamplingRepetition = cms.int32(10) # FS: Sampling repetition per ns [1/ns]
)
TofCharge_Test = cms.PSet(
	TofVector =  cms.vdouble(0.0, 25.2), 
	ChargeVect = cms.vdouble(3*6242, 8000),
	TestSensorSize= cms.double(0.04) # cm2
)

SiHitPulseShapeParam =cms.PSet(
	HitPulseParam =  cms.vdouble(0.6294422, 99.999855, 40.371655, 1.0, 3.5/2.2) # 0.6294422, 99.999855, 40.371655, 1.0, 3.5/2.2
) 

SiPadFrontEndBlock0 = cms.PSet(
	GoodForSizeRange = cms.vdouble(0.0,0.0255), # cm2, range from minimum size through maximum size
	
	MaxFEOutputVoltage = cms.double(700.0), # mV
	#LimmiterEdgeCorrFactor = cms.double(1.5), # unitless, by default should be 1.5 
	
	ZCComp_LowerTsh = cms.double(-5.0), # mV
	ZCComp_UpperTsh = cms.double(0.0), # mV
	ArmingComp_LowerTsh = cms.double(5.0), # mV
	ArmingComp_UpperTsh = cms.double(20.0), # mV
	
	TIA_Shaper_Gain = cms.double(28.0), # V/V (the amplifier gain after TIA and Shaper1)
	Tia_Rf = cms.double(5.0), # kOhm
	Tia_Cf = cms.double(0.25), # pF
	Tia_Cin_gs = cms.double(0.4), # pf (just the TIA input Capacitance), the SiPad Capacitance will be added to this
	Tia_Co = cms.double(0.4), # pf
	Tia_gin = cms.double(3.0), # mS
	SensorCouplingCapacitance = cms.double(315.0), # pF
	SensorCapPerCm2 = cms.double(86.207), # pF/cm2
	Shaper1_Tau = cms.double(0.9), # ns
	CFD_Delay  = cms.double(2.0), # ns
	CfdShaper_Gain = cms.double(1.5), # V/V
	CfdShaper_Tau = cms.double(0.25), # ns
	DelayModel = cms.string('FirstOrderPadeDelay'), # 'IdealDelay' or 'FirstOrderPadeDelay'
	CFD_Fraction = cms.double(0.5), # between 0-1, typically around 0.5
	
	lpGBT_AlignerDelay = cms.double(5.2), # ns
	Bx_Duration = cms.double(25.0), # ns
	
    ToAUpperCut = cms.double(30.0), # ns // for BIB study, more than one BX should be investigated
	ToALowerCut = cms.double(-30.0), # ns // for BIB study, more than one BX should be investigated
	BinLength = cms.double(6.26), # ns
	BinOffset = cms.double(0.0), # ns
	
)

SiPadFrontEndBlock1 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock2 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock3 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock4 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock5 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock6 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock7 = SiPadFrontEndBlock0.clone();
SiPadFrontEndBlock8 = SiPadFrontEndBlock0.clone();

#----------------------------

SiPadFrontEndBlock0.GoodForSizeRange = cms.vdouble(0.0   , 0.0255) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock1.GoodForSizeRange = cms.vdouble(0.0255, 0.0335) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock2.GoodForSizeRange = cms.vdouble(0.0335, 0.046) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock3.GoodForSizeRange = cms.vdouble(0.046 , 0.067) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock4.GoodForSizeRange = cms.vdouble(0.067 , 0.1065) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock5.GoodForSizeRange = cms.vdouble(0.1065, 0.1965) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock6.GoodForSizeRange = cms.vdouble(0.1965, 0.491) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock7.GoodForSizeRange = cms.vdouble(0.491 , 0.866) # cm2, range from minimum size through maximum size
SiPadFrontEndBlock8.GoodForSizeRange = cms.vdouble(0.866 , 100.0) # cm2, range from minimum size through maximum size
#-------------------------


SiPadFrontEndBlock1.CFD_Fraction = cms.double(0.55)  
SiPadFrontEndBlock1.lpGBT_AlignerDelay = cms.double(5.5) # ns
SiPadFrontEndBlock1.TIA_Shaper_Gain = cms.double(28.0)

SiPadFrontEndBlock2.CFD_Fraction = cms.double(0.6)  
SiPadFrontEndBlock2.lpGBT_AlignerDelay = cms.double(5.8) # ns
SiPadFrontEndBlock2.TIA_Shaper_Gain = cms.double(28.0)

SiPadFrontEndBlock3.CFD_Fraction = cms.double(0.65)  
SiPadFrontEndBlock3.lpGBT_AlignerDelay = cms.double(6.2) # ns
SiPadFrontEndBlock3.TIA_Shaper_Gain = cms.double(31.0)

SiPadFrontEndBlock4.CFD_Fraction = cms.double(0.7)  
SiPadFrontEndBlock4.lpGBT_AlignerDelay = cms.double(6.6) # ns
SiPadFrontEndBlock4.TIA_Shaper_Gain = cms.double(36.0)

SiPadFrontEndBlock5.CFD_Fraction = cms.double(0.75)  
SiPadFrontEndBlock5.lpGBT_AlignerDelay = cms.double(7.1) # ns
SiPadFrontEndBlock5.TIA_Shaper_Gain = cms.double(46.0)

SiPadFrontEndBlock6.CFD_Fraction = cms.double(0.8)  
SiPadFrontEndBlock6.lpGBT_AlignerDelay = cms.double(7.8) # ns
SiPadFrontEndBlock6.TIA_Shaper_Gain = cms.double(73.0)

SiPadFrontEndBlock7.CFD_Fraction = cms.double(0.85)  
SiPadFrontEndBlock7.lpGBT_AlignerDelay = cms.double(8.6) # ns
SiPadFrontEndBlock7.TIA_Shaper_Gain = cms.double(167.0)

SiPadFrontEndBlock8.CFD_Fraction = cms.double(0.88)  
SiPadFrontEndBlock8.lpGBT_AlignerDelay = cms.double(9.1) # ns
SiPadFrontEndBlock8.TIA_Shaper_Gain = cms.double(234.0) 


