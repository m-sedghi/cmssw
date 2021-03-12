import FWCore.ParameterSet.Config as cms

from SimFbcm.SiPadDigitizer.SiPadSimParameters_cfi import SiPadSimBlock
from SimFbcm.SiPadDigitizer.SiPadFrontEndParameters_cfi import *


SiPadDigitizer = cms.PSet(
    accumulatorType = cms.string("SiPadDigitizer"), 
    hitsProducer = cms.string('g4SimHits'), 
	SubdetName=cms.string('FBCMHits'), 
	GeometryType = cms.string('idealForDigi'), 
	InstanceName = cms.string('SiPad'), 
	SiPadSimParam = cms.PSet(SiPadSimBlock),
	FFT_SimParam = cms.PSet(fftSimParam),
	SiPadFrontEndParam = cms.VPSet( SiPadFrontEndBlock0,
									SiPadFrontEndBlock1,
									SiPadFrontEndBlock2,
									SiPadFrontEndBlock3,
									SiPadFrontEndBlock4,
									SiPadFrontEndBlock5,
									SiPadFrontEndBlock6,
									SiPadFrontEndBlock7,
									SiPadFrontEndBlock8 ),
	TofCharge_Test = cms.PSet(TofCharge_Test),
	SiHitPulseShapeParam = cms.PSet(SiHitPulseShapeParam)
	
)

