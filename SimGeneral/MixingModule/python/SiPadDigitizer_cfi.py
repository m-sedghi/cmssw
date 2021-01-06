import FWCore.ParameterSet.Config as cms

from SimGeneral.MixingModule.SiPadSimParameters_cfi import SiPadSimBlock
from SimGeneral.MixingModule.SiPadFrontEndParameters_cfi import *

SiPadDigitizer = cms.PSet(
    accumulatorType = cms.string("SiPadDigitizer"), 
    hitsProducer = cms.string('g4SimHits'), 
	SubdetName=cms.string('FBCMHits'), 
	GeometryType = cms.string('idealForDigi'), 
	InstanceName = cms.string('SiPad'), 
	#SimFBCM_Tag = cms.InputTag("SimFBCM", "SiPad"),
	#isReadoutAnalog = cms.bool(True),#set this to true if you want analog readout for OT
	#makeDigiSimLinks = cms.untracked.bool(True),
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

# from Configuration.ProcessModifiers.premix_stage1_cff import premix_stage1
# #premix_stage1.toModify(SiPadDigitizer, makeDigiSimLinks = False)

# # Customize here instead of SiPadSimBlock as the latter is imported
# # also to DataMixer configuration, and the original version is needed
# # there in stage2. Customize before phase2_tracker because this
# # customization applies only to phase0/1 pixel.
# from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
# premix_stage2.toModify(SiPadDigitizer,
    # AddPixelInefficiency = False # will be added in DataMixer
# )

# from SimTracker.SiPhase2Digitizer.phase2TrackerDigitizer_cfi import phase2TrackerDigitizer as _phase2TrackerDigitizer, _premixStage1ModifyDict
# from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
# #phase2_tracker.toReplaceWith(SiPadDigitizer, _phase2TrackerDigitizer.clone()) # have to clone here in order to not change the original with further customizations

# # Customize here instead of phase2TrackerDigitizer as the latter is
# # imported also to DataMixer configuration, and the original version
# # is needed there in stage2.
# #(premix_stage2 & phase2_tracker).toModify(SiPadDigitizer, **_premixStage1ModifyDict)
# from CalibTracker.SiPixelESProducers.PixelFEDChannelCollectionProducer_cfi import *
