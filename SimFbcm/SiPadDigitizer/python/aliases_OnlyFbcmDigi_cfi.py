import FWCore.ParameterSet.Config as cms

simFbcmDigis = cms.EDAlias(
    mix = cms.VPSet(
      cms.PSet(type = cms.string('SiPadDigiDataedmDetSetVector'))
    )
)

simCastorDigis = cms.EDAlias()

simEcalUnsuppressedDigis = cms.EDAlias()

simHcalUnsuppressedDigis = cms.EDAlias()

simSiPixelDigis = cms.EDAlias()

simSiStripDigis = cms.EDAlias()

simHGCalUnsuppressedDigis = cms.EDAlias()

simHFNoseUnsuppressedDigis = cms.EDAlias()

simAPVsaturation = cms.EDAlias()
