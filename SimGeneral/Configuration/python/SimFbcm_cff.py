import FWCore.ParameterSet.Config as cms

#Full Event content
SimFbcmDigiFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_simFbcmDigis_*_*')
)
