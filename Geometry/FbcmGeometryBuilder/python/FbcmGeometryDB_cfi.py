import FWCore.ParameterSet.Config as cms

#
# This cfi should be included to build the FBCM geometry model.
# reading From DB has not implemented yet
FbcmGeometryESModule = cms.ESProducer("FbcmGeometryESModule",
    useDDD = cms.bool(False),
    alignmentsLabel = cms.string(''),
    applyAlignment = cms.bool(False)
)
