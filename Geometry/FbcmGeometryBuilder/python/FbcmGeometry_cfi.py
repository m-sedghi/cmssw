import FWCore.ParameterSet.Config as cms

#
# This cfi should be included to build the FBCM geometry model.
#
FbcmGeometryESModule = cms.ESProducer("FbcmGeometryESModule",
    useDDD = cms.bool(True),
    alignmentsLabel = cms.string(''),
    applyAlignment = cms.bool(False)
)
