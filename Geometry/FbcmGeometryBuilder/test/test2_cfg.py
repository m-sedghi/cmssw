import FWCore.ParameterSet.Config as cms

process = cms.Process('Test2')

process.source = cms.Source('EmptySource')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )


process.load('Configuration.Geometry.GeometryExtended2026D60_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#Note: only use one of the following:
#process.load('Geometry.FbcmGeometryBuilder.FbcmGeometry_cfi')
process.FbcmGeometryESProducer = cms.ESProducer("FbcmGeometryESModule",
										useDDD = cms.bool(True),
										alignmentsLabel = cms.string('The1034!'),
										applyAlignment = cms.bool(False)
									)

process.TestAnalyzer = cms.EDAnalyzer("FbcmGeometryTest")

process.p = cms.Path(process.TestAnalyzer)
