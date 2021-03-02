import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2_cff import Phase2
from Configuration.Eras.Modifier_fbcmDigi_cff import fbcmDigi

process = cms.Process('SIM',Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D80Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D80_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
#process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


# randomeze the seeds every time cmsRun is invoked
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    #output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
#process.source = cms.Source("EmptySource")
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(
							'file:/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput/BeamHalo/GenBeamHalo_0.root',
							#'file:/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput/BeamGasOxygen/GenBeamGasOxygen_0.root',
							#'file:/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput/BeamGasHydrogen/GenBeamGasHydrogen_0.root',
							#'file:/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput/BeamGasCarbon/GenBeamGasCarbon_0.root',
							),
                            skipEvents = cms.untracked.uint32(0),
							#eventsToProcess = cms.untracked.VEventRange("1:{}-1:{}".format(range_min, range_max))
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    #SkipEvent = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('BIB particle G4 simulation similar to MinBias'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    #SelectEvents = cms.untracked.PSet(
    #    SelectEvents = cms.vstring('generation_step')
    #),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:BIB_SIM.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

import SimG4Core.Application.g4SimHits_cfi
process.g4SimHits.Generator.ApplyEtaCuts  = cms.bool(False)
process.g4SimHits.Generator.MinPCut  = cms.double(0.0001) #100keV
process.g4SimHits.Generator.BeamBkgdEvent = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SaveFirstLevelSecondary = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInCalo = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInMuon = cms.untracked.bool(True)



# Other statements
#process.genstepfilter.triggerConditions=cms.vstring("generation_step")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# tuning SensitiveDetector module of FBCM
#process.g4SimHits.FbcmSD.ZeroEnergyLoss = cms.bool(True) # by default it shoule be False
process.g4SimHits.FbcmSD.EnergyThresholdForPersistencyInGeV =  cms.double(0.0001)

# Path and EndPath definitions

#process.simulation_step = cms.Path(process.psim+process.mix)
process.simulation_step = cms.Path(process.psim)
#process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(#process.generation_step,
								#process.genfiltersummary_step,
								process.simulation_step,
								process.endjob_step,
								process.FEVTDEBUGoutput_step)

#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# filter all path with the production filter sequence
#for path in process.paths:
#	getattr(process,path).insert(0, process.ProductionFilterSequence)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 
#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
