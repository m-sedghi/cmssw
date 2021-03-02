# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --evt_type SingleNuE10_cfi -s GEN,SIM,DIGI --mc --fileout file:BIB_plus_GEN_SIM_DIGI.root --conditions auto:phase2_realistic --pileup_input file:MinBias_14TeV_pythia8_TuneCUETP8M1_GEN_SIM.root --pileup "AVE_200_BX_25ns,{'B':(-3,3),'N':1.5}" --era Phase2,fbcmDigi,OnlyfbcmDigi --datatier GEN-SIM-DIGI-RAW --geometry Extended2026D80 --eventcontent FEVTDEBUG --python_filename BIB_plus_GEN_SIM_DIGI_cfg.py --customise SimFbcm/SiPadDigitizer/aging.no_aging,Configuration/DataProcessing/Utils.addMonitoring --nThreads 2 -n 2 --no_exe
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2_cff import Phase2
from Configuration.Eras.Modifier_fbcmDigi_cff import fbcmDigi
from Configuration.Eras.Modifier_OnlyfbcmDigi_cff import OnlyfbcmDigi

process = cms.Process('DIGI',Phase2,fbcmDigi,OnlyfbcmDigi)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D80Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D80_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
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
    annotation = cms.untracked.string('SingleNuE10_cfi nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:BIB_plus_GEN_SIM_DIGI.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
#process.mix.input.nbPileupEvents.averageNumber = cms.double(1.500000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
#process.mix.input.fileNames = cms.untracked.vstring(['file:MinBias_14TeV_pythia8_TuneCUETP8M1_GEN_SIM.root'])

process.mix.input = cms.SecSource("EmbeddedRootSource",
    fileNames = cms.untracked.vstring('file:MinBias_14TeV_pythia8_TuneCUETP8M1_GEN_SIM.root'),
    nbPileupEvents = cms.PSet(
        averageNumber = cms.double(1.5)
    ),
    sequential = cms.untracked.bool(False),
    type = cms.string('poisson')
)

process.mix.cosmics = cms.SecSource("EmbeddedRootSource",
        nbPileupEvents = cms.PSet(
            averageNumber = cms.double(1.6625e-05) # 1.6625e-05
        ),
        seed = cms.int32(2345678),
        type = cms.string('poisson'),
	sequential = cms.untracked.bool(False),
        fileNames = cms.untracked.vstring('file:BIB_SIM.root')
)
	
process.mix.beamhalo_minus = cms.SecSource("EmbeddedRootSource",
        nbPileupEvents = cms.PSet(
            averageNumber = cms.double(0.00040503*10000.) # 0.00040503
        ),
        seed = cms.int32(3456789),
        type = cms.string('poisson'),
	sequential = cms.untracked.bool(False),
        fileNames = cms.untracked.vstring('file:BIB_SIM.root')
)

process.mix.beamhalo_plus = cms.SecSource("EmbeddedRootSource",
        nbPileupEvents = cms.PSet(
            averageNumber = cms.double(0.00040503*10000.) # 0.00040503
        ),
        seed = cms.int32(3456789),
        type = cms.string('poisson'),
	sequential = cms.untracked.bool(False),
        fileNames = cms.untracked.vstring('file:BIB_SIM.root')
)



process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxE = cms.double(10.01),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(9.99),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359),
        PartID = cms.vint32(12)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single Nu E 10')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# customisation of the process.

# Automatic addition of the customisation function from SimFbcm.SiPadDigitizer.aging
from SimFbcm.SiPadDigitizer.aging import no_aging 

#call to customisation function no_aging imported from SimFbcm.SiPadDigitizer.aging
process = no_aging(process)

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
