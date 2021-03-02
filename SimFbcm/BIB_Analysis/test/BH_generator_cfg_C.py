###########################################
#
# BH_generator.py
#
# Test script for MIB generation (generate 1000 B1 events
# using either FLUKA or MARS inputs)
#
# To use is just do
#
# cmsRun BH_generator.py
#
# SV: 20/12/2010
#
##########################################

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os
from Configuration.StandardSequences.Eras import eras

# In the line below 'analysis' is an instance of VarParsing object 
options = VarParsing ('analysis')

# Here we have defined our own two VarParsing options 
# add a list of strings for events to process
# this should be smaller than the total number of events in the file or otherwise get segfault when there are no more events in file
nevents = 1600
options.register ('nEvents',
                                 nevents,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of events to generate: {}".format(nevents))
options.register ('nThreads',
                                 1,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The number of threads to use: 1")
options.register ('jobId',
                                 0,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The job Id: 0")
options.register ('outputDirectory',
		  'file:/afs/cern.ch/work/m/msedghi/public/bibGeneratorOutput', 
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The output directory")
				  
options.register ('inFile',
		  'file:/afs/cern.ch/work/g/gauzinge/public/BeamHalo/run0001_hilumi_ir5_exp_SCO001_fort.30', 
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The input fileName")				  
options.register ('FileNo',
                                 0,
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.int,
                  "The Number: 0")

options.parseArguments()
#options.inputFiles = options.inFile
#print("Input file: %s" % (options.inputFiles))

# #specify input
#inputPathBeamHalo = "/afs/cern.ch/work/g/gauzinge/public/BeamHalo"
inputPathBeamGasCarbon = "/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_carbon"
# inputPathBeamGasHydrogen = "/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_hydrogen"
# inputPathBeamGasOxygen = "/afs/cern.ch/work/g/gauzinge/public/BeamGas/cms_conditioned_oxygen"

# #options.inputFiles= [inputPath + "/" + f for f in os.listdir(inputPath) if f[:3] == "run"][0]
#beamHaloInputFiles = [inputPathBeamHalo + "/" + f for f in os.listdir(inputPathBeamHalo) if f[:3] == "run"]
beamGasCarbonInputFiles = [inputPathBeamGasCarbon + "/" + f for f in os.listdir(inputPathBeamGasCarbon) if f[:3] == "run"]
# beamGasHydrogenInputFiles = [inputPathBeamGasHydrogen + "/" + f for f in os.listdir(inputPathBeamGasHydrogen) if f[:3] == "run"]
# beamGasOxygenInputFiles = [inputPathBeamGasOxygen + "/" + f for f in os.listdir(inputPathBeamGasOxygen) if f[:3] == "run"]

options.inputFiles = beamGasCarbonInputFiles[options.FileNo]
print("Input file: %s" % (options.inputFiles))
#options.nEvents = sum([len(set([int(line.split()[0]) for line in open(f) if line.split()[0].isdigit()])) for f in options.inputFiles])
print(options.nEvents)
#print(len(options.inputFiles))

# print("Number of files for BeamHalo: {}".format(len(beamHaloInputFiles)))
# print("Number of files for BeamGasCarbon: {}".format(len(beamGasCarbonInputFiles)))
# print("Number of files for BeamGasHydrogen: {}".format(len(beamGasHydrogenInputFiles)))
# print("Number of files for BeamGasOxygen: {}".format(len(beamGasOxygenInputFiles)))

#count number of events in all input files
#print("Number of events in BeamHalo files: {}".format(sum([len(set([int(line.split()[0]) for line in open(f) if line.split()[0].isdigit()])) for f in beamHaloInputFiles])))
#print("Number of events in BeamGasCarbon files: {}".format(sum([len(set([int(line.split()[0]) for line in open(f) if line.split()[0].isdigit()])) for f in beamGasCarbonInputFiles])))
#print("Number of events in BeamGasHydrogen files: {}".format(sum([len(set([int(line.split()[0]) for line in open(f) if line.split()[0].isdigit()])) for f in beamGasHydrogenInputFiles])))
#print("Number of events in BeamGasOxygen files: {}".format(sum([len(set([int(line.split()[0]) for line in open(f) if line.split()[0].isdigit()])) for f in beamGasOxygenInputFiles])))



#specify output name
BeamHaloOutputFile = options.outputDirectory+'/BeamHalo/GenBeamHalo_'+str(options.FileNo)+'.root'
BeamGasCarbonOutputFile = options.outputDirectory+'/BeamGasCarbon/GenBeamGasCarbon_'+str(options.FileNo)+'.root'
BeamGasHydrogenOutputFile = options.outputDirectory+'/BeamGasHydrogen/GenBeamGasHydrogen_'+str(options.FileNo)+'.root'
BeamGasOxygenOutputFile = options.outputDirectory+'/BeamGasOxygen/GenBeamGasOxygen_'+str(options.FileNo)+'.root'

options.outputFile = BeamGasCarbonOutputFile
print("Output File: %s" % (options.outputFile))

process = cms.Process('GEN', eras.Phase2) # Generation only

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D80_cff')
process.load('Configuration.Geometry.GeometryExtended2026D80Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
# process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.Reconstruction_cff')


#add message logger
"""
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.categories.append('BH_generator')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)
"""

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.168.2.1 $'),
    annotation = cms.untracked.string('Test script for MIB production'),
    name = cms.untracked.string('PyReleaseValidation'),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nEvents)
)

# Input source
process.source = cms.Source("EmptySource")


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string(options.outputFile),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Global tag clearly depends on the release you are working on
# To get the correct GT, look at this page:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?redirectedfrom=CMS.SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
#

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


# Here we choose the input, look into MIB_generator_cfi for more infos about the inputs
 
# process.load('BRIL_BIBGenerator.GeneratorInterface.BeamHaloGenerator.MIB_generator_cff')
process.load('GeneratorInterface.BeamHaloGenerator.MIB_generator_cff')

process.generator       = process.FLUKA_generator.clone()  # FLUKA
# process.generator       = process.MARS_generator   # MARS
# process.generator.InputFile = cms.string(options.inputFiles)
process.generator.FlukaFiles = cms.vstring(options.inputFiles)
#process.generator.InputFilej = cms.vstring(" ", "/afs/cern.ch/work/g/gauzinge/public/BeamHalo/run0001_hilumi_ir5_exp_SCO001_fort.30")

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.endjob_step     = cms.Path(process.endOfProcess)
process.out_step        = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.endjob_step,process.out_step)

# special treatment in case of production filter sequence  
for path in process.paths:
    getattr(process,path)._seq = process.generator*getattr(process,path)._seq
