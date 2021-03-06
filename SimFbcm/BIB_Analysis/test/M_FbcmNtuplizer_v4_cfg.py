import FWCore.ParameterSet.Config as cms
#this is a copy of H_FbcmNtuplizer_v3_cfg, but intended to be run locally 
# Mohammad
process = cms.Process("ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")
import os
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('PU',
                  0,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "puvalue")

options.register ('inputDir',
		  '/afs/cern.ch/work/m/msedghi/public/BeamInducedBackgrdFbcm/bibDIGI_SelfMixed', 
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The input directory")				  
options.register ('filein',
		  'file:BIB_SIM_DIGI_mixed.root', 
                                 VarParsing.multiplicity.singleton,
                                 VarParsing.varType.string,
                  "The input fileName")				  
options.parseArguments()
BaseDir=os.path.dirname(options.inputDir)
outputDir=BaseDir+'/nTuplizerOutput/'
FilesList = ["file:" + options.inputDir + "/" + f for f in os.listdir(options.inputDir) if f[:4] == "Digi"]
outFName = options.inputDir.split('/')[-1] +'Ntuple'+'Out_pu{0}_.root'.format( options.PU )
print("output is saved in {0}".format( outFName ) )
outFileName=outputDir+outFName

#print(outFileName)
#print(BaseDir)
#os._exit(0)



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.Geometry.GeometryExtended2026D80_cff')
process.load('Geometry.FbcmGeometryBuilder.FbcmGeometry_cfi')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(FilesList)
)


process.FbcmNtuple = cms.EDAnalyzer('FbcmNtuplizer_v3',
                                    FbcmDigiTag = cms.InputTag("simFbcmDigis", "SiPad"),
                                    #FbcmDigiTag = cms.InputTag("simFbcmDigis"),
                                    TreeName = cms.string( 'PU{0}'.format(options.PU) )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outFileName),
                                   closeFileFast = cms.untracked.bool(True) )

process.p = cms.Path(process.FbcmNtuple)
