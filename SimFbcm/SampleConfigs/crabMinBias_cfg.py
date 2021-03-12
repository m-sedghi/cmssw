from WMCore.Configuration import Configuration
import os,sys
config = Configuration()

config.section_('General')
config.General.requestName = 'FBCMMinBiasGeneration'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'MinBias_14TeV_pythia8_TuneCUETP8M1_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 3000
config.JobType.sendPythonFolder	 = True
config.JobType.numCores = 8
config.JobType.maxMemoryMB = 10000
config.JobType.maxJobRuntimeMin = 5000

config.section_('Data')
config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 2000
config.Data.totalUnits = 1000000 #config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'FBCMMinBias'
config.Data.outLFNDirBase = '/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
