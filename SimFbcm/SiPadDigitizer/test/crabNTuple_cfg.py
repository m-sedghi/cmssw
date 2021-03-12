from WMCore.Configuration import Configuration
import os,sys
config = Configuration()

reqNamedFromArg = [ arg for arg in sys.argv if arg.startswith( 'General.requestName=' ) ][0].split( '=' )[-1]
puFromArg = reqNamedFromArg[ reqNamedFromArg.find('PU')+2:]
generationInfo = {'0p5':[0.5 , '/NuGun/hbakhshi-FBCMNuGunPU0p5-bee0f6cca29d2ad52c10723372b5e551/USER'] ,
                  '1' : [1.0 , '/NuGun/hbakhshi-FBCMNuGunPU1-ecd2a523d2f78c811432ede5279aea28/USER'] ,
                  '1p5' : [1.5 , '/NuGun/hbakhshi-FBCMNuGunPU1p5-5574e5de471c21e9f6447b511c3f5447/USER' ] ,
                  '10' : [10 , '/NuGun/hbakhshi-FBCMNuGunPU10-26b3d8a1ea4fdb6b830c981c46852ea7/USER' ] ,
                  '50' : [50 , '/NuGun/hbakhshi-FBCMNuGunPU50-1b900cf6c3e80f0c6610b3f5adb26382/USER' ] , 
                  '100' : [100 , '/NuGun/hbakhshi-FBCMNuGunPU100-f65b27cb4f63546e3d4064b48622f535/USER'],
                  '140' : [140 , '/NuGun/hbakhshi-FBCMNuGunPU140-d33d56216e56dc2937a9e0b9b09425e6/USER' ] ,
                  '200' : [200 , '/NuGun/hbakhshi-FBCMNuGunPU200-79ef82d501ae7989114bc727a8340bda/USER' ] }
config.section_('General')
config.General.requestName = ''
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'H_FbcmNtuplizer_v3_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.sendPythonFolder	 = True
config.JobType.numCores = 1
config.JobType.pyCfgParams = ["PU={0}".format(puFromArg)]

config.section_('Data')
config.Data.inputDataset = generationInfo[puFromArg][1]
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/'
config.Data.inputDBS = 'phys03'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
config.Site.whitelist = ["T2_CH_CERN"]
