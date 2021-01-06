#include "SimG4CMS/Tracker/interface/TkAccumulatingSensitiveDetector.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveDetectorPluginFactory.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "SimG4CMS/Tracker/interface/FbcmSD.h"



typedef FbcmSD FBCMSensitiveDetector;
DEFINE_SENSITIVEDETECTOR(FBCMSensitiveDetector);
DEFINE_SENSITIVEDETECTOR(TkAccumulatingSensitiveDetector);
