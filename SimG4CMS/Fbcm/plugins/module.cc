#include "SimG4Core/SensitiveDetector/interface/SensitiveDetectorPluginFactory.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "SimG4CMS/Fbcm/interface/FbcmSD.h"



typedef FbcmSD FBCMSensitiveDetector;
DEFINE_SENSITIVEDETECTOR(FBCMSensitiveDetector);

