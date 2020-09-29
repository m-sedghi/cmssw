#include "SimG4CMS/Tracker/interface/FbcmSD.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h" //
#include "DataFormats/GeometryVector/interface/LocalVector.h" //

#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/G4TrackToParticleID.h"
#include "SimG4Core/Physics/interface/G4ProcessTypeEnumerator.h"

#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/SimHitMaker/interface/TrackingSlaveSD.h"

#include "SimG4Core/Notification/interface/TrackInformation.h"

#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Track.hh"
#include "G4SDManager.hh" //
#include "G4VProcess.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4VProcess.hh"

#include <string>
#include <vector>
#include <iostream>

#include "CLHEP/Units/GlobalSystemOfUnits.h"

FbcmSD::FbcmSD(const std::string& name,
                 const edm::EventSetup& es,
                 const SensitiveDetectorCatalog& clg,
                 edm::ParameterSet const& p,
                 const SimTrackManager* manager)
    : TimingSD(name, es, clg, p, manager) {
  edm::ParameterSet m_TrackerSD = p.getParameter<edm::ParameterSet>("FbcmSD");
  energyCut = m_TrackerSD.getParameter<double>("EnergyThresholdForPersistencyInGeV") * GeV;  //default must be 0.5 (?)
  energyHistoryCut =
      m_TrackerSD.getParameter<double>("EnergyThresholdForHistoryInGeV") * GeV;  //default must be 0.05 (?)

  setCuts(energyCut, energyHistoryCut);
}

FbcmSD::~FbcmSD() {}

uint32_t FbcmSD::setDetUnitId(const G4Step* aStep) {
  uint32_t detId = 0;
	
	
	
  //Find number of levels
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
//  int level = (touch) ? ((touch->GetHistoryDepth()) + 1) : 0;
  int level = 0;
    if (touch)
       level = ((touch->GetHistoryDepth()) + 1);

  //Get name and copy numbers
  if (level > 1) {
    G4String SensorPadName = touch->GetVolume(0)->GetName();
    G4String SensorRowName = touch->GetVolume(1)->GetName();
	G4String SiliconDieName = touch->GetVolume(2)->GetName();
	G4String StationName = touch->GetVolume(3)->GetName();
    G4String DetectorName = touch->GetVolume(4)->GetName();
    G4String volumeName = touch->GetVolume(5)->GetName();

	//edm::LogWarning("FBCM-w-levelNames") << SensorPadName << ", " << SensorRowName << ", " << SiliconDieName << ", " << StationName << ", " << DetectorName << ", " <<volumeName << "\n";
	//edm::LogInfo("FBCM-I-levelNames") << SensorPadName << ", " << SensorRowName << ", " << SiliconDieName << ", " << StationName << ", " << DetectorName << ", " <<volumeName << "\n";
	
	
    if (SensorPadName != "FBCM_SensorPad") {
      edm::LogWarning("FBCMSim") << "FbcmSD::setDetUnitId -w- SensorPadName is not FBCM_SensorPad ";
    }
	if (SensorRowName != "FBCM_SensorRow") {
      edm::LogWarning("FBCMSim") << "FbcmSD::setDetUnitId -w- SensorRowName is not FBCM_SensorRow ";
    }
	if (SiliconDieName != "FBCM_SiliconDie") {
      edm::LogWarning("FBCMSim") << "FbcmSD::setDetUnitId -w- SiliconDieName is not FBCM_SiliconDie ";
    }
	if (StationName != "FBCM_Station") {
      edm::LogWarning("FBCMSim") << "FbcmSD::setDetUnitId -w- StationName is not FBCM_Station ";
    }
    if (DetectorName != "FBCM") {
      edm::LogWarning("FBCMSim") << " FbcmSD::setDetUnitId -w- DetectorName is not FBCM ";
    }
	if (volumeName != "Phase2PixelEndcap") {
      edm::LogWarning("FBCMSim") << " FbcmSD::setDetUnitId -w- volumeName is not Phase2PixelEndcap ";
    }
	
	// get the copyNumbers in the Geom. XML
    int SensorPadNo = touch->GetReplicaNumber(0); // 0-1
    int SensorRowNo = touch->GetReplicaNumber(1); // 0-3
	int SilcionDieNo = touch->GetReplicaNumber(2); //0
	int StationNo = touch->GetReplicaNumber(3); //0-3
    int VolumeNo = touch->GetReplicaNumber(5); // 1-2

    // Detector ID definition
    // detId = ABC
    // A  = Phase2PixelEndcap_1 or _2,  1: +Z, 2: -Z
    // B = StationNo --> refers to the Station Number, 0 @45_deg, 
	//                      anticlockwise increment when looking toward to the IP
	//                      or eqivalenty clockwise when looking from the IP
	// SilcionDieNo=0
    // Z  = SensorPadID,  0-7 (increament like a 2D array)

	// Assuming 2 Sensors per Row. and 4 Rows
	int SensorPadID=SensorRowNo*2+SensorPadNo; // (0-7): 4(rows)x2(cols)
	//each Station comprises one SilcionDie
	
	//New FbcmsDetID assignment:
	FbcmDetId fbcmdet1(VolumeNo,StationNo,SilcionDieNo,SensorPadID);

    detId = fbcmdet1.rawId();
	//edm::LogVerbatim("FwkReport") << "*-FbcmG4Sim: A new G4SimHit occurred at: " << fbcmdet1 << "\n";
	//std::cout << "--**-- FbcmG4Sim: A new G4SimHit occurred at: " << fbcmdet1 ;
  }
  return detId;
}

bool FbcmSD::checkHit(const G4Step*, BscG4Hit* hit) {
  // 50 micron are allowed between the exit
  // point of the current hit and the entry point of the new hit
  static const float tolerance2 = (float)(0.0025 * CLHEP::mm * CLHEP::mm);
  return ((hit->getExitLocalP() - getLocalEntryPoint()).mag2() < tolerance2);
}
