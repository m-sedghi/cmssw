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
    G4String DetectorName = touch->GetVolume(3)->GetName();
    G4String volumeName = touch->GetVolume(4)->GetName();

    if (SensorPadName != "FBCM_SensorPad") {
      edm::LogWarning("FBCMSim") << "FbcmSD::setDetUnitId -w- Sensor name is not FBCM_SensorPad ";
    }
    if (DetectorName != "FBCM") {
      edm::LogWarning("FBCMSim") << " FbcmSD::setDetUnitId -w- Detector name is not FBCM ";
    }
	
	// get the copyNumbers in the Geom. XML
    int SensorPadNo = touch->GetReplicaNumber(0); // 0-3
    int SensorRowNo = touch->GetReplicaNumber(1); // 0-1
	int SilcionDieNo = touch->GetReplicaNumber(1); //1-4
    int VolumeNo = touch->GetReplicaNumber(4); // 1-2

    // Detector ID definition
    // detId = ABC
    // A  = Phase2PixelEndcap_1 or _2,  1: +Z, 2: -Z
    // B = SilcionDieNo --> refers to the Station Number, 0 @45_deg, 
	//                      anticlockwise increment when looking toward to the IP
	//                      or eqivalenty clockwise when looking from the IP
    // Z  = SensorPadID,  0-7 (increament like a 2D array)

	// Assuing 4 Sensors per Row. 
	int SensorPadID=SensorRowNo*4+SensorPadNo; // (0-7): 2(rows)x4(cols)

    detId = 100 * VolumeNo + 10 * SilcionDieNo + SensorPadID;
  }
  return detId;
}

bool FbcmSD::checkHit(const G4Step*, BscG4Hit* hit) {
  // 50 micron are allowed between the exit
  // point of the current hit and the entry point of the new hit
  static const float tolerance2 = (float)(0.0025 * CLHEP::mm * CLHEP::mm);
  return ((hit->getExitLocalP() - getLocalEntryPoint()).mag2() < tolerance2);
}
