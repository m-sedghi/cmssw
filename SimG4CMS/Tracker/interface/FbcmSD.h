#ifndef BRIL_FbcmSensitiveDetector_H
#define BRIL_FbcmSensitiveDetector_H

#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveTkDetector.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfJob.h"

#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/Records/interface/FbcmGeometryRecord.h"


#include "G4Step.hh"
#include "G4Track.hh"

#include <string>

class TrackInformation;
class SimTrackManager;
class TrackingSlaveSD;
class FrameRotation;
class UpdatablePSimHit;
class G4ProcessTypeEnumerator;
//class TrackerG4SimHitNumberingScheme;

class FbcmSD : public SensitiveTkDetector,
			public Observer<const BeginOfEvent *>,
			public Observer<const BeginOfTrack *>,
			public Observer<const BeginOfJob *> {
public:
  FbcmSD(const std::string &,
		  const edm::EventSetup &,
		  const SensitiveDetectorCatalog &,
		  edm::ParameterSet const &,
		  const SimTrackManager *);
  ~FbcmSD() override;
  bool ProcessHits(G4Step *, G4TouchableHistory *) override;
  uint32_t setDetUnitId(const G4Step *) override;
  void EndOfEvent(G4HCofThisEvent *) override;

  void fillHits(edm::PSimHitContainer &, const std::string &) override;
  void clearHits() override;

private:
  void createHit(const G4Step *);
  void sendHit();
  void updateHit(const G4Step *);
  bool newHit(const G4Step *);
  bool closeHit(const G4Step *);

protected:
  void update(const BeginOfEvent *) override;
  void update(const BeginOfTrack *) override;
  void update(const BeginOfJob *) override;

private:
  // data members initialised before run
  const SimTrackManager *theManager;
  std::unique_ptr<TrackingSlaveSD> _slaveSD;
  //std::unique_ptr<TrackingSlaveSD> slaveLowTof;
  //std::unique_ptr<TrackingSlaveSD> slaveHighTof;
  std::unique_ptr<FrameRotation> theRotation;
  std::unique_ptr<const G4ProcessTypeEnumerator> theG4ProcTypeEnumerator;
  //TrackerG4SimHitNumberingScheme *theNumberingScheme;  //does not own
  bool allowZeroEnergyLoss;
  bool printHits;
  bool neverAccumulate;
  double rMax2;  // tracker volume R^2
  double rMax;   // tracker volume R
  double rMin;
  double zMin;   // tracker volume Z
  double zMax;   // tracker volume Z
  
  //float theTofLimit;
  float energyCut;
  float energyHistoryCut;

  // run time cache
  UpdatablePSimHit *mySimHit;
  uint32_t lastId;
  int lastTrack;

  // -------- cache stuff for debugging and printout ---------
  Local3DPoint globalEntryPoint;
  Local3DPoint globalExitPoint;
  const G4VPhysicalVolume *oldVolume;
  float px, py, pz;
  int eventno;
  std::string pname;
  //----------------------------------------------------------
  
  edm::ESHandle<FbcmGeometry> FbcmGeom;
} ;

#endif
