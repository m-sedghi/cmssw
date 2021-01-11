///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "SimG4CMS/Fbcm/interface/FbcmSD.h"

//#define FAKEFRAMEROTATION


FbcmSD::FbcmSD(const std::string& name,
				 const edm::EventSetup& es,
				 const SensitiveDetectorCatalog& clg,
				 edm::ParameterSet const& p,
				 const SimTrackManager* manager)
    : SensitiveTkDetector(name, es, clg, p),
      theManager(manager),
      mySimHit(nullptr),
      lastId(0),
      lastTrack(0),
      oldVolume(nullptr),
      px(0.0f),
      py(0.0f),
      pz(0.0f),
      eventno(0),
      pname("") {
  edm::ParameterSet FbcmSD_Param = p.getParameter<edm::ParameterSet>("FbcmSD");
  allowZeroEnergyLoss = FbcmSD_Param.getParameter<bool>("ZeroEnergyLoss");
  neverAccumulate = FbcmSD_Param.getParameter<bool>("NeverAccumulate");
  printHits = FbcmSD_Param.getParameter<bool>("PrintHits");
  //theTofLimit = FbcmSD_Param.getParameter<double>("ElectronicSigmaInNanoSeconds") * 3 * CLHEP::ns;  // 3 sigma
  energyCut = FbcmSD_Param.getParameter<double>("EnergyThresholdForPersistencyInGeV") * CLHEP::GeV;  //default must be 0.5
  energyHistoryCut = FbcmSD_Param.getParameter<double>("EnergyThresholdForHistoryInGeV") * CLHEP::GeV;  //default must be 0.05
  zMin = FbcmSD_Param.getParameter<double>("zMin") * CLHEP::mm;  
  zMax = FbcmSD_Param.getParameter<double>("zMax") * CLHEP::mm;  
  rMin = FbcmSD_Param.getParameter<double>("rMin") * CLHEP::mm;  
  rMax = FbcmSD_Param.getParameter<double>("rMax") * CLHEP::mm;  
  
  rMax2 = rMax * rMax;

  // No Rotation given in input, automagically choose one based upon the name
  std::string rotType;
  theRotation.reset(new TrackerFrameRotation());
  rotType = "TrackerFrameRotation";

#ifdef FAKEFRAMEROTATION
  theRotation.reset(new FakeFrameRotation());
  rotType = "FakeFrameRotation";
#endif

  edm::LogInfo("TrackerSimInfo") << " FbcmSD: "
                                 << " Criteria for Saving Tracker SimTracks: \n"
                                 << " History: " << energyHistoryCut << " MeV; Persistency: " << energyCut
                              //   << " MeV;  TofLimit: " << theTofLimit << " ns"
                                 << "\n FrameRotation type " << rotType << " rMax(cm)= " << rMax / CLHEP::cm
                                 << " zMax(cm)= " << zMax / CLHEP::cm
                                 << " allowZeroEnergyLoss: " << allowZeroEnergyLoss
                                 << " neverAccumulate: " << neverAccumulate << " printHits: " << printHits;

  _slaveSD.reset(new TrackingSlaveSD(name));
  ///---- no need for the following 3 lines: ------------
  //std::vector<std::string> temp;
  //temp.push_back(_slaveSD.get()->name());
  //setNames(temp);
  ///-----------------------------------------------------

  theG4ProcTypeEnumerator.reset(new G4ProcessTypeEnumerator);
   
  es.get<FbcmGeometryRecord>().get(FbcmGeom);
    
}

FbcmSD::~FbcmSD() {}

bool FbcmSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  LogDebug("FbcmSimDebug") << " Entering a new Step " << aStep->GetTotalEnergyDeposit() << " "
                              << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();

  if (aStep->GetTotalEnergyDeposit() > 0. || allowZeroEnergyLoss) {
    if (!mySimHit) {
      createHit(aStep);
    } else if (neverAccumulate || newHit(aStep)) {
      sendHit();
      createHit(aStep);
    } else {
      updateHit(aStep);
    }
    return true;
  }
  return false;
}

uint32_t FbcmSD::setDetUnitId(const G4Step* aStep) {
  uint32_t detId = 0;
	
	//std::cout << "FbcmGeom in setDetUnitId, NumOfSiPads: " << FbcmGeom->SiPads().size() << "\n";
	
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
	//std::cout << "hit touched at: " << SensorPadName << ", " << SensorRowName << ", " << SiliconDieName << ", " << StationName << ", " << DetectorName << ", " <<volumeName << "\n";
	
	
	// get the copyNumbers in the Geom. XML
    int SensorPadNo = touch->GetReplicaNumber(0); // 
    int SensorColNo = touch->GetReplicaNumber(1); // 
	int SilcionDieNo = touch->GetReplicaNumber(2); //
	int StationNo = touch->GetReplicaNumber(3); //
    int VolumeNo = touch->GetReplicaNumber(4); // 

	//std:: cout << "SensorPadNo:" <<SensorPadNo << ", SensorColNo:" <<SensorColNo << ", SilcionDieNo: " << SilcionDieNo<< ", StationNo:" << StationNo<< ", FbcmNo:" << VolumeNo<< "\n" ;
	//New FbcmsDetID assignment:
	FbcmDetId TMPfbcmDet(VolumeNo,StationNo,SilcionDieNo,0);
	const FbcmSiliconDieGeom * DieGeomPtr=FbcmGeom->IdToSiliconDie(TMPfbcmDet);
	int nRows=0;
	if LIKELY(DieGeomPtr!=nullptr) {
		nRows=DieGeomPtr->NumOfRows() ;
		  //std::cout << "nRows: " << nRows <<"\tnCols: " << DieGeomPtr->NumOfCols() << ", ColNo: " << SensorColNo << ", RowNo: " << SensorPadNo << "\n"; 
		//std::cout << ",\tfor: " << DieGeomPtr->id();
	}
	else {edm::LogError("FbcmSD") << "Illegal Copy number, SiDie or Station \n"; }

	int SensorPadID=SensorColNo*nRows+SensorPadNo; 
		
	FbcmDetId fbcmdet1(VolumeNo,StationNo,SilcionDieNo,SensorPadID);
	
    detId = fbcmdet1.rawId();
	//edm::LogVerbatim("FwkReport") << "*-FbcmG4Sim: A new G4SimHit occurred at: " << fbcmdet1 << "\n";
	std::cout << "A new G4SimHit: " << fbcmdet1 ;
  }
  return detId;
}

void FbcmSD::update(const BeginOfTrack* bot) {
  const G4Track* gTrack = (*bot)();

#ifdef DUMPPROCESSES
  if (gTrack->GetCreatorProcess()) {
    edm::LogVerbatim("TrackerSimInfo") << " -> PROCESS CREATOR : " << gTrack->GetCreatorProcess()->GetProcessName();
  } else {
    edm::LogVerbatim("TrackerSimInfo") << " -> No Creator process";
  }
#endif

  //
  //Position
  //
  const G4ThreeVector& pos = gTrack->GetPosition();
		LogDebug("FbcmSimDebug") << " update(..) of " << gTrack->GetDefinition()->GetParticleName()
		                      << " trackID= " << gTrack->GetTrackID() << " E(MeV)= " << gTrack->GetKineticEnergy()
                              << " Ecut= " << energyCut << " R(mm)= " << pos.perp() << " Z(mm)= " << pos.z();

  //
  // Check if in Tracker Volume
  //
  if (pos.x() * pos.x() + pos.y() * pos.y() < rMax2 && std::abs(pos.z()) < zMax ) {
    	
	//
    // inside the Tracker
    //
    TrackInformation* info = nullptr;
    if (gTrack->GetKineticEnergy() > energyCut) {
      info = cmsTrackInformation(gTrack);
      info->storeTrack(true);
    }
    //
    // Save History?
    //
    if (gTrack->GetKineticEnergy() > energyHistoryCut) {
      if (!info) {
        info = cmsTrackInformation(gTrack);
      }
      info->putInHistory();
      LogDebug("FbcmSimDebug") << " Track inside the tracker selected for HISTORY"
                                  << " Track ID= " << gTrack->GetTrackID();
    }
  }
  //else
	  //std::cout << "Vaowww!\n";
}

void FbcmSD::sendHit() {
  if (mySimHit == nullptr)
    return;
  if (printHits) {
    TkSimHitPrinter thePrinter("FbcmHitPositionG4Sim.dat"); //TkHitPositionOSCAR.dat
    thePrinter.startNewSimHit(GetName(),
                              oldVolume->GetLogicalVolume()->GetName(),
                              mySimHit->detUnitId(),
                              mySimHit->trackId(),
                              lastTrack,
                              eventno);
    thePrinter.printLocal(mySimHit->entryPoint(), mySimHit->exitPoint());
    thePrinter.printGlobal(globalEntryPoint, globalExitPoint);
    thePrinter.printHitData(pname, mySimHit->pabs(), mySimHit->energyLoss(), mySimHit->timeOfFlight());
    thePrinter.printGlobalMomentum(px, py, pz);
    LogDebug("FbcmSimDebug") << " Storing PSimHit: " << mySimHit->detUnitId() << " " << mySimHit->trackId() << " "
                                << mySimHit->energyLoss() << " " << mySimHit->entryPoint() << " "
                                << mySimHit->exitPoint();
  }


  _slaveSD.get()->processHits(*mySimHit);
  
  //
  // clean up
  delete mySimHit;
  mySimHit = nullptr;
  lastTrack = 0;
  lastId = 0;
}

void FbcmSD::createHit(const G4Step* aStep) {
  // VI: previous hit should be already deleted
  //     in past here was a check if a hit is inside a sensitive detector,
  //     this is not needed, because call to senstive detector happens
  //     only inside the volume
  const G4Track* theTrack = aStep->GetTrack();
  Local3DPoint theExitPoint = theRotation.get()->transformPoint(LocalPostStepPosition(aStep));
  Local3DPoint theEntryPoint;
  //
  //  Check particle type - for gamma and neutral hadrons energy deposition
  //  should be local (VI)
  //
  if (0.0 == theTrack->GetDefinition()->GetPDGCharge()) {
    theEntryPoint = theExitPoint;
  } else {
    theEntryPoint = theRotation.get()->transformPoint(LocalPreStepPosition(aStep));
  }

  //
  //	This allows to send he skipEvent if it is outside!
  //
  const G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  float thePabs = preStepPoint->GetMomentum().mag() / GeV;
  float theTof = preStepPoint->GetGlobalTime() / nanosecond;
  float theEnergyLoss = aStep->GetTotalEnergyDeposit() / GeV;
  int theParticleType = G4TrackToParticleID::particleID(theTrack);
  uint32_t theDetUnitId = setDetUnitId(aStep);
  int theTrackID = theTrack->GetTrackID();
  if (theDetUnitId == 0) {
    edm::LogWarning("FbcmSD::createHit") << " theDetUnitId is not valid for " << GetName();
    throw cms::Exception("FbcmSD::createHit")
        << "cannot get theDetUnitId for G4Track " << theTrackID;
  }

  // To whom assign the Hit?
  // First iteration: if the track is to be stored, use the current number;
  // otherwise, get to the mother
  unsigned int theTrackIDInsideTheSimHit = theTrackID;

  const TrackInformation* temp = cmsTrackInformation(theTrack);
  if (!temp->storeTrack()) {
    // Go to the mother!
    theTrackIDInsideTheSimHit = theTrack->GetParentID();
    LogDebug("FbcmSimDebug") << " FbcmSD::createHit(): setting the TrackID from "
                                << theTrackIDInsideTheSimHit << " to the mother one " << theTrackIDInsideTheSimHit
                                << " " << theEnergyLoss;
  } else {
    LogDebug("FbcmSimDebug") << " FbcmSD:createHit(): leaving the current TrackID "
                                << theTrackIDInsideTheSimHit;
  }

  const G4ThreeVector& gmd = preStepPoint->GetMomentumDirection();
  // convert it to local frame
  G4ThreeVector lmd =
      ((G4TouchableHistory*)(preStepPoint->GetTouchable()))->GetHistory()->GetTopTransform().TransformAxis(gmd);
  Local3DPoint lnmd = theRotation.get()->transformPoint(ConvertToLocal3DPoint(lmd));
  float theThetaAtEntry = lnmd.theta();
  float thePhiAtEntry = lnmd.phi();

  mySimHit = new UpdatablePSimHit(theEntryPoint,
                                  theExitPoint,
                                  thePabs,
                                  theTof,
                                  theEnergyLoss,
                                  theParticleType,
                                  theDetUnitId,
                                  theTrackIDInsideTheSimHit,
                                  theThetaAtEntry,
                                  thePhiAtEntry,
                                  theG4ProcTypeEnumerator.get()->processId(theTrack->GetCreatorProcess()));
  lastId = theDetUnitId;
  lastTrack = theTrackID;

  // only for debugging
  if (printHits) {
    // point on Geant4 unit (mm)
    globalEntryPoint = ConvertToLocal3DPoint(preStepPoint->GetPosition());
    globalExitPoint = ConvertToLocal3DPoint(aStep->GetPostStepPoint()->GetPosition());
    // in CMS unit (GeV)
    px = preStepPoint->GetMomentum().x() / CLHEP::GeV;
    py = preStepPoint->GetMomentum().y() / CLHEP::GeV;
    pz = preStepPoint->GetMomentum().z() / CLHEP::GeV;
    oldVolume = preStepPoint->GetPhysicalVolume();
    pname = theTrack->GetDefinition()->GetParticleName();
    LogDebug("FbcmSimDebug") << " Created PSimHit: " << pname << " " << mySimHit->detUnitId() << " "
                                << mySimHit->trackId() << " " << theTrackID
                                << " p= " << aStep->GetPreStepPoint()->GetMomentum().mag() << " "
                                << mySimHit->energyLoss() << " " << mySimHit->entryPoint() << " "
                                << mySimHit->exitPoint();
  }
}

void FbcmSD::updateHit(const G4Step* aStep) {
  // VI: in past here was a check if a hit is inside a sensitive detector,
  //     this is not needed, because call to senstive detector happens
  //     only inside the volume
  Local3DPoint theExitPoint = theRotation.get()->transformPoint(LocalPostStepPosition(aStep));
  float theEnergyLoss = aStep->GetTotalEnergyDeposit() / GeV;
  mySimHit->setExitPoint(theExitPoint);
  mySimHit->addEnergyLoss(theEnergyLoss);
  if (printHits) {
    globalExitPoint = ConvertToLocal3DPoint(aStep->GetPostStepPoint()->GetPosition());
    LogDebug("FbcmSimDebug") << " updateHit: for " << aStep->GetTrack()->GetDefinition()->GetParticleName()
                                << " trackID= " << aStep->GetTrack()->GetTrackID() << " deltaEloss= " << theEnergyLoss
                                << "\n Updated PSimHit: " << mySimHit->detUnitId() << " " << mySimHit->trackId() << " "
                                << mySimHit->energyLoss() << " " << mySimHit->entryPoint() << " "
                                << mySimHit->exitPoint();
  }
}

bool FbcmSD::newHit(const G4Step* aStep) {
  const G4Track* theTrack = aStep->GetTrack();

  // for neutral particles do not merge hits (V.I.)
  if (0.0 == theTrack->GetDefinition()->GetPDGCharge())
    return true;

  uint32_t theDetUnitId = setDetUnitId(aStep);
  int theTrackID = theTrack->GetTrackID();

  LogDebug("FbcmSimDebug") << "newHit: OLD(detID,trID) = (" << lastId << "," << lastTrack << "), NEW = ("
                              << theDetUnitId << "," << theTrackID << ") Step length(mm)= " << aStep->GetStepLength()
                              << " Edep= " << aStep->GetTotalEnergyDeposit()
                              << " p= " << aStep->GetPreStepPoint()->GetMomentum().mag();
  return ((theTrackID == lastTrack) && (lastId == theDetUnitId) && closeHit(aStep)) ? false : true;
}

bool FbcmSD::closeHit(const G4Step* aStep) {
  const float tolerance2 = 0.0025f;  // (0.5 mm)^2 are allowed between entry and exit
  Local3DPoint theEntryPoint = theRotation.get()->transformPoint(LocalPreStepPosition(aStep));
  LogDebug("FbcmSimDebug") << " closeHit: distance = " << (mySimHit->exitPoint() - theEntryPoint).mag();

  return ((mySimHit->exitPoint() - theEntryPoint).mag2() < tolerance2) ? true : false;
}

void FbcmSD::EndOfEvent(G4HCofThisEvent*) {
  LogDebug("FbcmSimDebug") << " Saving the last hit in a ROU " << GetName();
  if (mySimHit != nullptr)
    sendHit();
}

void FbcmSD::update(const BeginOfEvent* i) {
  clearHits();
  eventno = (*i)()->GetEventID();
  delete mySimHit;
  mySimHit = nullptr;
}

void FbcmSD::update(const BeginOfJob* i) {
  edm::ESHandle<GeometricDet> pDD;
  const edm::EventSetup* es = (*i)();
  es->get<IdealGeometryRecord>().get(pDD);
	// Here I can load the FbcmGeometry rather than in the initializer
}

void FbcmSD::clearHits() {
  _slaveSD.get()->Initialize();
  
}

void FbcmSD::fillHits(edm::PSimHitContainer& cc, const std::string& hname) {
  if (_slaveSD.get()->name() == hname)
    cc = _slaveSD.get()->hits();
}
