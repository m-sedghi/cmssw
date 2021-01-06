///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

FbcmGeometry::FbcmGeometry(): nStations_(0) {
	SetNumOfStations(0);
	theMap.clear();
	theMapToStations.clear();
  theMapToSiDies.clear();
  theMapToSiPads.clear();
}
FbcmGeometry::~FbcmGeometry() {}

const FbcmGeometry::DetTypeContainer& FbcmGeometry::detTypes() const { return theSiPadTypes; }

const FbcmGeometry::DetContainer& FbcmGeometry::detUnits() const { return theSiPadGeoms; }

const FbcmGeometry::DetContainer& FbcmGeometry::dets() const { return theDets; }

const FbcmGeometry::DetIdContainer& FbcmGeometry::detUnitIds() const { return theSiPadIds; }

const FbcmGeometry::DetIdContainer& FbcmGeometry::detIds() const { return theDetIds; }

const GeomDet* FbcmGeometry::idToDetUnit(DetId id) const { 
	 //return nullptr;
	return dynamic_cast<const GeomDet*>(idToDet(id)); 
}

const GeomDet* FbcmGeometry::idToDet(DetId id) const {
	//return dynamic_cast<const GeomDet*>(IdToSiPad(id));
	//return nullptr;
   mapIdToDet::const_iterator i = theMap.find(id); // it is Ok, due to the == operator defined in DetId.h
   // mapIdToDet::const_iterator i = theMap.find(id.rawId()); // it is also OK
   return (i != theMap.end()) ? i->second : nullptr;
}

const std::vector<FbcmStationGeom const*>& FbcmGeometry::Stations() const { return allStationGeoms; }

const std::vector<FbcmSiliconDieGeom const*>& FbcmGeometry::SiliconDies() const { return allSiliconDieGeoms; }

const std::vector<FbcmSiPadGeom const*>& FbcmGeometry::SiPads() const { return allSiPadGeoms; }

const FbcmSiPadGeom* FbcmGeometry::IdToSiPad(FbcmDetId id) const {
	mapIdToSiPad::const_iterator i = theMapToSiPads.find(id); 
	return (i != theMapToSiPads.end()) ? i->second : nullptr;
}

const FbcmSiliconDieGeom* FbcmGeometry::IdToSiliconDie(FbcmDetId id) const {
  	mapIdToSiDie::const_iterator i = theMapToSiDies.find(id.SiliconDieDetId()); 
	return (i != theMapToSiDies.end()) ? i->second : nullptr;
}

const FbcmStationGeom* FbcmGeometry::IdToStation(FbcmDetId id) const {
    mapIdToStation::const_iterator i = theMapToStations.find(id.StationDetId()); 
	return (i != theMapToStations.end()) ? i->second : nullptr;
}

void FbcmGeometry::add(FbcmSiPadGeom* SiPadGeom) {
  allSiPadGeoms.emplace_back(SiPadGeom);
  theSiPadGeoms.emplace_back(SiPadGeom);
  theSiPadIds.emplace_back(SiPadGeom->geographicalId());
  theDets.emplace_back(SiPadGeom);
  theDetIds.emplace_back(SiPadGeom->geographicalId());
  theSiPadTypes.emplace_back(&SiPadGeom->type());
  theMapToSiPads.insert(std::pair<DetId, FbcmSiPadGeom*>(SiPadGeom->id(), SiPadGeom));
  theMap.insert(std::pair<DetId, GeomDet*>(SiPadGeom->geographicalId(), SiPadGeom));
}

void FbcmGeometry::add(FbcmSiliconDieGeom* SiliconDieGeom) {
  allSiliconDieGeoms.emplace_back(SiliconDieGeom); 
  theDets.emplace_back(SiliconDieGeom);
  theDetIds.emplace_back(SiliconDieGeom->geographicalId());
  theSiPadTypes.emplace_back(&SiliconDieGeom->type()); // ?????
  theMapToSiDies.insert(std::pair<DetId, FbcmSiliconDieGeom*>(SiliconDieGeom->id().SiliconDieDetId(), SiliconDieGeom));
  //theMap.insert(std::pair<DetId, GeomDet*>(SiliconDieGeom->geographicalId(), SiliconDieGeom));
}

void FbcmGeometry::add(FbcmStationGeom* StationGeom) {
  allStationGeoms.emplace_back(StationGeom);
  theDets.emplace_back(StationGeom);
  theDetIds.emplace_back(StationGeom->geographicalId());
  theMapToStations.insert(std::pair<DetId, FbcmStationGeom*>(StationGeom->id().StationDetId(), StationGeom));
//  theMap.insert(std::pair<DetId, GeomDet*>(StationGeom->geographicalId(), StationGeom));
}
