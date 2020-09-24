///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

FbcmGeometry::FbcmGeometry() {}
FbcmGeometry::~FbcmGeometry() {}

const FbcmGeometry::DetTypeContainer& FbcmGeometry::detTypes() const { return theSiPadTypes; }

const FbcmGeometry::DetContainer& FbcmGeometry::detUnits() const { return theSiPadGeoms; }

const FbcmGeometry::DetContainer& FbcmGeometry::dets() const { return theDets; }

const FbcmGeometry::DetIdContainer& FbcmGeometry::detUnitIds() const { return theSiPadIds; }

const FbcmGeometry::DetIdContainer& FbcmGeometry::detIds() const { return theDetIds; }

const GeomDet* FbcmGeometry::idToDetUnit(DetId id) const { return dynamic_cast<const GeomDet*>(idToDet(id)); }

const GeomDet* FbcmGeometry::idToDet(DetId id) const {
  mapIdToDet::const_iterator i = theMap.find(id); // it is Ok, due to the == operator defined in DetId.h
   // mapIdToDet::const_iterator i = theMap.find(id.rawId()); // it is also OK
  return (i != theMap.end()) ? i->second : nullptr;
}

const std::vector<FbcmStationGeom const*>& FbcmGeometry::Stations() const { return allStationGeoms; }

const std::vector<FbcmSiliconDieGeom const*>& FbcmGeometry::SiliconDies() const { return allSiliconDieGeoms; }

const std::vector<FbcmSiPadGeom const*>& FbcmGeometry::SiPads() const { return allSiPadGeoms; }

const FbcmSiPadGeom* FbcmGeometry::SiPad(FbcmDetId id) const {
  return dynamic_cast<const FbcmSiPadGeom*>(idToDetUnit(id));
}

const FbcmSiliconDieGeom* FbcmGeometry::SiliconDie(FbcmDetId id) const {
  return dynamic_cast<const FbcmSiliconDieGeom*>(idToDetUnit(id.SiliconDieId()));
}

const FbcmStationGeom* FbcmGeometry::Station(FbcmDetId id) const {
  return dynamic_cast<const FbcmStationGeom*>(idToDetUnit(id.StationId()));
}

void FbcmGeometry::add(FbcmSiPadGeom* SiPadGeom) {
  allSiPadGeoms.emplace_back(SiPadGeom);
  theSiPadGeoms.emplace_back(SiPadGeom);
  theSiPadIds.emplace_back(SiPadGeom->geographicalId());
  theDets.emplace_back(SiPadGeom);
  theDetIds.emplace_back(SiPadGeom->geographicalId());
  theSiPadTypes.emplace_back(&SiPadGeom->type());
  theMap.insert(std::pair<DetId, GeomDet*>(SiPadGeom->geographicalId(), SiPadGeom));
}

void FbcmGeometry::add(FbcmSiliconDieGeom* SiliconDieGeom) {
  allSiliconDieGeoms.emplace_back(SiliconDieGeom); 
  theDets.emplace_back(SiliconDieGeom);
  theDetIds.emplace_back(SiliconDieGeom->geographicalId());
  theSiPadTypes.emplace_back(&SiliconDieGeom->type()); // ?????
  theMap.insert(std::pair<DetId, GeomDet*>(SiliconDieGeom->geographicalId(), SiliconDieGeom));
}

void FbcmGeometry::add(FbcmStationGeom* StationGeom) {
  allStationGeoms.emplace_back(StationGeom);
  theDets.emplace_back(StationGeom);
  theDetIds.emplace_back(StationGeom->geographicalId());
  theMap.insert(std::pair<DetId, GeomDet*>(StationGeom->geographicalId(), StationGeom));
}
