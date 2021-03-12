///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "Geometry/FbcmGeometry/interface/FbcmStationGeom.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiliconDieGeom.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadGeom.h"
#include <iostream>



FbcmStationGeom::FbcmStationGeom(FbcmDetId id, const ReferenceCountingPointer<BoundPlane>& plane, unsigned int nDiesPerRing, unsigned int nRings) : GeomDet(plane),
	detId_(id),
	nDiesPerRing_(nDiesPerRing),
	nRings_(nRings) {
  setDetId(id);
}

FbcmStationGeom::~FbcmStationGeom() {}

FbcmDetId FbcmStationGeom::id() const { return detId_; }

bool FbcmStationGeom::operator==(const FbcmStationGeom& ch) const { return this->id() == ch.id(); }

void FbcmStationGeom::add(FbcmSiliconDieGeom* SiliconDieGeom) { SiliconDies_.emplace_back(SiliconDieGeom); }

std::vector<const GeomDet*> FbcmStationGeom::components() const {
  return std::vector<const GeomDet*>(SiliconDies_.begin(), SiliconDies_.end());
}

const GeomDet* FbcmStationGeom::component(DetId id) const { return SiliconDie(FbcmDetId(id.rawId())); }

const std::vector<const FbcmSiliconDieGeom*>& FbcmStationGeom::SiliconDies() const { return SiliconDies_; }

int FbcmStationGeom::nSiliconDies() const { return SiliconDies_.size(); }

const FbcmSiliconDieGeom* FbcmStationGeom::SiliconDie(FbcmDetId id) const {
  if (id.StationDetId() != detId_) // double check ?!!
    return nullptr;  // not in this SiliconDie! //  or not in this Staiotn?
  return SiliconDie(id.SiliconDie());
}

const FbcmSiliconDieGeom* FbcmStationGeom::SiliconDie(unsigned int SiDieNo) const {
  for (auto SiliconDieItem : SiliconDies_) 
  {
    if (SiliconDieItem->id().SiliconDie() == SiDieNo)
      return SiliconDieItem;
  }
  return nullptr;
}


// we may need to maintain this for a while
//void FbcmStationGeom::add(FbcmSiPadGeom* SiPadGeom) { SiPads_.emplace_back(SiPadGeom); }

const std::vector<const FbcmSiPadGeom*>& FbcmStationGeom::SiPads() const { return SiPads_; }

int FbcmStationGeom::nSiPads() const { return SiPads_.size(); }

const FbcmSiPadGeom* FbcmStationGeom::SiPad(FbcmDetId id) const {
  if (id.StationDetId() != detId_)
    return nullptr;  // not in this SiliconDie! //  or not in this Staiotn?
  return SiPad(id.SiPad());
}

const FbcmSiPadGeom* FbcmStationGeom::SiPad(unsigned int SiPadNo) const {
  for (auto SiPadItem : SiPads_) {
    if (SiPadItem->id().SiPad() == SiPadNo)
      return SiPadItem;
  }
  return nullptr;
}

/* ////??????
// it was assumed that the layers are paraller inside the chember
// However, the siliconDies are not paraller inside the Station. They are in the same palne
// Thus, this functions is besides the point and will not used in aour imlementation
float FbcmStationGeom::computeDeltaPhi(const LocalPoint& position, const LocalVector& direction) const {
  auto extrap = [](const LocalPoint& point, const LocalVector& dir, double extZ) -> LocalPoint 
  {
    double extX = point.x() + extZ * dir.x() / dir.z();
    double extY = point.y() + extZ * dir.y() / dir.z();
    return LocalPoint(extX, extY, extZ);
  } ;
  
  if (nSiliconDies() < 2) {
    return 0;
  }

  const float beginOfChamber = SiliconDie(1)->position().z();
  const float centerOfChamber = this->position().z();
  const float endOfChamber = SiliconDie(nSiliconDies())->position().z();

  LocalPoint projHigh =
      extrap(position, direction, (centerOfChamber < 0 ? -1.0 : 1.0) * (endOfChamber - centerOfChamber));
  LocalPoint projLow =
      extrap(position, direction, (centerOfChamber < 0 ? -1.0 : 1.0) * (beginOfChamber - centerOfChamber));
  auto globLow = toGlobal(projLow);
  auto globHigh = toGlobal(projHigh);
  return globHigh.phi() - globLow.phi();  //Geom::phi automatically normalizes to [-pi, pi]
}
 */