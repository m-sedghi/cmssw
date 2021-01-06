///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "Geometry/FbcmGeometry/interface/FbcmSiliconDieGeom.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadGeom.h"
#include <iostream>



FbcmSiliconDieGeom::FbcmSiliconDieGeom(FbcmDetId id, const ReferenceCountingPointer<BoundPlane>& plane, unsigned int nCols, unsigned int nRows) :
 GeomDet(plane),
 detId_(id),
 nCols_(nCols),
 nRows_(nRows) {
  setDetId(id);
}

FbcmSiliconDieGeom::~FbcmSiliconDieGeom() {}

FbcmDetId FbcmSiliconDieGeom::id() const { return detId_; }

bool FbcmSiliconDieGeom::operator==(const FbcmSiliconDieGeom& ch) const { return this->id() == ch.id(); }

void FbcmSiliconDieGeom::add(const FbcmSiPadGeom* SiPadGeom) { SiPads_.emplace_back(SiPadGeom); }

std::vector<const GeomDet*> FbcmSiliconDieGeom::components() const {
  return std::vector<const GeomDet*>(SiPads_.begin(), SiPads_.end());
}

const GeomDet* FbcmSiliconDieGeom::component(DetId id) const { return SiPad(FbcmDetId(id.rawId())); }

const std::vector<const FbcmSiPadGeom*>& FbcmSiliconDieGeom::SiPads() const { return SiPads_; }

int FbcmSiliconDieGeom::nSiPads() const { return SiPads_.size(); }

const FbcmSiPadGeom* FbcmSiliconDieGeom::SiPad(FbcmDetId id) const {
  if (id.SiliconDieDetId() != detId_)
    return nullptr;  // not in this SiPad! // Not in this SiliconDie
  return SiPad(id.SiPad());
}

const FbcmSiPadGeom* FbcmSiliconDieGeom::SiPad(unsigned int SiPadNo) const {
  for (auto SiPadItem : SiPads_) {
    if (SiPadItem->id().SiPad() == SiPadNo)
      return SiPadItem;
  }
  return nullptr;
}
