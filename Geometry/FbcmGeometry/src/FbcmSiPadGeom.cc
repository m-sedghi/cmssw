///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#include "Geometry/FbcmGeometry/interface/FbcmSiPadGeom.h"


FbcmSiPadGeom::FbcmSiPadGeom(FbcmDetId id, const BoundPlane::BoundPlanePointer& bp, FbcmSiPadSpecs* rrs)
    : GeomDet(bp), id_(id), specs_(rrs) {
  setDetId(id);
}

FbcmSiPadGeom::~FbcmSiPadGeom() {
  delete specs_;  //Assume the SiPad owns it specs (specs are not shared)
}

const Topology& FbcmSiPadGeom::topology() const { return specs_->topology(); }

const FbcmSiPadTopology& FbcmSiPadGeom::SiPadTopology() const { return specs_->SiPadTopology(); }

const GeomDetType& FbcmSiPadGeom::type() const { return (*specs_); }
//SiPadSensorGeom->surface().position()