///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#ifndef Geometry_FbcmSiPadGeom_H
#define Geometry_FbcmSiPadGeom_H

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h" 
#include "Geometry/FbcmGeometry/interface/FbcmSiPadTopology.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadSpecs.h"


class FbcmSiPadGeom : public GeomDet {
public:
  FbcmSiPadGeom(FbcmDetId id, const BoundPlane::BoundPlanePointer& bp, FbcmSiPadSpecs* rrs);
  ~FbcmSiPadGeom() override;

  const FbcmSiPadSpecs* specs() const { return specs_; }
  FbcmDetId id() const { return id_; }

  const Topology& topology() const override;
  const FbcmSiPadTopology& SiPadTopology() const;
  const GeomDetType& type() const override;
  

private:
  FbcmDetId id_;
  FbcmSiPadSpecs* specs_;
};

#endif
