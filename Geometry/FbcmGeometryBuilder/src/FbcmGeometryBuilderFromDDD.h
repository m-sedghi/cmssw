///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#ifndef Geometry_FbcmGeometry_FbcmGeometryBuilderFromDDD_H
#define Geometry_FbcmGeometry_FbcmGeometryBuilderFromDDD_H

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include <string>
#include <map>
#include <vector>

class DDCompactView;
class DDFilteredView;

class FbcmGeometry;
class FbcmDetId;
class FbcmStationGeom;
class FbcmSiliconDieGeom;
class FbcmSiPadGeom;

//class FbcmGeometryConstants;

class FbcmGeometryBuilderFromDDD {
public:
  FbcmGeometryBuilderFromDDD();

  ~FbcmGeometryBuilderFromDDD();

  void build(FbcmGeometry& theGeometry, const DDCompactView* cview);

private:
  void buildGeometry(FbcmGeometry& theGeometry, DDFilteredView& fview);
  
  std::map<FbcmDetId, std::vector<FbcmDetId>> childs; /// ?

  typedef ReferenceCountingPointer<BoundPlane> FbcmBoundPlane;

  FbcmBoundPlane boundPlane(const DDFilteredView& fv, Bounds* bounds) const;

  FbcmStationGeom* buildStation(DDFilteredView& fv, FbcmDetId detId) const;

  FbcmSiliconDieGeom* buildSiliconDie(DDFilteredView& fv, FbcmDetId detId) const;

  FbcmSiPadGeom* buildSiPad(DDFilteredView& fv, FbcmDetId detId) const;
  
  //FbcmGeometryConstants FbcmGeomConstants;
  
};

#endif
