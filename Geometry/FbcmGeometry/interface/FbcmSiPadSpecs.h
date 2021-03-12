///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#ifndef Geometry_FbcmGeometry_FbcmSiPadSpecs_H
#define Geometry_FbcmGeometry_FbcmSiPadSpecs_H

#include <vector>
#include <string>
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadTopology.h"

class FbcmSiPadSpecs : public GeomDetType {
public:
  typedef std::vector<float> FbcmSpecs;

  FbcmSiPadSpecs(SubDetector rss, const std::string& name, const FbcmSpecs& pars);

  ~FbcmSiPadSpecs() override;

  const Topology& topology() const override;
  
  const FbcmSiPadTopology& SiPadTopology() const ;
  
  const std::string& detName() const;

  const FbcmSpecs& parameters() const;

private:
  FbcmSiPadTopology* _topology;

  std::vector<float> _p;
  std::string _n; ///
};
#endif
