///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "Geometry/FbcmGeometry/interface/FbcmSiPadSpecs.h"

using namespace GeomDetEnumerators;

FbcmSiPadSpecs::FbcmSiPadSpecs(SubDetector rss, const std::string& name, const FbcmSpecs& pars)
    : GeomDetType(name, rss), _p(pars), _n(name) {
  if (rss == FBCM) {
    float pitchX = _p[0];
    float pitchY = _p[1];
    _topology = new FbcmSiPadTopology(pitchX, pitchY);
  } else {
    _topology = nullptr;
  }
}

FbcmSiPadSpecs::~FbcmSiPadSpecs() {
  if (_topology)
    delete _topology;
}

const Topology& FbcmSiPadSpecs::topology() const { return *_topology; }

const FbcmSiPadTopology& FbcmSiPadSpecs::SiPadTopology() const { return *_topology; }

const std::string& FbcmSiPadSpecs::detName() const { return _n; }

const std::vector<float>& FbcmSiPadSpecs::parameters() const { return _p; }
