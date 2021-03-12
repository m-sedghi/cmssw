///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#ifndef Geometry_FbcmSiliconDieGeom_h
#define Geometry_FbcmSiliconDieGeom_h

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h" 

class FbcmSiPadGeom; // ??

class FbcmSiliconDieGeom : public GeomDet {
public:
  /// Constructor
  FbcmSiliconDieGeom(FbcmDetId id, const ReferenceCountingPointer<BoundPlane>& plane, unsigned int nCols, unsigned int nRows);

  /// Destructor
  ~FbcmSiliconDieGeom() override;

  /// Return the FbcmDetId of this SiliconDie
  FbcmDetId id() const;

  // Which subdetector
  SubDetector subDetector() const override { return GeomDetEnumerators::FBCM; }

  /// equal if the id is the same
  bool operator==(const FbcmSiliconDieGeom& ch) const;

  /// Add SiPadGeom to the SiliconDie which takes ownership
  void add(const FbcmSiPadGeom* SiPadGeom);

  /// Return the SiPads in the SiliconDie
  std::vector<const GeomDet*> components() const override;

  /// Return the sub-component (SiPad) with a given id in this SiliconDie
  const GeomDet* component(DetId id) const override;

  /// Return the SiPadGeom corresponding to the given id
  const FbcmSiPadGeom* SiPad(FbcmDetId id) const;

  const FbcmSiPadGeom* SiPad(unsigned int SiPadNo) const;

  /// Return the SiPads
  const std::vector<const FbcmSiPadGeom*>& SiPads() const;

  /// Retunr numbers of SiPads
  int nSiPads() const;
   unsigned int NumOfCols() const {return nCols_; }
   unsigned int NumOfRows()  const {return nRows_;}
   
private:
  FbcmDetId detId_;
   unsigned int nCols_;
   unsigned int nRows_;

  // vector of SiPadGeoms for a SiliconDie
  std::vector<const FbcmSiPadGeom*> SiPads_;
};
#endif
