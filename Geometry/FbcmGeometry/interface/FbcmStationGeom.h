///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#ifndef Geometry_FbcmStationGeom_h
#define Geometry_FbcmStationGeom_h

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h" 

class FbcmSiliconDieGeom;
class FbcmSiPadGeom;

class FbcmStationGeom : public GeomDet {
public:
  /// Constructor
  FbcmStationGeom(FbcmDetId id, const ReferenceCountingPointer<BoundPlane>& plane, unsigned int nDiesPerRing, unsigned int nRings);

  /// Destructor
  ~FbcmStationGeom() override;

  /// Return the FbcmDetId of this Station
  FbcmDetId id() const;

  // Which subdetector
  SubDetector subDetector() const override { return GeomDetEnumerators::FBCM; }

  /// equal if the id is the same
  bool operator==(const FbcmStationGeom& ch) const;

  /// Add SiliconDie to the Station which takes ownership
  void add(FbcmSiliconDieGeom* SiliconDieGeom);

  /// Return the SiPadGeoms in the Station
  std::vector<const GeomDet*> components() const override;

  /// Return the sub-component (SiPadGeoms) with a given id in this Station
  const GeomDet* component(DetId id) const override;

  /// Return the SilcionDie corresponding to the given id
  const FbcmSiliconDieGeom* SiliconDie(FbcmDetId id) const;

  const FbcmSiliconDieGeom* SiliconDie(unsigned int SiDieNo) const;

  /// Return the SiliconDies
  const std::vector<const FbcmSiliconDieGeom*>& SiliconDies() const;

  /// Retun numbers of SiliconDie
  int nSiliconDies() const;

// Add SiPadGeon to the Station taking direct ownership
//  void add(FbcmSiPadGeom* SiPadGeom); // No need!

  /// To support the old
    const FbcmSiPadGeom* SiPad(FbcmDetId id) const;

  const FbcmSiPadGeom* SiPad(unsigned int SiPadNo) const;

  /// To support the old 
  const std::vector<const FbcmSiPadGeom*>& SiPads() const;

    /// Retunr numbers of SiPads
  int nSiPads() const;
  unsigned int NumOfDiesPerRing(void)  const {return nDiesPerRing_;}
  unsigned int NumOfRings(void)  const {return nRings_;}

  //For a line fit in the Station frame, compute: global phi position extrapolated
  //to the last SiliconDie - that extrapolated to the inner SiliconDie
  //float computeDeltaPhi(const LocalPoint& position, const LocalVector& direction) const;
  //?????

private:
  FbcmDetId detId_;
  unsigned int nDiesPerRing_;
  unsigned int nRings_;
  // vector of SiliconDies for a Station
  std::vector<const FbcmSiliconDieGeom*> SiliconDies_;
  // vector of SiPads for a Station
  std::vector<const FbcmSiPadGeom*> SiPads_;
};
#endif
