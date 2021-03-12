
#ifndef Geometry_FbcmGeometry_h
#define Geometry_FbcmGeometry_h

///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------
/**
FbcmGeometry Class: this class stores the geometries of FBCM
*/


#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h" 
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadGeom.h" 
#include "Geometry/FbcmGeometry/interface/FbcmStationGeom.h" 
#include "Geometry/FbcmGeometry/interface/FbcmSiliconDieGeom.h" 
#include <vector>
#include <map>

/**
Note: the followig type was defined in  TrackingGeometry.h
  using DetTypeContainer = std::vector<const GeomDetType*>;
  using DetContainer = std::vector<const GeomDet*>;
  using DetIdContainer = std::vector<DetId>;
  using mapIdToDetUnit = std::unordered_map<unsigned int, const GeomDet*>;
  using mapIdToDet = std::unordered_map<unsigned int, const GeomDet*>;
  */
  
  
  using mapIdToStation = std::unordered_map<unsigned int, const FbcmStationGeom*>;
  using mapIdToSiDie = std::unordered_map<unsigned int, const FbcmSiliconDieGeom*>;
  using mapIdToSiPad = std::unordered_map<unsigned int, const FbcmSiPadGeom*>;
  

class FbcmGeometry : public TrackingGeometry {
public:
  /// Default constructor
  FbcmGeometry();

  /// Destructor
  ~FbcmGeometry() override;

  // Return a vector of all det types
  const DetTypeContainer& detTypes() const override;

  // Return a vector of all GeomDetUnit
  const DetContainer& detUnits() const override;

  // Return a vector of all GeomDet
  const DetContainer& dets() const override;

  // Return a vector of all GeomDetUnit DetIds
  const DetIdContainer& detUnitIds() const override;

  // Return a vector of all GeomDet DetIds
  const DetIdContainer& detIds() const override;

  // Return the pointer to the GeomDetUnit corresponding to a given DetId
  const GeomDet* idToDetUnit(DetId) const override;

  // Return the pointer to the GeomDet corresponding to a given DetId
  const GeomDet* idToDet(DetId) const override;

  //---- Extension of the interface

  /// Return a FbcmStationGeom given its FbcmDetId
  const FbcmStationGeom* IdToStation(FbcmDetId id) const;

  /// Return a FbcmSiliconDieGeom given its FbcmDetId
  const FbcmSiliconDieGeom* IdToSiliconDie(FbcmDetId id) const;

  /// Return a FbcmSiPadGeom given its FbcmDetId
  const FbcmSiPadGeom* IdToSiPad(FbcmDetId id) const;

  /// Return a vector of all FbcmSiPad Geometries
  const std::vector<FbcmSiPadGeom const*>& SiPads() const;

  /// Return a vector of all FbcmStation Geometries
  const std::vector<const FbcmStationGeom*>& Stations() const;

  /// Return a vector of all FbcmSiliconDie Geometries
  const std::vector<const FbcmSiliconDieGeom*>& SiliconDies() const;

  /// Add a FbcmSiPadGeom to the Geometry
  void add(FbcmSiPadGeom* SiPadGeom);

  /// Add a FbcmStationGeom to the Geometry
  void add(FbcmStationGeom* StationGeom);

  /// Add a FbcmSiliconDieGeom to the Geometry
  void add(FbcmSiliconDieGeom* SiliconDieGeom);
  
  void SetNumOfStations(int nStations) {nStations_=nStations;}
  
  int NumOfStations(void) const {return nStations_;}

private:

  DetContainer theSiPadGeoms;
  DetTypeContainer theSiPadTypes;
  DetIdContainer theSiPadIds;
  DetIdContainer theDetIds;
  DetContainer theDets;

  
	mapIdToDet theMap;
  
  mapIdToStation theMapToStations;
  mapIdToSiDie theMapToSiDies;
  mapIdToSiPad theMapToSiPads;

  std::vector<FbcmStationGeom const*> allStationGeoms; 
  std::vector<FbcmSiliconDieGeom const*> allSiliconDieGeoms;  
  std::vector<FbcmSiPadGeom const*> allSiPadGeoms;  
  int nStations_;
  
};

#endif
