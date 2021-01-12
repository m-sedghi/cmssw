#ifndef _FBCM_SiPadDigitizerAlgorithm_h
#define _FBCM_SiPadDigitizerAlgorithm_h

///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------


#include <map>
#include <memory>
#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadGeom.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "SimTracker/Common/interface/SiG4UniversalFluctuation.h"

#include "SimFbcm/SiPadDigitizer/plugins/CommonDigiUtility.h"
#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"
#include "SimFbcm/SiPadDigitizer/interface/FftPreparation.h"
#include "SimFbcm/SiPadDigitizer/interface/FbcmFrontEndChip.h"
#include "SimFbcm/SiPadDigitizer/interface/SiHitPulseShape.h"
#include "SimFbcm/SiPadDigitizer/interface/FeConfigSelector.h"
//#include "SimFbcm/SiPadDigitizer/interface/HitAnalysisInfo.h"
#include "DataFormats/FbcmDigi/interface/HitAnalysisInfo.h"

using SiPadDigi = PixelDigi;
using namespace FbcmFE;

// REMEMBER CMS conventions:
// -- Energy: GeV
// -- momentum: GeV/c
// -- mass: GeV/c^2
// -- Distance, position: cm
// -- Time: ns
// -- Angles: radian
// Some constants in convenient units
//constexpr double c_cm_ns = CLHEP::c_light * CLHEP::ns / CLHEP::cm;
//constexpr double c_inv = 1.0 / c_cm_ns;
constexpr double c_inv = 1.0 / 29.9792458 ;

class SiPadDigitizerAlgorithm {
public:
	SiPadDigitizerAlgorithm(const edm::ParameterSet& conf);
  
  ~SiPadDigitizerAlgorithm();

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es);
  void initializeEvent(CLHEP::HepRandomEngine& eng);

  // run the algorithm to digitize a single det
  void accumulateSimHits(const std::vector<PSimHit>::const_iterator inputBegin,
                                 const std::vector<PSimHit>::const_iterator inputEnd,
                                 const size_t inputBeginGlobalIndex,
                                 const unsigned int tofBin,
                                 const FbcmSiPadGeom* SiPadGeom,
                                 const GlobalVector& bfield);
 
  void GetDigiResults(const FbcmSiPadGeom* SiPadGeom, std::map<int, SiPadDigiData>& SiPadDigilMap);	
  // For premixing
  void loadAccumulator(unsigned int detId, const std::map<int, float>& accumulator);

protected:
  //edm::ESHandle<FbcmGeometry> geom_;

  // Accessing Lorentz angle from DB:
  //edm::ESHandle<SiPixelLorentzAngle> SiPixelLorentzAngle_;

  // Internal type aliases
  using signal_map_type = std::map<int, CommonDigiUtility::Amplitude, std::less<int> >;  // from Digi.Skel.
  using signal_map_iterator = signal_map_type::iterator;                                // from Digi.Skel.
  using signal_map_const_iterator = signal_map_type::const_iterator;                    // from Digi.Skel.
  using simlink_map = std::map<unsigned int, std::vector<float>, std::less<unsigned int> >;
  using signalMaps = std::map<uint32_t, signal_map_type>;
  using Frame = GloballyPositioned<double>;
  using Parameters = std::vector<edm::ParameterSet>;

	const edm::ParameterSet& conf_specific;
	const edm::ParameterSet& FFT_SimParam ; // = SiPadDigitizerParam.getParameter<edm::ParameterSet>("FFT_SimParam");
	const edm::ParameterSet& SiHitPulseShapeParam ; // = SiPadDigitizerParam.getParameter<edm::ParameterSet>("SiHitPulseShapeParam");
	const std::vector< edm::ParameterSet > & SiPadFrontEndParamVect; //  = SiPadDigitizerParam.getParameter< std::vector< edm::ParameterSet > >("SiPadFrontEndParam");

	FftPreparation FftPrep; //(FFT_SimParam); 
	SiHitPulseShape HitPulse; // ( SiHitPulseShapeParam.getParameter< std::vector<double> >("HitPulseParam") );
	FbcmFrontEndChip FrontEnd; //(FftPrep); // the FrontEnd parameters will be set just before running for each sensor size
	FeConfigSelector FeParamSelector; //(SiPadFrontEndParamVect);
	int FirstBxSlotNo;
	int LastBxSlotNo;
  
  // Contains the accumulated hit info.
  signalMaps _signal;

  //const bool makeDigiSimLinks_;
  //-- charge fluctuation
  const double tMax;  // The delta production cut, should be as in OSCAR = 30keV
  const bool use_LorentzAngle_DB_;  // if we want to get Lorentz angle from the DataBase.
  const float GeVperElectron;  // 3.7E-09 // go from Geant energy GeV to number of electrons
  //-- drift
  const bool alpha2Order;  // Switch on/off of E.B effect
  const float Sigma0;      //=0.0007  // Charge diffusion in microns for 300 micron Si
  const float SigmaCoeff;  // delta in the diffusion across the strip pitch
  //-- induce_signal
  const float ClusterWidth;  // Gaussian charge cutoff width in sigma units

  const float theTailNoiseInElectrons; 
  const float theReadoutNoise;      // Noise of the readount chain in elec,

  const int hitSelectionMode_;
  const float theTofLowerCut;                  // Cut on the particle TOF
  const float theTofUpperCut;                  // Cut on the particle TOF
  const float tanLorentzAnglePerTesla_;  // Lorentz angle tangent per Tesla at Fbcm location
  // -- add_noise
  const bool addNoise;
  const bool fluctuateCharge;
  const std::unique_ptr<SiG4UniversalFluctuation> fluctuate;  // make a pointer
  //const std::unique_ptr<GaussianTailNoiseGenerator> theNoiser;
  // pseudoRadDamage
  const double pseudoRadDamage;        // Decrease the amount off freed charge that reaches the collector
  const double pseudoRadDamageRadius;  // Only apply pseudoRadDamage to pixels with radius<=pseudoRadDamageRadius

  bool FilterHit(const PSimHit& hit, double tCorr); 
  
  //-- additional member functions
  void primary_ionization(const PSimHit& hit,
						  std::vector<CommonDigiUtility::EnergyDepositUnit>& ionization_points) const;
  void drift(const PSimHit& hit, 
             const FbcmSiPadGeom* SiPadGeom,
             const GlobalVector& bfield,
             const std::vector<CommonDigiUtility::EnergyDepositUnit>& ionization_points,
             std::vector<CommonDigiUtility::SignalPoint>& collection_points) const;
  void induce_signal(const PSimHit& hit,
                     const size_t hitIndex,
                     const unsigned int tofBin,
                     const FbcmSiPadGeom* SiPadGeom,
                     const std::vector<CommonDigiUtility::SignalPoint>& collection_points);
  void fluctuateEloss(int particleId,
                      float momentum,
                      float eloss,
                      float length,
                      int NumberOfSegments,
                      std::vector<float>& elossVector) const;
   void add_noise(const FbcmSiPadGeom* SiPadGeom);
   
   LocalVector DriftDirection(const FbcmSiPadGeom* SiPadGeom,
                             const GlobalVector& bfield,
                             const DetId& detId) const;

  // For random numbers
  std::unique_ptr<CLHEP::RandGaussQ> gaussDistribution_;
  //for engine passed into the constructor from Digitizer
  CLHEP::HepRandomEngine* rengine_;

  double calcQ(float x) const {
    auto xx = std::min(0.5f * x * x, 12.5f);
    return 0.5 * (1.0 - std::copysign(std::sqrt(1.f - unsafe_expf<4>(-xx * (1.f + 0.2733f / (1.f + 0.147f * xx)))), x));
  }
  
};
#endif
