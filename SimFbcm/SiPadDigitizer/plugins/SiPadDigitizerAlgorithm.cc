///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <typeinfo>
#include <iostream>
#include <cmath>


#include "SiPadDigitizerAlgorithm.h"

using namespace edm;

SiPadDigitizerAlgorithm::SiPadDigitizerAlgorithm(const edm::ParameterSet& conf) :
	conf_specific(conf.getParameter<ParameterSet>("SiPadSimParam")),
	FFT_SimParam(conf.getParameter<edm::ParameterSet>("FFT_SimParam") ),
	SiHitPulseShapeParam(conf.getParameter<edm::ParameterSet>("SiHitPulseShapeParam")),
	SiPadFrontEndParamVect(conf.getParameter< std::vector< edm::ParameterSet > >("SiPadFrontEndParam")),
	FftPrep(FFT_SimParam),
	HitPulse( SiHitPulseShapeParam.getParameter< std::vector<double> >("HitPulseParam") ),
    FrontEnd(FftPrep), // the FrontEnd parameters will be set just before running each sensor size
	FeParamSelector(SiPadFrontEndParamVect),
	FirstBxSlotNo(conf_specific.getParameter< int >("FirstBxSlotNo")),
	LastBxSlotNo(conf_specific.getParameter< int >("LastBxSlotNo")),
	_signal(),
    tMax(conf_specific.getParameter<double> ("DeltaProductionCut")), // delta cutoff in MeV, has to be same as in OSCAR(0.030/cmsim=1.0 MeV  // tMax(0.030) In MeV.
    use_LorentzAngle_DB_(false),
    GeVperElectron(3.61E-09),                                      // 1 electron(3.61eV, 1keV(277e, mod 9/06 d.k.
    alpha2Order(conf_specific.getParameter<bool>("Alpha2Order")),  // switch on/off of E.B effect
    Sigma0(conf_specific.getParameter<double>("SigmaZero")),       // Charge diffusion constant 7->3.7
    SigmaCoeff(conf_specific.getParameter<double>("SigmaCoeff")),  // delta in the diffusion across the strip pitch
    // (set between 0 to 0.9,  0-->flat Sigma0, 1-->Sigma_min=0 & Sigma_max=2*Sigma0
    ClusterWidth(conf_specific.getParameter<double>("ClusterWidth")), // Charge integration spread on the collection plane
   // Noise in electrons:
   //theNoiseInElectrons(conf_specific.getParameter<double>("NoiseInElectrons")), // Pixel cell noise, relevant for generating noisy pixels
	theTailNoiseInElectrons(conf_specific.getParameter<double>("GaussianTailNoise")),
    theReadoutNoise(conf_specific.getParameter<double>("ReadoutNoiseInElec")), // Fill readout noise, including all readout chain
    hitSelectionMode_(conf_specific.getParameter<int>("HitSelectionMode")),
    // theTofCut 12.5, cut in particle TOD +/- 12.5ns
    theTofLowerCut(conf_specific.getParameter<double>("TofLowerCut")),
    theTofUpperCut(conf_specific.getParameter<double>("TofUpperCut")),
     // Get the Lorentz angle from the cfg file:
    tanLorentzAnglePerTesla_( use_LorentzAngle_DB_ ? 0.0 : conf_specific.getParameter<double>("TanLorentzAnglePerTesla_Fbcm")),
    addNoise(conf_specific.getParameter<bool>("AddNoise")), // Add noise
    // Fluctuate charge in track subsegments
    fluctuateCharge(conf_specific.getUntrackedParameter<bool>("FluctuateCharge", true)),
    fluctuate(fluctuateCharge ? new SiG4UniversalFluctuation() : nullptr),
	// Add some pseudo-red damage	
    pseudoRadDamage(conf_specific.getUntrackedParameter<double>("PseudoRadDamage", 0.0)),
	pseudoRadDamageRadius(conf_specific.getUntrackedParameter<double>("PseudoRadDamageRadius", 0.0)),
	chargeCollectionEff(conf_specific.getParameter<double>("chargeCollectionEfficiency"))
	{

	LogInfo("SiPadDigitizerAlgorithm") << "SiPadDigitizerAlgorithm constructed\n" ;
}

void SiPadDigitizerAlgorithm::init(const edm::EventSetup& es) 
   { 
   // if the FbcmGeometry is needed, then it can be retrieved here.
   //es.get<FbcmGeometryRecord>().get(geom_); 
   }


SiPadDigitizerAlgorithm::~SiPadDigitizerAlgorithm() {
  LogDebug("SiPadDigitizerAlgorithm") << "SiPadDigitizerAlgorithm deleted";
}

void SiPadDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
                                              std::vector<PSimHit>::const_iterator inputEnd,
                                              const size_t inputBeginGlobalIndex,
                                              const unsigned int tofBin,
                                              const FbcmSiPadGeom* SiPadGeom,
                                              const GlobalVector& bfield) 
{
	
  // produce SignalPoint's for all SimHit's in detector
  // Loop over hits
  uint32_t detId = SiPadGeom->geographicalId().rawId();
  size_t simHitGlobalIndex = inputBeginGlobalIndex;  // This needs to be stored to create the digi-sim link later. No need any more!
  for (auto it = inputBegin; it != inputEnd; ++it, ++simHitGlobalIndex) {
    // skip hits not in this detector.
    if ((*it).detUnitId() != detId)
      continue;

    LogDebug("PSPDigitizerAlgorithm") << (*it).particleType() << " " << (*it).pabs() << " " << (*it).energyLoss() << " "
                                      << (*it).tof() << " " << (*it).trackId() << " " << (*it).processType() << " "
                                      << (*it).detUnitId() << (*it).entryPoint() << " " << (*it).exitPoint();
									  
	std::vector<CommonDigiUtility::EnergyDepositUnit> ionization_points;
    std::vector<CommonDigiUtility::SignalPoint> collection_points;

    // fill collection_points for this SimHit, indpendent of topology
    double tCorr = SiPadGeom->surface().toGlobal((*it).localPosition()).mag() * c_inv;
	
	if (FilterHit((*it), tCorr)) {
	  primary_ionization(*it, ionization_points);  // fills _ionization_points
      // transforms _ionization_points to collection_points
      drift(*it, SiPadGeom, bfield, ionization_points, collection_points);
      // compute induced signal on readout elements and add to _signal
      induce_signal(*it, simHitGlobalIndex, tofBin, SiPadGeom, collection_points);
    }
  }
}

// =================================================================
// Generate primary ionization along the track segment.
// Divide the track into small sub-segments
void SiPadDigitizerAlgorithm::primary_ionization(
    const PSimHit& hit, std::vector<CommonDigiUtility::EnergyDepositUnit>& ionization_points) const {
  // Straight line approximation for trajectory inside active media
  const float SegmentLength = 0.0010;  // in cm (10 microns)
  float energy;

  // Get the 3D segment direction vector
  LocalVector direction = hit.exitPoint() - hit.entryPoint();

  float eLoss = hit.energyLoss();  // Eloss in GeV
  float length = direction.mag();  // Track length in Silicon

  int NumberOfSegments = int(length / SegmentLength);  // Number of segments
  if (NumberOfSegments < 1)
    NumberOfSegments = 1;
  LogDebug("SiPadDigitizerAlgorithm")
      << "enter primary_ionzation " << NumberOfSegments << " shift = " << hit.exitPoint().x() - hit.entryPoint().x()
      << " " << hit.exitPoint().y() - hit.entryPoint().y() << " " << hit.exitPoint().z() - hit.entryPoint().z() << " "
      << hit.particleType() << " " << hit.pabs();

  std::vector<float> elossVector;  // Eloss vector
  elossVector.reserve(NumberOfSegments);
  if (fluctuateCharge) {
    int pid = hit.particleType();
    // int pid=211;  // assume it is a pion

    float momentum = hit.pabs();
    // Generate fluctuated charge points
    fluctuateEloss(pid, momentum, eLoss, length, NumberOfSegments, elossVector);
  }
  ionization_points.reserve(NumberOfSegments);  // set size

  // loop over segments
  for (int i = 0; i != NumberOfSegments; ++i) {
    // Divide the segment into equal length subsegments
    Local3DPoint point = hit.entryPoint() + float((i + 0.5) / NumberOfSegments) * direction;
    if (fluctuateCharge)
      energy = elossVector[i] / GeVperElectron;  // Convert charge to elec.
    else
      energy = hit.energyLoss() / GeVperElectron / float(NumberOfSegments);

    CommonDigiUtility::EnergyDepositUnit edu(energy, point);  // define position,energy point
    ionization_points.push_back(edu);                        // save
    LogDebug("SiPadDigitizerAlgorithm")
        << i << " " << ionization_points[i].x() << " " << ionization_points[i].y() << " " << ionization_points[i].z()
        << " " << ionization_points[i].energy();
  }
}
//==============================================================================
// Fluctuate the charge comming from a small (10um) track segment.
// Use the G4 routine. For mip pions for the moment.
void SiPadDigitizerAlgorithm::fluctuateEloss(int pid,
                                                     float particleMomentum,
                                                     float eloss,
                                                     float length,
                                                     int NumberOfSegs,
                                                     std::vector<float>& elossVector) const {
  // Get dedx for this track
  //float dedx;
  //if( length > 0.) dedx = eloss/length;
  //else dedx = eloss;

  double particleMass = 139.6;  // Mass in MeV, Assume pion
  pid = std::abs(pid);
  if (pid != 211) {  // Mass in MeV
    if (pid == 11)
      particleMass = 0.511;
    else if (pid == 13)
      particleMass = 105.7;
    else if (pid == 321)
      particleMass = 493.7;
    else if (pid == 2212)
      particleMass = 938.3;
  }
  // What is the track segment length.
  float segmentLength = length / NumberOfSegs;

  // Generate charge fluctuations.
  float de = 0.;
  float sum = 0.;
  double segmentEloss = (1000. * eloss) / NumberOfSegs;  //eloss in MeV
  for (int i = 0; i < NumberOfSegs; ++i) {
    //       material,*,   momentum,energy,*, *,  mass
    //myglandz_(14.,segmentLength,2.,2.,dedx,de,0.14);
    // The G4 routine needs momentum in MeV, mass in Mev, delta-cut in MeV,
    // track segment length in mm, segment eloss in MeV
    // Returns fluctuated eloss in MeV
    double deltaCutoff = tMax;  // the cutoff is sometimes redefined inside, so fix it.
    de = fluctuate->SampleFluctuations(static_cast<double>(particleMomentum * 1000.),
                                       particleMass,
                                       deltaCutoff,
                                       static_cast<double>(segmentLength * 10.),
                                       segmentEloss,
                                       rengine_) /
         1000.;  //convert to GeV
    elossVector.push_back(de);
    sum += de;
  }
  if (sum > 0.) {  // If fluctuations give eloss>0.
    // Rescale to the same total eloss
    float ratio = eloss / sum;
    for (int ii = 0; ii < NumberOfSegs; ++ii)
      elossVector[ii] = ratio * elossVector[ii];
  } else {  // If fluctuations gives 0 eloss
    float averageEloss = eloss / NumberOfSegs;
    for (int ii = 0; ii < NumberOfSegs; ++ii)
      elossVector[ii] = averageEloss;
  }
}

// ======================================================================
// Drift the charge segments to the sensor surface (collection plane)
// Include the effect of E-field and B-field
void SiPadDigitizerAlgorithm::drift(const PSimHit& hit,
                                            const FbcmSiPadGeom* SiPadGeom,
                                            const GlobalVector& bfield,
                                            const std::vector<CommonDigiUtility::EnergyDepositUnit>& ionization_points,
                                            std::vector<CommonDigiUtility::SignalPoint>& collection_points) const {
  LogDebug("SiPadDigitizerAlgorithm") << "enter drift ";

  collection_points.resize(ionization_points.size());                      // set size
  LocalVector driftDir = DriftDirection(SiPadGeom, bfield, hit.detUnitId());  // get the charge drift direction
  if (driftDir.z() == 0.) {
    LogWarning("SiPadDigitizerAlgorithm") << " pxlx: drift in z is zero ";
    return;
  }

  float TanLorenzAngleX, TanLorenzAngleY, dir_z, CosLorenzAngleX, CosLorenzAngleY;
  if (alpha2Order) {
    TanLorenzAngleX = driftDir.x();  // tangen of Lorentz angle
    TanLorenzAngleY = driftDir.y();
    dir_z = driftDir.z();                                                      // The z drift direction
    CosLorenzAngleX = 1. / std::sqrt(1. + TanLorenzAngleX * TanLorenzAngleX);  // cosine
    CosLorenzAngleY = 1. / std::sqrt(1. + TanLorenzAngleY * TanLorenzAngleY);  // cosine;
  } else {
    TanLorenzAngleX = driftDir.x();
    TanLorenzAngleY = 0.;                                                      // force to 0, driftDir.y()/driftDir.z();
    dir_z = driftDir.z();                                                      // The z drift direction
    CosLorenzAngleX = 1. / std::sqrt(1. + TanLorenzAngleX * TanLorenzAngleX);  // cosine to estimate the path length
    CosLorenzAngleY = 1.;
  }

  float moduleThickness = SiPadGeom->specificSurface().bounds().thickness();
  float stripPitch = SiPadGeom->SiPadTopology().pitch().first;

  LogDebug("SiPadDigitizerAlgorithm")
      << " Lorentz Tan " << TanLorenzAngleX << " " << TanLorenzAngleY << " " << CosLorenzAngleX << " "
      << CosLorenzAngleY << " " << moduleThickness * TanLorenzAngleX << " " << driftDir;

  float Sigma_x = 1.;  // Charge spread
  float Sigma_y = 1.;
  float DriftDistance;  // Distance between charge generation and collection
  float DriftLength;    // Actual Drift Lentgh
  float Sigma;

  for (unsigned int i = 0; i != ionization_points.size(); ++i) {
    float SegX, SegY, SegZ;  // position
    SegX = ionization_points[i].x();
    SegY = ionization_points[i].y();
    SegZ = ionization_points[i].z();

    // Distance from the collection plane
    // DriftDistance = (moduleThickness/2. + SegZ); // Drift to -z
    // Include explixitely the E drift direction (for CMS dir_z=-1)
    DriftDistance = moduleThickness / 2. - (dir_z * SegZ);  // Drift to -z

    if (DriftDistance < 0.)
      DriftDistance = 0.;
    else if (DriftDistance > moduleThickness)
      DriftDistance = moduleThickness;

    // Assume full depletion now, partial depletion will come later.
    float XDriftDueToMagField = DriftDistance * TanLorenzAngleX;
    float YDriftDueToMagField = DriftDistance * TanLorenzAngleY;

    // Shift cloud center
    float CloudCenterX = SegX + XDriftDueToMagField;
    float CloudCenterY = SegY + YDriftDueToMagField;

    // Calculate how long is the charge drift path
    DriftLength = std::sqrt(DriftDistance * DriftDistance + XDriftDueToMagField * XDriftDueToMagField +
                            YDriftDueToMagField * YDriftDueToMagField);

    // What is the charge diffusion after this path
    // Sigma0=0.00037 is for 300um thickness (make sure moduleThickness is in [cm])
    Sigma = std::sqrt(DriftLength / moduleThickness) * (Sigma0 * moduleThickness / 0.0300);
    // D.B.: sigmaCoeff=0 means no modulation
    if (SigmaCoeff)
      Sigma *= (SigmaCoeff * cos(SegX * M_PI / stripPitch) * cos(SegX * M_PI / stripPitch) + 1);
    // NB: divided by 4 to get a periodicity of stripPitch

    // Project the diffusion sigma on the collection plane
    Sigma_x = Sigma / CosLorenzAngleX;
    Sigma_y = Sigma / CosLorenzAngleY;

    // Insert a charge loss due to Rad Damage here
    float energyOnCollector = ionization_points[i].energy();  // The energy that reaches the collector

    // pseudoRadDamage
    if (pseudoRadDamage >= 0.001) {
      float moduleRadius = SiPadGeom->surface().position().perp();
      if (moduleRadius <= pseudoRadDamageRadius) {
        float kValue = pseudoRadDamage / (moduleRadius * moduleRadius);
        energyOnCollector = energyOnCollector * exp(-1 * kValue * DriftDistance / moduleThickness);
      }
    }
    LogDebug("SiPadDigitizerAlgorithm")
        << "Dift DistanceZ = " << DriftDistance << " module thickness = " << moduleThickness
        << " Start Energy = " << ionization_points[i].energy() << " Energy after loss= " << energyOnCollector;
    CommonDigiUtility::SignalPoint sp(CloudCenterX, CloudCenterY, Sigma_x, Sigma_y, hit.tof(), energyOnCollector);
    // Load the Charge distribution parameters
    collection_points[i] = sp;
  }
}

// ====================================================================
// Induce the signal on the collection plane of the active sensor area.
void SiPadDigitizerAlgorithm::induce_signal(
    const PSimHit& hit,
    const size_t hitIndex,
    const unsigned int tofBin,
    const FbcmSiPadGeom* SiPadGeom,
    const std::vector<CommonDigiUtility::SignalPoint>& collection_points) {
    const FbcmSiPadTopology* topol = &SiPadGeom->SiPadTopology();
    uint32_t detID = SiPadGeom->geographicalId().rawId();
    signal_map_type& theSignal = _signal[detID];

  LogDebug("SiPadDigitizerAlgorithm")
      << " enter induce_signal, " << topol->pitch().first << " " << topol->pitch().second;  //OK
	
  // local map to store pixels hit by 1 Hit.
  using hit_map_type = std::map<int, float, std::less<int> >;
  hit_map_type hit_signal;

  // map to store pixel integrals in the x and in the y directions
  std::map<int, float, std::less<int> > x, y;

  // Assign signals to readout channels and store sorted by channel number
  int iseg = 0;
  float ESum = 0.0;

  // Iterate over collection points on the collection plane
  for (auto const& v : collection_points) {
    iseg++;
    float CloudCenterX = v.position().x();  // Charge position in x
    float CloudCenterY = v.position().y();  //                 in y
    float SigmaX = v.sigma_x();             // Charge spread in x
    float SigmaY = v.sigma_y();             //               in y
    float Charge = v.amplitude();           // Charge amplitude

    LogDebug("SiPadDigitizerAlgorithm") << " cloud " << v.position().x() << " " << v.position().y() << " "
                                                << v.sigma_x() << " " << v.sigma_y() << " " << v.amplitude();

    // Find the maximum cloud spread in 2D plane , assume 3*sigma
    float CloudRight = CloudCenterX + ClusterWidth * SigmaX;
    float CloudLeft = CloudCenterX - ClusterWidth * SigmaX;
    float CloudUp = CloudCenterY + ClusterWidth * SigmaY;
    float CloudDown = CloudCenterY - ClusterWidth * SigmaY;

    // Define 2D cloud limit points
    LocalPoint PointRightUp = LocalPoint(CloudRight, CloudUp);
    LocalPoint PointLeftDown = LocalPoint(CloudLeft, CloudDown);

    // This points can be located outside the sensor area.
    // The conversion to measurement point does not check for that
    // so the returned pixel index might be wrong (outside range).
    // We rely on the limits check below to fix this.
    // But remember whatever we do here THE CHARGE OUTSIDE THE ACTIVE
    // PIXEL ARE IS LOST, it should not be collected.


    // Convert the 2D points to pixel indices
    MeasurementPoint mp = topol->measurementPosition(PointRightUp);  //OK
	
    int IPixRightUpX = int(floor(mp.x()));
    int IPixRightUpY = int(floor(mp.y()));

    LogDebug("SiPadDigitizerAlgorithm")
        << " right-up " << PointRightUp << " " << mp.x() << " " << mp.y() << " " << IPixRightUpX << " " << IPixRightUpY;

    mp = topol->measurementPosition(PointLeftDown);  // OK

    int IPixLeftDownX = int(floor(mp.x()));
    int IPixLeftDownY = int(floor(mp.y()));

    LogDebug("SiPadDigitizerAlgorithm") << " left-down " << PointLeftDown << " " << mp.x() << " " << mp.y()
                                                << " " << IPixLeftDownX << " " << IPixLeftDownY;

    // Check detector limits to correct for pixels outside range.
    int numColumns = topol->ncolumns();  // det module number of cols&rows
    int numRows = topol->nrows();

    IPixRightUpX = numRows > IPixRightUpX ? IPixRightUpX : numRows - 1;
    IPixRightUpY = numColumns > IPixRightUpY ? IPixRightUpY : numColumns - 1;
    IPixLeftDownX = 0 < IPixLeftDownX ? IPixLeftDownX : 0;
    IPixLeftDownY = 0 < IPixLeftDownY ? IPixLeftDownY : 0;

    x.clear();  // clear temporary integration array
    y.clear();

    // First integrate cahrge strips in x
    int ix;                                               // TT for compatibility
    for (ix = IPixLeftDownX; ix <= IPixRightUpX; ix++) {  // loop over x index
      float xUB, xLB, UpperBound, LowerBound;

      if (ix == 0 || SigmaX == 0.)  // skip for surface segemnts
        LowerBound = 0.;
      else {
        mp = MeasurementPoint(float(ix), 0.0);
        xLB = topol->localPosition(mp).x();
        LowerBound = 1 - calcQ((xLB - CloudCenterX) / SigmaX);
      }

      if (ix == numRows - 1 || SigmaX == 0.)
        UpperBound = 1.;
      else {
        mp = MeasurementPoint(float(ix + 1), 0.0);
        xUB = topol->localPosition(mp).x();
        UpperBound = 1. - calcQ((xUB - CloudCenterX) / SigmaX);
      }
      float TotalIntegrationRange = UpperBound - LowerBound;  // get strip
      x[ix] = TotalIntegrationRange;                          // save strip integral
    }

    // Now integarte strips in y
    int iy;                                               // TT for compatibility
    for (iy = IPixLeftDownY; iy <= IPixRightUpY; iy++) {  //loope over y ind
      float yUB, yLB, UpperBound, LowerBound;

      if (iy == 0 || SigmaY == 0.)
        LowerBound = 0.;
      else {
        mp = MeasurementPoint(0.0, float(iy));
        yLB = topol->localPosition(mp).y();
        LowerBound = 1. - calcQ((yLB - CloudCenterY) / SigmaY);
      }

      if (iy == numColumns - 1 || SigmaY == 0.)
        UpperBound = 1.;
      else {
        mp = MeasurementPoint(0.0, float(iy + 1));
        yUB = topol->localPosition(mp).y();
        UpperBound = 1. - calcQ((yUB - CloudCenterY) / SigmaY);
      }

      float TotalIntegrationRange = UpperBound - LowerBound;
      y[iy] = TotalIntegrationRange;  // save strip integral
    }

    // Get the 2D charge integrals by folding x and y strips
    int chan;
    for (ix = IPixLeftDownX; ix <= IPixRightUpX; ix++) {    // loop over x index
      for (iy = IPixLeftDownY; iy <= IPixRightUpY; iy++) {  //loope over y ind
        float ChargeFraction = Charge * x[ix] * y[iy];
		
        if (ChargeFraction > 0.) {
		 chan=SiPadDigi::pixelToChannel(ix, iy);
          // Load the amplitude
          hit_signal[chan] += ChargeFraction;
        }

        mp = MeasurementPoint(float(ix), float(iy));
        LocalPoint lp = topol->localPosition(mp);
        chan = topol->channel(lp);

        //LogDebug("SiPadDigitizerAlgorithm")
		    // << " Ix Iy" << ix << " " << iy << " -chan: "
            // << chan << ", ChargeFraction:" << ChargeFraction << "mp.x,y: " << mp.x() << "," << mp.y() << "lp.x,y: " << lp.x() << "," << lp.y()
            // << "ch: "    // givex edge position
            // << chan << "\n";  // edge belongs to previous ?
        ESum += ChargeFraction;
      }
    }
  }
  
  // Fill the global map with all hit pixels from this event
  float corr_time = hit.tof() - SiPadGeom->surface().toGlobal(hit.localPosition()).mag() * c_inv;
  for (auto const& hit_s : hit_signal) {
    int chan = hit_s.first;
	theSignal[chan] += CommonDigiUtility::Amplitude(hit_s.second, &hit, hit_s.second, corr_time, hitIndex, tofBin); 
  }
}

// ======================================================================
//  Add electronic noise to pixel charge
void SiPadDigitizerAlgorithm::add_noise(const FbcmSiPadGeom* SiPadGeom) {
  uint32_t detID = SiPadGeom->geographicalId().rawId();
  signal_map_type& theSignal = _signal[detID];
  for (auto& s : theSignal) {
    float noise = gaussDistribution_->fire();
    if ((s.second.ampl() + noise) < 0.)
      s.second.set(0);
    else
      s.second += noise;
  }
}

void SiPadDigitizerAlgorithm::initializeEvent(CLHEP::HepRandomEngine& eng) {
  if (addNoise || fluctuateCharge) 
  {
    gaussDistribution_ = std::make_unique<CLHEP::RandGaussQ>(eng, 0., theReadoutNoise);
  } 
  rengine_ = (&eng);
  _signal.clear();
}

// =======================================================================================
// Set the drift direction accoring to the Bfield in local det-unit frame
// Works for both barrel and forward pixels.
// Replace the sign convention to fit M.Swartz's formulaes.
LocalVector SiPadDigitizerAlgorithm::DriftDirection(const FbcmSiPadGeom* SiPadGeom,
                                                            const GlobalVector& bfield,
                                                            const DetId& detId) const {
  Frame detFrame(SiPadGeom->surface().position(), SiPadGeom->surface().rotation());
  LocalVector Bfield = detFrame.toLocal(bfield);
  float alpha2_Endcap;
  //float alpha2;

  float dir_x = 0.0;
  float dir_y = 0.0;
  float dir_z = 0.0;
  float scale = 0.0;

  uint32_t detID = SiPadGeom->geographicalId().rawId();
  unsigned int Sub_detid = DetId(detID).subdetId();

  // Read Lorentz angle from cfg file:
  if (!use_LorentzAngle_DB_) {
    if (alpha2Order) {
      alpha2_Endcap = tanLorentzAnglePerTesla_ * tanLorentzAnglePerTesla_;
    } else {
      alpha2_Endcap = 0.0;
    }
	if (Sub_detid == FbcmSubdetId::FbcmModule) { // similar to Forward disks in the Tracker-Endcap.
		  dir_x = -(tanLorentzAnglePerTesla_ * Bfield.y() + alpha2_Endcap * Bfield.z() * Bfield.x());
      dir_y = +(tanLorentzAnglePerTesla_ * Bfield.x() - alpha2_Endcap * Bfield.z() * Bfield.y());
      dir_z = -(1 + alpha2_Endcap * Bfield.z() * Bfield.z());
      scale = (1 + alpha2_Endcap * Bfield.z() * Bfield.z());
	}
	else {
		throw cms::Exception("SiPadDigitizer Error") << "Wrong Sub_detid\n"; 
		}
  }

/*
  // Read Lorentz angle from DB:
  if (use_LorentzAngle_DB_) {
    float lorentzAngle = SiPixelLorentzAngle_->getLorentzAngle(detId);
    alpha2 = lorentzAngle * lorentzAngle;

    dir_x = -(lorentzAngle * Bfield.y() + alpha2 * Bfield.z() * Bfield.x());
    dir_y = +(lorentzAngle * Bfield.x() - alpha2 * Bfield.z() * Bfield.y());
    dir_z = -(1 + alpha2 * Bfield.z() * Bfield.z());
    scale = (1 + alpha2 * Bfield.z() * Bfield.z());
  }
*/
  LocalVector theDriftDirection = LocalVector(dir_x / scale, dir_y / scale, dir_z / scale);

  LogDebug("SiPadDigitizerAlgorithm") << " The drift direction in local coordinate is " << theDriftDirection;
  return theDriftDirection;
}


// For premixing
void SiPadDigitizerAlgorithm::loadAccumulator(unsigned int detId, const std::map<int, float>& accumulator) {
  auto& theSignal = _signal[detId];
  for (const auto& elem : accumulator) {
    auto inserted = theSignal.emplace(elem.first, CommonDigiUtility::Amplitude(elem.second, nullptr));
    if (!inserted.second) {
      throw cms::Exception("LogicError") << "Signal was already set for DetId " << detId;
    }
  }
}


void SiPadDigitizerAlgorithm::GetDigiResults(const FbcmSiPadGeom* SiPadGeom,  std::map<int, SiPadDigiData>& SiPadDigilMap) {

uint32_t detID = SiPadGeom->geographicalId().rawId();
FbcmDetId SiPadDetId(detID);
  auto it = _signal.find(detID);
  if (it == _signal.end())
    return;

   if (addNoise)
    add_noise(SiPadGeom);  // generate noise

  const signal_map_type& theSignal = _signal[detID];
	HitAnalysisInfo HitTotToaInfo;
	std::vector<std::pair<float, CommonDigiUtility::PSimHitInfo> > BXC_CahrgePSim_Vect;
	std::vector< TofChargePair > Tof_Q_pairVect;
	std::pair<float, float> SiPadDimension;
	std::pair<float, const edm::ParameterSet * > Area_FeParamPtr;
	std::vector < HitAnalysisInfo > HitAnalysisVect;
  for (auto const& s : theSignal) // loop over channels ?? // by default one SiPad has one channel 
  {
	const CommonDigiUtility::Amplitude& sig_data = s.second;
	float signalInElectrons = sig_data.ampl()*chargeCollectionEff;
	BXC_CahrgePSim_Vect.clear();
	Tof_Q_pairVect.clear();
	for (auto const& l : sig_data.simInfoList()) {
		float PsimCharge_= l.first*chargeCollectionEff;
		BXC_CahrgePSim_Vect.push_back({PsimCharge_,*(l.second)}); // first: charge, second: PSimHit
		Tof_Q_pairVect.emplace_back(std::make_pair( l.second->time() , PsimCharge_ ));
        }
		
	SiPadDimension = SiPadGeom->SiPadTopology().pitch();
	float SiPadArea = SiPadDimension.first * SiPadDimension.second ;
	float RRadius = SiPadGeom->surface().position().perp();
	float PhiDegrees = SiPadGeom->surface().position().phi().degrees();
		
	HitPulse.GetPulseSeriesShape(FftPrep, Tof_Q_pairVect); // vector for charge amplitude
	Area_FeParamPtr = FeParamSelector.SelectFrontEndConfig(SiPadArea);
	FrontEnd.RunFECircuit(Area_FeParamPtr);
	
	HitAnalysisVect.clear();
	for (int BxSlotNo = FirstBxSlotNo ; BxSlotNo <= LastBxSlotNo ; BxSlotNo++) {
		HitTotToaInfo.clear();
		FrontEnd.GetHitAnalysisInfo(BxSlotNo, HitTotToaInfo);
		HitAnalysisVect.emplace_back(HitTotToaInfo);
		// if (HitTotToaInfo.nbrOfRecognizedHitsInBx() >= 1)  {
			// FrontEnd.printInfo_with_AlignedTime();
			// std::cin >> Teststp;						
		// } 
	}
	
	SiPadDigiData SiPadDigiRes( SiPadDetId.Side(),
								SiPadDetId.Station(),
								SiPadDetId.SiliconDie(),
								SiPadDetId.SiPad(),
								RRadius,
								PhiDegrees,
								SiPadArea,
								signalInElectrons  ,
								BXC_CahrgePSim_Vect,
								HitAnalysisVect);
	//std::cout << SiPadDigiRes;

	SiPadDigilMap.insert({s.first,SiPadDigiRes});
}

}

bool SiPadDigitizerAlgorithm::FilterHit(const PSimHit& hit, double tCorr) {
	double toa = hit.tof() - tCorr;
	//std::cout << "hitSelectionMode:" << hitSelectionMode_ << "\n";
	if (hitSelectionMode_ == 1) // all hits
		return true;
	else if (hitSelectionMode_ == 0) // filter with LowerToFCut and UpperToFCut
		return (toa >= theTofLowerCut && toa <= theTofUpperCut);
	else
		throw cms::Exception("Wrong hitSelectionMode") << "hitSelectionMode should be either '1' (select all hits) or '0' (filter with LowerToFCut and UpperToFCut) \n";	
}
