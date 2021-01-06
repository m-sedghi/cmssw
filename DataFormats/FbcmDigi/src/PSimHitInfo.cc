///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "DataFormats/FbcmDigi/interface/PSimHitInfo.h"
namespace CommonDigiUtility {

    PSimHitInfo::PSimHitInfo():
		eventId_(),
		trackId_(0),
		hitIndex_(0),
		tofBin_(0),
		time_(0.0),
		thePabs(0.0),
		theTof(0.0),
		theEnergyLoss(0.0),
		theParticleType(0),
		theDetUnitId(0),
		theProcessType(0)
		{}

std::ostream& operator<<(std::ostream& o, const PSimHitInfo& Ps) {

  o << "RawId:" << Ps.DetUnitId() << ", " 
	<< "Time:" << Ps.time() << ", " 
	<< "ToF:" << Ps.Tof() << ", " 
	<< "BxNo:" << Ps.BunchCrossing() << ", " 
	<< "Pabs:" << Ps.Pabs() << ", " 
	<< "EnLoss:" << Ps.EnergyLoss() << ", " 
	<< "PdgId:" << Ps.ParticleType() << ", " 
	<< "trackId:" << Ps.trackId() << ", " 
	<< "ProcessType:" << Ps.ProcessType() << "\n" ;

  return o ;
}

 std::ostream& operator<<(std::ostream& o, const CahrgePsimVect& ChP) {
for (auto d : ChP) {
	o << "\tCharge: " << d.first << ", ";
	o << d.second;
}
  return o ;
}



}