///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_PSimHitInfo_h
#define _FbcmFE_PSimHitInfo_h

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

#include <iostream>
#include <cstdint>
#include <utility>
#include <vector>

namespace CommonDigiUtility {

 class PSimHitInfo {
  public:
    explicit PSimHitInfo(const PSimHit* hitp, float corrTime, size_t hitIndex, uint32_t tofBin):
		eventId_(hitp->eventId()),
		trackId_(hitp->trackId()),
		hitIndex_(hitIndex),
		tofBin_(tofBin),
		time_(corrTime),
		thePabs(hitp->pabs()),
		theTof(hitp->tof()),
		theEnergyLoss(hitp->energyLoss()),
		theParticleType(hitp->particleType()),
		theDetUnitId(hitp->detUnitId()),
		theProcessType(hitp->processType())
		{}

	PSimHitInfo();
	~PSimHitInfo(){};
	
	EncodedEventId eventId() const { return eventId_; };
	uint32_t trackId() const { return trackId_; };
    uint32_t hitIndex() const { return hitIndex_; };
    uint32_t tofBin() const { return tofBin_; };
    float time() const { return time_; };
	float Pabs() const { return thePabs; };
	float Tof() const { return theTof; };
	float EnergyLoss() const { return theEnergyLoss; }; 
	int ParticleType() const { return theParticleType; };
	uint32_t DetUnitId() const { return theDetUnitId; };
	unsigned short ProcessType() const { return theProcessType; };
	int BunchCrossing() const { return eventId_.bunchCrossing(); };

	inline bool operator<(const PSimHitInfo& other) const { return BunchCrossing() < other.BunchCrossing(); }

  private:
	EncodedEventId eventId_;
    uint32_t trackId_;
    uint32_t hitIndex_; 
    uint32_t tofBin_;
    float time_;
	float thePabs;
	float theTof;
	float theEnergyLoss; 
	int theParticleType;
	uint32_t theDetUnitId;
	unsigned short theProcessType;
	//int TheBunchCrossing;
  };


typedef std::vector<std::pair<float,PSimHitInfo> > CahrgePsimVect;  

std::ostream& operator<<(std::ostream& o, const PSimHitInfo& Ps) ;
std::ostream& operator<<(std::ostream& o, const CahrgePsimVect& ChP);

}

#endif
