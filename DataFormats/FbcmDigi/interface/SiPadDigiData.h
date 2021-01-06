#ifndef DataFormats_FbcmDigi_SiPadDigiData_H
#define DataFormats_FbcmDigi_SiPadDigiData_H
///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <cstdint>
#include <utility>
#include <cassert>
#include "DataFormats/FbcmDigi/interface/PSimHitInfo.h"
#include "DataFormats/FbcmDigi/interface/HitAnalysisInfo.h"

using namespace CommonDigiUtility;
using namespace FbcmFE;
	
class SiPadDigiData {

public:
 explicit SiPadDigiData(
  uint8_t SideNo,
  uint8_t StationNo,
  uint8_t SiliconDieNo,
  uint16_t SiPadNo,
  float radius,
  float Phi_Deg,
  float SiPadArea,
  float AccumulatedCharge  ,
  std::vector<std::pair<float, PSimHitInfo> > CahrgePsimVect ,
  std::vector < HitAnalysisInfo > BxSlotHitAnalysisVect  
  ) :
		  SideNo_(SideNo),
		  StationNo_(StationNo),
		  SiliconDieNo_(SiliconDieNo),
		  SiPadNo_(SiPadNo),
		  Radius_(radius),
		  Phi_Deg_(Phi_Deg),
		  SiPadArea_(SiPadArea),
		  AccumulatedCharge_(AccumulatedCharge) ,
		  CahrgePsimVect_(CahrgePsimVect) ,
		  BxSlotHitAnalysisVect_(BxSlotHitAnalysisVect) 
		  {        };

  SiPadDigiData( ):
		  SideNo_(0),
		  StationNo_(0),
		  SiliconDieNo_(0),
		  SiPadNo_(0),
		  Radius_(0.0),
		  Phi_Deg_(0.0),
		  SiPadArea_(0.0),
		  AccumulatedCharge_(0.0) ,
		  CahrgePsimVect_() ,
		  BxSlotHitAnalysisVect_()
		  { 	 };

  ~SiPadDigiData(){};
  
  
  
  int SideIndex() const { return (int)SideNo_; }
  int StationIndex() const { return (int)StationNo_; }
  int SiliconDieIndex() const { return (int)SiliconDieNo_; }
  int SiPadIndex() const { return (int)SiPadNo_; }
  float Radius() const { return Radius_ ; }
  float Phi_Degrees() const { return Phi_Deg_ ; }
  float Area() const { return SiPadArea_ ; }
  float ChargeSum() const { return AccumulatedCharge_; }
  std::vector<std::pair<float, PSimHitInfo> > CahrgePsimVector() const { return CahrgePsimVect_; }
  std::vector < HitAnalysisInfo > BxSlotHitAnalysisVector() const { return BxSlotHitAnalysisVect_; }
  

  inline bool operator<(const SiPadDigiData& other) const { return ChargeSum() < other.ChargeSum(); }

private:
  uint8_t SideNo_;
  uint8_t StationNo_;
  uint8_t SiliconDieNo_;
  uint16_t SiPadNo_;
  float Radius_ ;
  float Phi_Deg_ ;
  float SiPadArea_ ;
  float AccumulatedCharge_;
  std::vector<std::pair<float, PSimHitInfo> > CahrgePsimVect_;
  std::vector < HitAnalysisInfo > BxSlotHitAnalysisVect_;
  
};

#include <iostream>
inline std::ostream& operator<<(std::ostream& o, const SiPadDigiData& D) {
	o << "=====================================================================================\n" ;
  o << "SideNo:" << (int) (D.SideIndex()) << ", " 
	<< "StationNo:" << (int) (D.StationIndex()) << ", " 
	<< "SiDieNo:" << (int) (D.SiliconDieIndex()) << ", " 
	<< "SiPadNo:" << (int) (D.SiPadIndex()) << ", " 
	<< "Radius:" << D.Radius() << ", " 
	<< "PhiDeg:" << D.Phi_Degrees() << ", " 
	<< "Area:" << D.Area() << ", " 
	<< "ChargeSum:" << D.ChargeSum()  << ", " 
	<< "QPsVect_Size:" << D.CahrgePsimVector().size()  << ", " 
	<< "BxHitAnalVect_Size:" << D.BxSlotHitAnalysisVector().size()  << "\n"; 
  o << D.CahrgePsimVector() ; 
  for (auto h : D.BxSlotHitAnalysisVector()) {
	    o << h;
  }
  

  return o ;
}

#endif