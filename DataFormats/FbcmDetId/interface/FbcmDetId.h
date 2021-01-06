#ifndef DataFormats_FbcmDetId_h
#define DataFormats_FbcmDetId_h

///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <iosfwd>
#include <iostream>


class FbcmSubdetId {
public:
  static constexpr int FbcmModule = 1;
  // for future usage!
  //very important: the value of subdet should be excatlcy eaqual to the FBCM enum in subDetId[20] in
  // Geometry/CommonDetUnit/interface/GeomDetEnumerators.h
  
};



class FbcmDetId : public DetId {
public:
  FbcmDetId();

  FbcmDetId(uint32_t id);
  FbcmDetId(DetId id);

  /// Construct from fully qualified identifier.
  FbcmDetId(int Side, int Station, int SiliconDie, int SiPad);

 
  /// Side id: 1 or 2 for +/- Z 
  unsigned int Side() const { return  int((id_ >> SideStartBit_) & SideMask_); }
  //unsigned int Side() const { return ((int(id_ & CopyNoMask)/100)%10); }

  /// Station id: each Side has four Station: 0-3
  unsigned int Station() const { return ((id_ >> StationStartBit_) & StationMask_); }
  //unsigned int Station() const { return ((int(id_ & CopyNoMask)/10)%10); }

  /// SiliconDie id: it identifies a SiliconDie in a Station: generally each station has one SiliconDie: Zero by default
  unsigned int SiliconDie() const { return ((id_ >> SiliconDieStartBit_) & SiliconDieMask_); }
  //unsigned int SiliconDie() const { return 0 ; }

  /// SiPad 
  unsigned int SiPad() const {  return ((id_ >> SiPadStartBit_) & SiPadMask_);   }
  //unsigned int SiPad() const {    return ((int(id_ & CopyNoMask))%10);    }
  
  /// Return the corresponding StationId (mask SiliconDieIds)
  FbcmDetId StationDetId() const { return FbcmDetId(id_ & StationIdMask_); }
  /// Return the corresponding SiliconDieId (mask SiPadIds)
  FbcmDetId SiliconDieDetId() const { return FbcmDetId(id_ & SiliconDieIdMask_); }

 

	static const int maxSideId = 2; // only 1 or 2 // Zero is not valid
	static const int maxStationId = 63; // 0-63
	static const int maxSiliconDieId = 127; //0-127 
	static const int maxSiPadId = 1023; // 0-1023

private:
/**
The DetId is a 32-bit unsigned integer.
The four most significant bits (calc) identify the large-scale detector (e.g. Tracker or Ecal)
The next three bits ([27:25]) identify a part of the detector (such as HcalBarrel (HB) for Hcal).
*/


  static const unsigned int SideStartBit_ = 23;
  static const unsigned int StationStartBit_ = 17;
  static const unsigned int SiliconDieStartBit_ = 10;
  static const unsigned int SiPadStartBit_ = 0;
  

  static const unsigned int SideMask_ = 0x3; // 2 bits
  static const unsigned int StationMask_ = 0x3F; // 6 bits
  static const unsigned int SiliconDieMask_ = 0x7F; // 7 bits
  static const unsigned int SiPadMask_ = 0x3FF; // 10 bits
  
   static const uint32_t StationIdMask_ = ~((SiliconDieMask_ << SiliconDieStartBit_) | (SiPadMask_ << SiPadStartBit_));
  static const uint32_t SiliconDieIdMask_ = ~(SiPadMask_ << SiPadStartBit_);


private:
  void init(int Side, int Station, int SiliconDie, int SiPad);
  
};  // FbcmDetId

//inline std::ostream& operator<<(std::ostream& os, const FbcmDetId& id) { return os << " " << id.SiPad(); }
std::ostream& operator<<(std::ostream& os, const FbcmDetId& id);

#endif
