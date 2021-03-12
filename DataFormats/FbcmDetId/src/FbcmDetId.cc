
///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"



FbcmDetId::FbcmDetId() : DetId(DetId::BRIL, FbcmSubdetId::FbcmModule)
{}

FbcmDetId::FbcmDetId(uint32_t id) : DetId(id) {
  if (det() != DetId::BRIL || subdetId() != FbcmSubdetId::FbcmModule) {
    throw cms::Exception("InvalidDetId") << "FbcmDetId ctor:"
                                         << " det: " << det() << " subdet: " << subdetId() << " is not a valid FbcmDetId id";
  }
}

FbcmDetId::FbcmDetId(DetId id) : DetId(id) {
  if (det() != DetId::BRIL || subdetId() != FbcmSubdetId::FbcmModule) {
    throw cms::Exception("InvalidDetId") << "FbcmDetId ctor:"
                                         << " det: " << det() << " subdet: " << subdetId() << " is not a valid FbcmDetId id";
  }
}

FbcmDetId::FbcmDetId(int Side, int Station, int SiliconDie, int SiPad) : DetId(DetId::BRIL, FbcmSubdetId::FbcmModule) {
  this->init(Side, Station, SiliconDie, SiPad);
}

void FbcmDetId::init(int Side, int Station, int SiliconDie, int SiPad) {
  if ( Side < 1 || Side > maxSideId || Station < 0 || Station > maxStationId ||
      SiliconDie < 0 || SiliconDie > maxSiliconDieId || SiPad < 0 || SiPad > maxSiPadId) {
    throw cms::Exception("InvalidDetId") << "FbcmDetId ctor:"
                                         << " Invalid parameters: "
                                         " Side " << Side << " Station " << Station << " SiliconDie " << SiliconDie
                                         << " SiPad " << SiPad << std::endl;
  }
  

  id_ |= (Side & SideMask_) << SideStartBit_ | (Station & StationMask_) << StationStartBit_ |
         (SiliconDie & SiliconDieMask_) << SiliconDieStartBit_ | (SiPad & SiPadMask_) << SiPadStartBit_;
		 
	// id_	 has been defined in DetId class!
}

std::ostream& operator<<(std::ostream& os, const FbcmDetId& id) {
  os << " Side " << id.Side() << " Station " << id.Station() << " SiliconDie " << id.SiliconDie() << " SiPad "
     << id.SiPad() << " \n";

  return os;
}
