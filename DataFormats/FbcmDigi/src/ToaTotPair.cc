///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "DataFormats/FbcmDigi/interface/ToaTotPair.h"
namespace FbcmFE {
std::ostream& operator<<(std::ostream& o, ToaTotPair& TTpair) {
			return o 
			<< "\tToA: " << TTpair.ToA() << ", "
            << "\tToT: " << TTpair.ToT() << ", "
            << "\tToAToAStatus: " << TTpair.ToAToAStatus() << ", "
            << "\tSubBxBinNumber: " << TTpair.SubBxBinNumber() << ", "
            << "\tIsToTValid: " << TTpair.IsToTValid() << "\n" ;
}

}