///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_ToaTotPair_h
#define _FbcmFE_ToaTotPair_h
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

namespace FbcmFE {

#define ToAStatus_BeforeBx_NonOverlapping   -2
#define ToAStatus_BeforeBx_Overlapping      -1
#define ToAStatus_FullyWithinBx              0
#define ToAStatus_WithinBx_LastsAfter        1
#define ToAStatus_AfterBx                    2
#define ToAStatus_UnKnown                    3

enum  ToAStatus : char
{
    BeforeBx_NonOverlapping  = ToAStatus_BeforeBx_NonOverlapping,
    BeforeBx_Overlapping     = ToAStatus_BeforeBx_Overlapping,
    FullyWithinBx            = ToAStatus_FullyWithinBx,
    WithinBx_LastsAfter      = ToAStatus_WithinBx_LastsAfter,
    AfterBx                  = ToAStatus_AfterBx,
    UnKnown                  = ToAStatus_UnKnown
};

class ToaTotPair {
public:
    explicit ToaTotPair():
        TimeOfArrival(0.0),
        TimeOverThreshold(0.0),
        IsToTValid_(true),
        //ToA_Status_(ToAStatus::UnKnown),
		ToA_Status_(ToAStatus_UnKnown),
        SubBxBinNo_(0)   {};

    ~ToaTotPair() {};

    void SetPairInfo(float ToA_,float ToT_, bool IsToTValid, ToAStatus ToA_Status, int16_t SubBxBinNo) {
        TimeOfArrival=ToA_;
        TimeOverThreshold=ToT_;
        IsToTValid_=IsToTValid;
        ToA_Status_=(char)ToA_Status;
        SubBxBinNo_=SubBxBinNo; }

    float ToA() const {return TimeOfArrival;}
    float ToT() const {return TimeOverThreshold;}
    ToAStatus ToAToAStatus() const {return (ToAStatus)ToA_Status_;}
	int ToAToAStatusInt() const {return (int)ToA_Status_;}
    int SubBxBinNumber() const {return SubBxBinNo_;}
    bool IsToTValid() {return IsToTValid_;}
	
	inline bool operator<(const ToaTotPair& other) const { return ToA() < other.ToA(); }
	
private:
    float TimeOfArrival;
    float TimeOverThreshold;
    bool IsToTValid_; // if ToT is valid over the time span of the simulator, this flag will rise up. otherwise, this means that ToT is too long and exceeds the time span of the simulator.
    //ToAStatus ToA_Status_;
	char ToA_Status_;
    int16_t SubBxBinNo_;
};

std::ostream& operator<<(std::ostream& o, ToaTotPair& TTpair);


}

#endif
