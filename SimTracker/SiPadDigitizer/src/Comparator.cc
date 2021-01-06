///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/Comparator.h"

namespace FbcmFE {


    Comparator::Comparator(SignalType & InputSignal, LogicSignalType & OutputSignal):
		InputSignal_(InputSignal),
		OutputSignal_(OutputSignal),
		CompStatus(CompLastStatus::UNKNOWN)
		{ };

	Comparator::~Comparator(){};

    void Comparator::RunComparator(){
		bool OuputLogic=0;
        for (unsigned int i=0 ; i < InputSignal_.size() ; i++) {
            switch (CompStatus) {
            case LOW:
                if (InputSignal_[i] >= UpperTsh) {
                    OuputLogic=HIGH;
                    CompStatus=HIGH;
                }
                else {
                    OuputLogic=LOW;
                    CompStatus=LOW;
                }
                break;

            case HIGH:
                if (InputSignal_[i] <= LowerTsh) {
                    OuputLogic=LOW;
                    CompStatus=LOW;
                }
                else {
                    OuputLogic=HIGH;
                    CompStatus=HIGH;
                }
                break;
            case UNKNOWN:
                if (InputSignal_[i] >= UpperTsh) {
                    OuputLogic=HIGH;
                    CompStatus=HIGH;
                }
                else if (InputSignal_[i] <= LowerTsh) {
                    OuputLogic=LOW;
                    CompStatus=LOW;
                }
                else {
                    // assuming low by default for no signal
                    // it could be "high" for a different circuit
                    OuputLogic=LOW;
                    CompStatus=LOW;
                }
                break;
            default:
                break;
            }
            OutputSignal_[i]=OuputLogic;
    }
}
}
