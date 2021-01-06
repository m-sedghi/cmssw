///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/LinearCircuitFreqDomain.h"
#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

	LinearCircuitFreqDomain::LinearCircuitFreqDomain(SignalType & InputSignal, SignalType & OutputSignal, Filter & FilterTF):
		InputSignal_(InputSignal),
		OutputSignal_(OutputSignal),
		FilterTF_(FilterTF){	};

	LinearCircuitFreqDomain::~LinearCircuitFreqDomain(){};

	void LinearCircuitFreqDomain::CalculateOutputSignal(){
		InternalSig.init(InputSignal_);
		InternalSig.ApplyFilter(FilterTF_);
		InternalSig.GetRealSignal(OutputSignal_);
	}

}


