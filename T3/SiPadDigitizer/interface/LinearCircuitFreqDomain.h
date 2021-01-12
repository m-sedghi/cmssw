///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_LinearCircuitFreqDomain_h
#define _FbcmFE_LinearCircuitFreqDomain_h

#include "SimFbcm/SiPadDigitizer/interface/Complx.h"
#include "SimFbcm/SiPadDigitizer/interface/FourierTransformableSignal.h"
#include "SimFbcm/SiPadDigitizer/interface/Filter.h"
#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

	class LinearCircuitFreqDomain {
	public:

		LinearCircuitFreqDomain(SignalType & InputSignal, SignalType & OutputSignal, Filter & FilterTF);
		~LinearCircuitFreqDomain();
		void CalculateOutputSignal();

	private:
		SignalType & InputSignal_;
		SignalType & OutputSignal_;
		Filter & FilterTF_;
		FourierTransformableSignal InternalSig;
	};

}
#endif
