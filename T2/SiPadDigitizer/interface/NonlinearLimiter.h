///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_NonlinearLimiter_h
#define _FbcmFE_NonlinearLimiter_h

#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"


namespace FbcmFE {

class NonlinearLimiter {
	
	public:
		NonlinearLimiter(SignalType & InputSignal, SignalType & OutputSignal);
		~NonlinearLimiter();
		void SetAbsMaxOut(double MaxLimAbs_) { MaxLimAbs = MaxLimAbs_ ;}
		void RunLimiter();
	private:
		void Clamp(SignalType & Sig_, double lowerlimit, double upperlimit);
		
		SignalType & InputSignal_;
		SignalType & OutputSignal_;
		double MaxLimAbs; // mV
	};

}
#endif
