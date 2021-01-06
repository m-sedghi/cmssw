///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_CFDCore_h
#define _FbcmFE_CFDCore_h

#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

	class CFDCore {
	public:
		CFDCore(SignalType & InputSignalOrg, SignalType & SigDelayedInput, SignalType & OutputSignal);
		~CFDCore();
		void RunCFD();
		void SetFraction(double cfd_fraction) { Fraction = cfd_fraction; }

	private:
		SignalType & InputSignalOrg_;
		SignalType & InputSignalDelayed_;
		SignalType & OutputSignal_;
		double Fraction;
	};

}
#endif
