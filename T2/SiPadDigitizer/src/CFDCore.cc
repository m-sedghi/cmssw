///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimFbcm/SiPadDigitizer/interface/CFDCore.h"

namespace FbcmFE {


	CFDCore::CFDCore(SignalType & InputSignalOrg, SignalType & SigDelayedInput, SignalType & OutputSignal):
		InputSignalOrg_(InputSignalOrg),
		InputSignalDelayed_(SigDelayedInput),
		OutputSignal_(OutputSignal)
		{ };

	CFDCore::~CFDCore(){};

    void CFDCore::RunCFD(){
        OutputSignal_=InputSignalDelayed_-(Fraction*InputSignalOrg_);
	}

}


