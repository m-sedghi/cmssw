///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_Comparator_h
#define _FbcmFE_Comparator_h

#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

    enum CompLastStatus	{
		LOW,
		HIGH,
		UNKNOWN
	};

	class Comparator {


	public:
		Comparator(SignalType & InputSignal, LogicSignalType & OutputSignal);
		~Comparator();
		void RunComparator();
		void SetComparatorThresholds(double lowerTSH_, double upperTSH_) { LowerTsh = lowerTSH_;
																		   UpperTsh = upperTSH_ ;}

	private:
		SignalType & InputSignal_;
		LogicSignalType & OutputSignal_;
		double LowerTsh; // the threshold that makes low-level signal
		double UpperTsh; // the threshold that makes high-level signal
		CompLastStatus CompStatus;

	};

}
#endif
