///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_LogicCircuit_h
#define _FbcmFE_LogicCircuit_h

#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

	class LogicCircuit {

    enum LogicLevel
	{
		LOW,
		HIGH,
		UNKNOWN
	};


	public:
		LogicCircuit(LogicSignalType & RisingEdgeClkInput, LogicSignalType & ActiveLowReset, LogicSignalType & OutputSignal);
		~LogicCircuit();
		void RunLogic();

	private:
		LogicSignalType & RisingEdgeClkInput_;
		LogicSignalType & ActiveLowReset_;
		LogicSignalType & OutputSignal_;
	};

}
#endif
