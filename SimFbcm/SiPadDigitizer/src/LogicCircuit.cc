///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimFbcm/SiPadDigitizer/interface/LogicCircuit.h"

namespace FbcmFE {

	LogicCircuit::LogicCircuit(LogicSignalType & RisingEdgeClkInput, LogicSignalType & ActiveLowReset, LogicSignalType & OutputSignal):
		RisingEdgeClkInput_(RisingEdgeClkInput),
		ActiveLowReset_(ActiveLowReset),
		OutputSignal_(OutputSignal)
		{ 		};

	LogicCircuit::~LogicCircuit(){};

    void LogicCircuit::RunLogic(){

    bool lastClk=LogicLevel::HIGH;
    bool clk;
    bool out=LogicLevel::LOW;

    //int start, end_;
    for (unsigned int i=0; i < RisingEdgeClkInput_.size() ; i++) {
            clk=RisingEdgeClkInput_[i];
            if (ActiveLowReset_[i]==LogicLevel::LOW)
                out=LogicLevel::LOW;
            else{
                if (clk!=lastClk && clk==LogicLevel::HIGH)
                    out=LogicLevel::HIGH;
            }
            lastClk=clk;
            OutputSignal_[i]=out;
    }
    }
}


