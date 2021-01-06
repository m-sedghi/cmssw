///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_FbcmFrontEndChip_h
#define _FbcmFE_FbcmFrontEndChip_h

#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"
#include "SimTracker/SiPadDigitizer/interface/FftPreparation.h"
#include "SimTracker/SiPadDigitizer/interface/Complx.h"
#include "SimTracker/SiPadDigitizer/interface/FourierTransformableSignal.h"
#include "SimTracker/SiPadDigitizer/interface/Filter.h"
#include "SimTracker/SiPadDigitizer/interface/LinearCircuitFreqDomain.h"
#include "SimTracker/SiPadDigitizer/interface/NonlinearLimiter.h"
#include "SimTracker/SiPadDigitizer/interface/CFDCore.h"
#include "SimTracker/SiPadDigitizer/interface/Comparator.h"
#include "SimTracker/SiPadDigitizer/interface/LogicCircuit.h"
#include "SimTracker/SiPadDigitizer/interface/HitAnalyzer.h"
//#include "SimTracker/SiPadDigitizer/interface/HitAnalysisInfo.h"
#include "DataFormats/FbcmDigi/interface/HitAnalysisInfo.h"

#include <fstream>

namespace FbcmFE {

    class FbcmFrontEndChip  {
    public:
        FbcmFrontEndChip(FftPreparation & FFtPrep); 

        ~FbcmFrontEndChip();
        void printInfo();
        void printInfo_with_AlignedTime();
        void PrintInputCurrent_uA();
        void Print_TIAOutput_Voltage_mV();
        void Print_ShaperOutput_Voltage_mV();
        void PrintTiaOutput();

        void RunFECircuit( std::pair<float, const edm::ParameterSet * > Area_FeParamPtr );
        void GetHitAnalysisInfo(int BXC_SlotNo, HitAnalysisInfo & HitTotToaInfo);

    private:

		void SetSubmoduleParameters();
        void Set_TIA_InputSignal(SignalType & PulseShape);

        FftPreparation & fftPrep;
		const edm::ParameterSet * ActiveFrontEndParamPtr;
		float SensorArea;
		unsigned int nFFT_;
		SignalType & freqVect_GHz;
		SignalType SiPadReceivedSig_;
		SignalType TIA_OutputSig_;
		SignalType Limiter_OutputSig_;
		SignalType Shaper_OutputSig_;
		SignalType DelayedSig_;
		SignalType DiscrimSig_;
		SignalType DiscrimShaperOutSig_;
		LogicSignalType CFD_ZeroCCompOutSig_;
		LogicSignalType CFD_ArmingCompOutSig_;
		LogicSignalType CFD_LogicOutput_;
		Filter TIA_Hf_;
		Filter Shaper_Hf_;
		Filter Delay_Hf_;
		Filter CFD_ShaperHf_;
		LinearCircuitFreqDomain TIA_;
		NonlinearLimiter Limiter_;
		LinearCircuitFreqDomain Shaper_;
		LinearCircuitFreqDomain DelayCircuit_;
		CFDCore Cfd_;
		LinearCircuitFreqDomain CFD_Shaper_;
		Comparator ZC_Comp;
		Comparator Arming_Comp;
		LogicCircuit OutputLogicCirc;
		HitAnalyzer Hit_Analyzer;
		
		
};


}
#endif
