///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_Filter_h
#define _FbcmFE_Filter_h


#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "SimFbcm/SiPadDigitizer/interface/Complx.h"


namespace FbcmFE {

	const static double PI=3.14159265359;
	enum FilterType
	{
		TIA_Filter,
		Shaper_Filter,
		TIA_Shaper_Filter,
		IdealDelay,
		FirstOrderPadeDelay,
		CDFShaper,
		NotSet
	};

	class Filter {
	public:
		Filter(SignalType & FreqGHz);
		~Filter();
		ComplexSignalType &  GetTransferFunction() {return Filter_Tf;}
		void RunFilter(FilterType type_);
		void SetParameters( const edm::ParameterSet * FEParamPtr );
		void SetSensorArea( float Area ){ SensorArea = Area; }

	private:
		void RunTIAFilter();
		void RunShaperFilter();
		void RunIdealDelayFilter();
		void Run1stPadeDelayFilter();
		void RunCDFShaper();
        FilterType Type_;
		ComplexSignalType Filter_Tf;
		SignalType & freq_GHz_ ;
		const edm::ParameterSet * FrontEndParamPtr; 
		
		double TIA_ShaperGain;
		double rf; // kOhm
		double cf; // pF
        double cgs; // pF
        double cc; // pF // coupling capacitance
		double C_per_cm2; // pF/cm2
		float SensorArea;
		double co; // pF
		double gin; // mS
		double Tau_sh; // ns //Shaper1 Tau
		double CDF_Delay; // ns
		double CFD_Gain; // gain [V/V]
		double Tau_CFDsh; // ns
		
	};

}
#endif

