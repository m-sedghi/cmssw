///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimFbcm/SiPadDigitizer/interface/Filter.h"

namespace FbcmFE {

	Filter::Filter(SignalType & freqGHz):Type_(FilterType::NotSet),freq_GHz_(freqGHz){
		SensorArea = 0.0;
	};

	Filter::~Filter(){};

	void Filter::RunFilter(FilterType type_){
		Type_=type_;
		switch (Type_){
		case FbcmFE::TIA_Filter:
			RunTIAFilter();
			break;
		case FbcmFE::Shaper_Filter:
			RunShaperFilter();
			break;
        case FbcmFE::IdealDelay:
            RunIdealDelayFilter();
			break;
		case FbcmFE::FirstOrderPadeDelay:
            Run1stPadeDelayFilter();
			break;
        case FbcmFE::CDFShaper:
            RunCDFShaper();
			break;
		default:
		throw cms::Exception("Wrong FilterType")
          << "FilterType is wrong or has not been set.\n";
			break;
		}
	}


	void Filter::RunTIAFilter(){
		Complx HfVal_tia;
		Complx d1,d2;
		double omega; // rad per ns = f_Ghz*2*Pi
        // current : uA
        // Voltage : mV
        //double TIA_ShaperGain=28;// 1.8*90.0;
		//double rf=5.0; // kOhm
		//double cf=0.25; // pF
        //double cgs=0.4; // pF
        //double cc=315.0 ; // pF // coupling capacitance
		double cd = SensorArea*C_per_cm2;
		double cin_=cgs+(cc*cd/(cd+cc)); // pF
		double i_norton_factor=1.0/(1.0+(cd/cc));
		//double co=0.4; // pF
		//double gin=3; // mS
		double Tau1=rf*cf;
		double Tau2=(cin_*co)/(gin*cf);
        Filter_Tf.clear();
		for (auto f_ghz: freq_GHz_) {
			omega=2*f_ghz*(FbcmFE::PI);
			d1.Set(1,omega*Tau1);
			d2.Set(1,omega*Tau2);
			HfVal_tia.Set(d1*d2);
			HfVal_tia.Reciprocal();
			HfVal_tia=HfVal_tia*rf*TIA_ShaperGain*i_norton_factor;
			Filter_Tf.emplace_back(HfVal_tia);
		}
		
		//std::cout <<"Area: " << SensorArea << ", Cd: " << cd << ", " << "gain: " << TIA_ShaperGain << "\n"  ;
	}


	void Filter::RunShaperFilter(){

		Complx HfVal_Shap;
		Complx d3;
        double Shaper1Gain=1.0;; // gain [V/V]
		double omega; // rad per ns = f_Ghz*2*Pi
        //double Tau_sh= 0.9; // ns
		Filter_Tf.clear();
		for (auto f_ghz: freq_GHz_) {
			omega=2*f_ghz*(FbcmFE::PI);
			d3.Set(1,omega*Tau_sh);
			HfVal_Shap.Set(d3*d3);
			HfVal_Shap.Reciprocal();
			HfVal_Shap=HfVal_Shap*Shaper1Gain;
			Filter_Tf.emplace_back(HfVal_Shap);
		}
	}

    void Filter::Run1stPadeDelayFilter(){
        Complx HfVal;
		Complx n1,d1;
		//double CDF_Delay=2.0; // ns
        double omega; // rad per ns = f_Ghz*2*Pi
        Filter_Tf.clear();
		for (auto f_ghz: freq_GHz_) {
			omega=2*f_ghz*(FbcmFE::PI);
			n1.Set(2, omega*CDF_Delay);
			d1.Set(2, -omega*CDF_Delay);
			HfVal=n1/d1;
			Filter_Tf.emplace_back(HfVal);
		}
    }
    void Filter::RunIdealDelayFilter(){
        Complx HfVal;
		//double CDF_Delay=2.0; // ns
        double omega; // rad per ns = f_Ghz*2*Pi
        double re, im;
        Filter_Tf.clear();
		for (auto f_ghz: freq_GHz_) {
			omega=2*f_ghz*(FbcmFE::PI);
			re=cos(omega*CDF_Delay);
			im=-sin(omega*CDF_Delay);
			HfVal.Set(re,im);
			Filter_Tf.emplace_back(HfVal);
		}
    }

    void Filter::RunCDFShaper(){
    Complx HfVal_Shap;
		Complx d3;
        //double CFD_Gain=1.5; // gain [V/V]
		double omega; // rad per ns = f_Ghz*2*Pi
        //double Tau_CFDsh= 0.25; // ns
		Filter_Tf.clear();
		for (auto f_ghz: freq_GHz_) {
			omega=2*f_ghz*(FbcmFE::PI);
			d3.Set(1,omega*Tau_CFDsh);
			HfVal_Shap.Set(d3*d3);
			HfVal_Shap.Reciprocal();
			HfVal_Shap=HfVal_Shap*CFD_Gain;
			Filter_Tf.emplace_back(HfVal_Shap);
		}

    }
	
	
	void Filter::SetParameters( const edm::ParameterSet * FEParamPtr ){ 
	FrontEndParamPtr = FEParamPtr;
	TIA_ShaperGain = FrontEndParamPtr->getParameter< double >("TIA_Shaper_Gain"); // gain [V/V]
	rf = FrontEndParamPtr->getParameter< double >("Tia_Rf"); // kOhm
	cf = FrontEndParamPtr->getParameter< double >("Tia_Cf"); // pF
	cgs = FrontEndParamPtr->getParameter< double >("Tia_Cin_gs"); // pF
	cc = FrontEndParamPtr->getParameter< double >("SensorCouplingCapacitance"); // pF // coupling capacitance
	C_per_cm2 = FrontEndParamPtr->getParameter< double >("SensorCapPerCm2"); // pF/cm2
	co = FrontEndParamPtr->getParameter< double >("Tia_Co"); // pF
	gin = FrontEndParamPtr->getParameter< double >("Tia_gin"); // mS
	Tau_sh = FrontEndParamPtr->getParameter< double >("Shaper1_Tau"); // ns //Shaper1 Tau
	CDF_Delay = FrontEndParamPtr->getParameter< double >("CFD_Delay"); // ns
	CFD_Gain = FrontEndParamPtr->getParameter< double >("CfdShaper_Gain"); // gain [V/V]
	Tau_CFDsh = FrontEndParamPtr->getParameter< double >("CfdShaper_Tau"); // ns

// std::cout << "TIA_Shaper_Gain: " << TIA_ShaperGain << "\n"
// << "Tia_Rf: " << rf << "\n"
// << "Tia_Cf: " << cf << "\n"
// << "Tia_Cin_gs: " << cgs << "\n"
// << "SensorCouplingCapacitance: " << cc << "\n"
// << "SensorCapPerCm2: " << C_per_cm2 << "\n"
// << "Tia_Co: " << co << "\n"
// << "Tia_gin: " << gin << "\n"
// << "Shaper1_Tau: " << Tau_sh << "\n"
// << "CDF_Delay: " << CDF_Delay << "\n"
// << "CfdShaper_Gain: " << CFD_Gain << "\n"
// << "CfdShaper_Tau: " << Tau_CFDsh << "\n"
// << "SensorArea: " << SensorArea << "\n";

	}
	
}


