///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#include "SimFbcm/SiPadDigitizer/interface/FbcmFrontEndChip.h"

namespace FbcmFE {

FbcmFrontEndChip::FbcmFrontEndChip(FftPreparation & FFtPrep):
        fftPrep(FFtPrep),
		nFFT_(fftPrep.nbrFFTpoints()),
        freqVect_GHz(FFtPrep.FreqVectGHz()),
		SiPadReceivedSig_(nFFT_),
		TIA_OutputSig_(nFFT_),
		Limiter_OutputSig_(nFFT_),
		Shaper_OutputSig_(nFFT_),
		DelayedSig_(nFFT_),
        DiscrimSig_(nFFT_),
        DiscrimShaperOutSig_(nFFT_),
        CFD_ZeroCCompOutSig_(nFFT_),
		CFD_ArmingCompOutSig_(nFFT_),
		CFD_LogicOutput_(nFFT_),
		TIA_Hf_(freqVect_GHz),
		Shaper_Hf_(freqVect_GHz),
		Delay_Hf_(freqVect_GHz), // IdealDelay or FirstOrderPadeDelay
		CFD_ShaperHf_(freqVect_GHz),
		TIA_(SiPadReceivedSig_,TIA_OutputSig_,TIA_Hf_),
		Limiter_(TIA_OutputSig_,Limiter_OutputSig_),
		Shaper_(Limiter_OutputSig_,Shaper_OutputSig_,Shaper_Hf_),
		DelayCircuit_(Shaper_OutputSig_,DelayedSig_,Delay_Hf_),
		Cfd_(Shaper_OutputSig_,DelayedSig_,DiscrimSig_),
		CFD_Shaper_(DiscrimSig_,DiscrimShaperOutSig_,CFD_ShaperHf_),
		ZC_Comp(DiscrimShaperOutSig_,CFD_ZeroCCompOutSig_),
		Arming_Comp(Shaper_OutputSig_,CFD_ArmingCompOutSig_),
		OutputLogicCirc(CFD_ZeroCCompOutSig_,CFD_ArmingCompOutSig_,CFD_LogicOutput_),
		Hit_Analyzer(FFtPrep, CFD_LogicOutput_)
    {

} ;
	FbcmFrontEndChip::~FbcmFrontEndChip(){};

	void FbcmFrontEndChip::printInfo(){
      std::ofstream myfile;
      myfile.open ("Result.txt");

        for (unsigned int i=0; i < nFFT_; i++ )
            myfile << SiPadReceivedSig_[i] << "\t"
                    << TIA_OutputSig_[i] << "\t"
                    << Shaper_OutputSig_[i] << "\t"
                    << DelayedSig_[i] << "\t"
                    << DiscrimShaperOutSig_[i] << "\t"
                    << CFD_ZeroCCompOutSig_[i] << "\t"
                    << CFD_ArmingCompOutSig_[i] << "\t"
                    << CFD_LogicOutput_[i] << "\n";

        myfile.close();
        std::cout << "done! \n";
	}



	void FbcmFrontEndChip::printInfo_with_AlignedTime(){
       SignalType & timeVectAligned_=Hit_Analyzer.GetAlignedTimeVect();
        std::ofstream myfile;
      myfile.open ("/afs/cern.ch/work/m/msedghi/CMSSW_11_0_0/src/SimTracker/SiPadDigitizer/Result.txt");

        for (unsigned int i=0; i < nFFT_; i++ )
            myfile << timeVectAligned_[i] << "\t"
                    << SiPadReceivedSig_[i] << "\t"
                    << TIA_OutputSig_[i] << "\t"
                    << Shaper_OutputSig_[i] << "\t"
                    << DelayedSig_[i] << "\t"
                    << DiscrimShaperOutSig_[i] << "\t"
                    << CFD_ZeroCCompOutSig_[i] << "\t"
                    << CFD_ArmingCompOutSig_[i] << "\t"
                    << CFD_LogicOutput_[i] << "\n";

        myfile.close();
        std::cout << "done! \n";
	}


    void FbcmFrontEndChip::Print_ShaperOutput_Voltage_mV(){
		for (auto i:Shaper_OutputSig_)
			 std::cout << i << "\n";
    }

    void FbcmFrontEndChip::Print_TIAOutput_Voltage_mV(){
		for (auto i:TIA_OutputSig_)
			 std::cout << i << "\n";
    }

	void FbcmFrontEndChip::PrintInputCurrent_uA(){
		for (auto i:SiPadReceivedSig_)
			 std::cout << i << "\n";
    }



	void FbcmFrontEndChip::Set_TIA_InputSignal(SignalType & PulseShape){
	    double Dimention=e_q/_ns/_uA; // the current will be in micro Ampere (uA)
		for (unsigned int i=0; i< nFFT_ ; i++ )
			SiPadReceivedSig_[i]=PulseShape[i]*Dimention; // uA
	}

	void FbcmFrontEndChip::PrintTiaOutput(){

	    double max_=0.0;
	    for (uint16_t i=0 ; i< nFFT_ ; i++ )
            if (Shaper_OutputSig_[i] > max_)
                max_=Shaper_OutputSig_[i];
        std::cout << max_ << "\n";


	}

	void FbcmFrontEndChip::SetSubmoduleParameters() {
		
	TIA_Hf_.SetParameters( ActiveFrontEndParamPtr ) ;
	Shaper_Hf_.SetParameters( ActiveFrontEndParamPtr ) ;
	Delay_Hf_.SetParameters( ActiveFrontEndParamPtr ) ;
	CFD_ShaperHf_.SetParameters( ActiveFrontEndParamPtr ) ;
	Limiter_.SetAbsMaxOut(ActiveFrontEndParamPtr->getParameter< double >("MaxFEOutputVoltage"));
	Cfd_.SetFraction(ActiveFrontEndParamPtr->getParameter< double >("CFD_Fraction"));
	//ZC_Comp.SetComparatorThresholds(-5.0,0.0);
	ZC_Comp.SetComparatorThresholds(ActiveFrontEndParamPtr->getParameter< double >("ZCComp_LowerTsh"),
									ActiveFrontEndParamPtr->getParameter< double >("ZCComp_UpperTsh"));
	
	//Arming_Comp.SetComparatorThresholds(5.0,20.0);
	Arming_Comp.SetComparatorThresholds(ActiveFrontEndParamPtr->getParameter< double >("ArmingComp_LowerTsh"),
										ActiveFrontEndParamPtr->getParameter< double >("ArmingComp_UpperTsh"));

	Hit_Analyzer.SetParameters( ActiveFrontEndParamPtr ) ;
	
	}

	void FbcmFrontEndChip::RunFECircuit(std::pair<float, const edm::ParameterSet * > Area_FeParamPtr){
		
		SensorArea =  Area_FeParamPtr.first;
		ActiveFrontEndParamPtr=Area_FeParamPtr.second;
		TIA_Hf_.SetSensorArea( SensorArea ); 
		SetSubmoduleParameters();
	
		Set_TIA_InputSignal(fftPrep.TimeDomainSignal());
		TIA_Hf_.RunFilter(FilterType::TIA_Filter);
		TIA_.CalculateOutputSignal();
		Limiter_.RunLimiter();
		Shaper_Hf_.RunFilter(FilterType::Shaper_Filter);
		Shaper_.CalculateOutputSignal();
		
		std::string DelayModel(ActiveFrontEndParamPtr->getParameter< std::string >("DelayModel"));
		FilterType DelayType;
		if (DelayModel == "FirstOrderPadeDelay")
			DelayType=FilterType::FirstOrderPadeDelay;
		else if ( DelayModel == "IdealDelay")
			DelayType=FilterType::IdealDelay;
		else
			throw cms::Exception("Wrong Delay type")
				<< "the delay type in CFD should be either IdealDelay or FirstOrderPadeDelay\n";

		Delay_Hf_.RunFilter(DelayType);
		
		DelayCircuit_.CalculateOutputSignal();
		Cfd_.RunCFD();
		CFD_ShaperHf_.RunFilter(FilterType::CDFShaper);
		CFD_Shaper_.CalculateOutputSignal();
		ZC_Comp.RunComparator();
		Arming_Comp.RunComparator();
        OutputLogicCirc.RunLogic();

	}

	void FbcmFrontEndChip::GetHitAnalysisInfo(int BXC_SlotNo, HitAnalysisInfo & HitTotToaInfo){
	    Hit_Analyzer.RunHitAnalyzer(BXC_SlotNo);
	    HitTotToaInfo = Hit_Analyzer.GetHitAnalysisInfo();
    }


}


