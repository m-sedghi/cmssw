///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/FftPreparation.h"

namespace FbcmFE {
    FftPreparation::FftPreparation(const edm::ParameterSet& FFT_SimParam) {
             FS = FFT_SimParam.getParameter<int>("SamplingRepetition");;
             FixNfft_as_pwr2(FFT_SimParam.getParameter<int>("NumOfFFT_Points"));
             CreateTimeVect_ZeroFirst();
             CreateTimeVect_ZeroCenter();
             CreateFreqVect();
             SignalZeroInit();

         };

	FftPreparation::~FftPreparation() {};


		void FftPreparation::CreateTimeVect_ZeroFirst(){
		timeVect_Zerofirst.clear();
		double T=1.0/((double)FS),t;
		for (unsigned int i=0; i < Nfft ; i++) {
			t=(double)i*T;
			timeVect_Zerofirst.emplace_back(t);
		}
		timeVect_Zerofirst.shrink_to_fit();
	}

	void FftPreparation::CreateTimeVect_ZeroCenter(){
		timeVect_ZeroCenter.clear();
		double T=1.0/((double)FS),t;
		for (unsigned int i=0; i < Nfft ; i++) {
			t=(double)i*T-((double)(Nfft/2)/FS);
			timeVect_ZeroCenter.emplace_back(t);
		}
		timeVect_ZeroCenter.shrink_to_fit();
	}


	void FftPreparation::CreateFreqVect(){
		freqVect_GHz.clear();
		double FN=((double)FS)/Nfft;
		double ftmp;
		for (unsigned int i=0; i < Nfft ; i++) {
			ftmp=FN*i-((Nfft)/2)*FN;
			freqVect_GHz.emplace_back(ftmp);
		}
		freqVect_GHz.shrink_to_fit();
		// freqVect_GHz is in the general frequency order.
		// we need to perform fftshift to be compatible with fft components
		//FbcmFE::fftshift(freqVect_GHz);
		fftshift(freqVect_GHz);
	}

    void FftPreparation::SignalZeroInit(){
    SignalType temp(Nfft,0.0);
    TD_Signal=temp;
    }

    void FftPreparation::FixNfft_as_pwr2(int nbrFFT_points){
        unsigned int  m, n=1, i;

		for (m=0; nbrFFT_points>0; m++)
			nbrFFT_points>>=1;
		m--;

		for (i=0; i<m; i++)
			n <<= 1;
        Nfft=n;
        pwr2=m;

    }


    template <class VectType>
	void FftPreparation::fftshift(VectType & src){
		VectType tmp(src.size());
		unsigned int midIndx= src.size()/2;

		for (unsigned int i=0; i < midIndx ; i++)
			tmp[i+midIndx]=src[i];

		for (unsigned int i=midIndx; i < src.size() ; i++)
			tmp[i-midIndx]=src[i];

		src.swap(tmp);
	};


}


