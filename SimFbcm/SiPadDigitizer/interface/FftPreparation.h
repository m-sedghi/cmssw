///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_FftPreparation_h
#define _FbcmFE_FftPreparation_h

#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"


namespace FbcmFE {

class FftPreparation {
	public:
		FftPreparation(const edm::ParameterSet& FFT_SimParam);
		~FftPreparation();


        unsigned int nbrFFTpoints() {return Nfft;}
        unsigned int SamplingRepetition() {return FS;}
		unsigned int pwr2OfNfft(){return pwr2;}
        SignalType & TimeVectZerofirst() {return timeVect_Zerofirst;}
        SignalType & TimeVectZeroCenter() {return timeVect_ZeroCenter;}
		SignalType & FreqVectGHz() {return freqVect_GHz;}
		SignalType & TimeDomainSignal() {return TD_Signal;}
		SignalType & SetSignal(SignalType & sig) {TD_Signal=sig; return TD_Signal;}


	private:
	    void FixNfft_as_pwr2(int nbrFFT_points);
	    void CreateTimeVect_ZeroFirst();
	    void CreateTimeVect_ZeroCenter();
        void CreateFreqVect();
        void SignalZeroInit();

        template <class VectType>
        void fftshift(VectType & src);



        unsigned int Nfft; // Length of signal
        unsigned int FS; //FS: Sampling repetition per ns [1/ns]
		unsigned int pwr2; // Nfft=2^pwr2
		///T = 1/Fs;     // Sampling period [ns]
		///t = (0:Nfft-1)*T; =(0:Nfft-1)*(1/FS);        // Time vector [ns] (Zero first)
		///f = Fs*(0:(L/2))/L; =(0:(Nfft/2))*FS/Nfft = (0:(Nfft/2))*(FS/Nfft) [GHz]
		/// FN=(FS/Nfft) [GHz]
        SignalType timeVect_Zerofirst;
        SignalType timeVect_ZeroCenter;
		SignalType freqVect_GHz;
		SignalType TD_Signal;

	};

}
#endif
