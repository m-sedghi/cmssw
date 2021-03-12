///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_FourierTransformableSignal_h
#define _FbcmFE_FourierTransformableSignal_h


#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "SimFbcm/SiPadDigitizer/interface/Complx.h"
#include "SimFbcm/SiPadDigitizer/interface/Filter.h"



namespace FbcmFE {


	enum FFTStatus
	{
		NOTHING,
		FWD,
		INV
	};

	enum ShiftStatus
	{
		NA,
		DC_centered,
		DC_Side
	};

	class FourierTransformableSignal {
	public:
		FourierTransformableSignal();
		~FourierTransformableSignal();
		void ApplyFilter(Filter &);
		void init(ComplexSignalType & SrcCplxSig);
		void init(SignalType & SrcRealSig);
		void GetRealSignal (SignalType & SrcRealSig);
	private:
		void fft();
		void ifft();
		void ClearImagParts();
		void MS_FFT(int dir,  ComplexSignalType &);

		ComplexSignalType CplxSignal_;
		FFTStatus fft_Status_;
		ShiftStatus shift_Status_;
		unsigned int pwr2;
		unsigned int nfft;
	};



}
#endif
