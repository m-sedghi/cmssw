///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/NonlinearLimiter.h"

namespace FbcmFE {
    NonlinearLimiter::NonlinearLimiter(SignalType & InputSignal, SignalType & OutputSignal):
		InputSignal_(InputSignal),
		OutputSignal_(OutputSignal)
	{ 
	
    };

	NonlinearLimiter::~NonlinearLimiter() {};
		
	void NonlinearLimiter::RunLimiter(){

		unsigned int len=InputSignal_.size();
        SignalType ScaledSig(len);
	    //double MaxLimAbs=700; // mV

	    double LimmiterCorr=1.5; // fix for S1 SmoothStep function
	    double edge0=-LimmiterCorr*MaxLimAbs;
	    double edge1=-edge0;
	    ScaledSig= (InputSignal_ - edge0) / (edge1-edge0);
        Clamp(ScaledSig, 0.0, 1.0);

        for (unsigned int i=0; i < len ; i++ )
            OutputSignal_[i]= MaxLimAbs*2*((ScaledSig[i] * ScaledSig[i] * (3.0- 2.0*ScaledSig[i]))-0.5);

	}
	void NonlinearLimiter::Clamp(SignalType & Sig_, double lowerlimit, double upperlimit){

	    for (unsigned int i=0; i < Sig_.size() ; i++)
	    {
	        if (Sig_[i] < lowerlimit)
                Sig_[i]=lowerlimit;
            else if (Sig_[i] > upperlimit)
                Sig_[i]=upperlimit;

	    }


	}

}


