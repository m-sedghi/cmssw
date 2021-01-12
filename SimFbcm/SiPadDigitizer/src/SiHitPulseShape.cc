///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimFbcm/SiPadDigitizer/interface/SiHitPulseShape.h"

namespace FbcmFE {

	SiHitPulseShape::SiHitPulseShape(double _tau, double _r, double _theta, double _nTerms, double _maxCharge){
		xOffset=0.0;
		tau=_tau;
		r=_r;
		theta=_theta;
		nTerms=_nTerms;
		maxCharge=_maxCharge;
	};

	SiHitPulseShape::SiHitPulseShape(const std::vector<double> & param){
		xOffset=0.0;
		tau= param[0];
		r= param[1];
		theta= param[2];
		nTerms= param[3];
		maxCharge= param[4];
	};

	double SiHitPulseShape::nFactorial(int n) {
		return std::tgamma(n + 1);
	}
	double SiHitPulseShape::aScalingConstant(int N, int i) {
		return std::pow(-1, (double)i) * nFactorial(N) * nFactorial(N + 2) /  (nFactorial(N - i) * nFactorial(N + 2 - i) * nFactorial(i));
	}

	SiHitPulseShape::~SiHitPulseShape(){};
	void SiHitPulseShape::GetPulseShape(SignalType & timeVect, double _xOffset, SignalType & Pulse){

	    xOffset=_xOffset;
		Pulse.clear();
		for (unsigned int i=0 ; i < timeVect.size(); i++)
			Pulse.emplace_back(SiHitPulseValue(timeVect[i]));

		Pulse.shrink_to_fit();

	}

	void SiHitPulseShape::GetPulseSeriesShape(FftPreparation & FftPrepRef, std::vector<TofChargePair> Tof_Charge){

	    SignalType TempPulse(FftPrepRef.nbrFFTpoints(), 0.0);
        SignalType & TimeVectRef=FftPrepRef.TimeVectZeroCenter();

        for (TofChargePair tof_q:Tof_Charge){
                xOffset=tof_q.first;
                for (unsigned int i=0 ; i < FftPrepRef.nbrFFTpoints(); i++)
                    TempPulse[i]+=(SiHitPulseValue(TimeVectRef[i])*tof_q.second);
        }

        FftPrepRef.SetSignal(TempPulse);
	}

	double SiHitPulseShape::SiPulseExpansion(double x) {

		double fN = 0;
		double xx = x - xOffset;
		if (xx < 0)
			return 0;

		for (int i = 0; i < nTerms; i++) {
			double angularTerm = 0;
			double temporalTerm = 0;
			double rTerm = std::pow(r, i) / (std::pow(tau, 2. * i) * nFactorial(i + 2));
			for (int j = 0; j <= i; j++) {
				angularTerm += std::pow(std::cos(theta), (double)(i - j)) * std::pow(std::sin(theta), (double)j);
				temporalTerm += aScalingConstant(i, j) * std::pow(xx, (double)(i - j)) * std::pow(tau, (double)j);
			}
			double fi = rTerm * angularTerm * temporalTerm;

			fN += fi;
		}
		return fN;
	}

	double SiHitPulseShape::SiHitPulseValue(double x) {
		double xx = x - xOffset;
		return maxCharge * (std::exp(-xx / tau) * std::pow(xx / tau, 2.) * SiPulseExpansion(x));
	}

}


