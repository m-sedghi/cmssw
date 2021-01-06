///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_SiHitPulseShape_h
#define _FbcmFE_SiHitPulseShape_h


#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"
#include "SimTracker/SiPadDigitizer/interface/FftPreparation.h"

using namespace std;
namespace FbcmFE {

	class SiHitPulseShape {
	public:
		
		SiHitPulseShape(const std::vector<double> & param);
		SiHitPulseShape(double _tau, double _r, double _theta, double _nTerms, double _maxCharge);
		~SiHitPulseShape();
		void GetPulseShape(SignalType & timeVect, double _xOffset, SignalType & Pulse);
		void GetPulseSeriesShape(FftPreparation & FftPrepRef, std::vector<TofChargePair> Tof_Charge);
	private:
		double nFactorial(int n);
		double aScalingConstant(int N, int i);
		double SiPulseExpansion(double x);
		double SiHitPulseValue(double x);

		double xOffset;
		double tau;
		double r;
		double theta;
		double nTerms;
		double maxCharge;
	};



}
#endif
