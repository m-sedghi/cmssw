///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_GeneralUtilities_h
#define _FbcmFE_GeneralUtilities_h

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace FbcmFE {

    const double e_q=1.602176634e-19;
    const double _ns=1.0e-9;
    const double _uA=1.0e-6;

	typedef std::vector<double> SignalType;
	typedef std::vector<bool> LogicSignalType;
	typedef std::pair<float, float> TofChargePair;

	SignalType operator-(const SignalType& lhs, const SignalType& rhs);
    SignalType operator-(const double& lhs, const SignalType& rhs);
	SignalType operator-(const SignalType& lhs, const double& rhs);

	SignalType operator+(const SignalType& lhs, const SignalType& rhs);
    SignalType operator+(const double& lhs, const SignalType& rhs);
	SignalType operator+(const SignalType& lhs, const double& rhs);

	SignalType operator*(const double& lhs, const SignalType& rhs);
	SignalType operator*(const SignalType& lhs, const double& rhs);


    SignalType operator/(const SignalType& lhs, const double& rhs);

}

#endif
