///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

    SignalType operator+(const SignalType& lhs, const SignalType& rhs) {
        if(lhs.size() != rhs.size())// Vectors must be the same size
            throw "Can't add two vectors of different sizes!";
        SignalType result(lhs.size());
        for(unsigned int i=0; i < lhs.size(); i++)
            result[i]= lhs[i] + rhs[i];
        return result;
	}

	SignalType operator-(const SignalType& lhs, const SignalType& rhs) {
        if(lhs.size() != rhs.size())// Vectors must be the same size
            throw "Can't add two vectors of different sizes!";
        SignalType result(lhs.size());
        for(unsigned int i=0; i < lhs.size(); i++)
            result[i]= lhs[i] - rhs[i];
        return result;
	}

	SignalType operator*(const double& lhs, const SignalType& rhs){
	 SignalType result(rhs.size());
        for(unsigned int i=0; i < rhs.size(); i++)
            result[i]= lhs * rhs[i];
        return result;
	}

	SignalType operator*(const SignalType& lhs, const double& rhs){
     SignalType result(lhs.size());
    for(unsigned int i=0; i < lhs.size(); i++)
        result[i]= lhs[i] * rhs ;
    return result;
	}


	SignalType operator/(const SignalType& lhs, const double& rhs){
     SignalType result(lhs.size());
    for(unsigned int i=0; i < lhs.size(); i++)
        result[i]= lhs[i] / rhs ;
    return result;
	}

    SignalType operator+(const double& lhs, const SignalType& rhs){
	 SignalType result(rhs.size());
        for(unsigned int i=0; i < rhs.size(); i++)
            result[i]= lhs + rhs[i];
        return result;
	}

	SignalType operator+(const SignalType& lhs, const double& rhs){
     SignalType result(lhs.size());
    for(unsigned int i=0; i < lhs.size(); i++)
        result[i]= lhs[i] + rhs ;
    return result;
	}

	    SignalType operator-(const double& lhs, const SignalType& rhs){
	 SignalType result(rhs.size());
        for(unsigned int i=0; i < rhs.size(); i++)
            result[i]= lhs - rhs[i];
        return result;
	}

	SignalType operator-(const SignalType& lhs, const double& rhs){
     SignalType result(lhs.size());
    for(unsigned int i=0; i < lhs.size(); i++)
        result[i]= lhs[i] - rhs ;
    return result;
	}

}
