///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#ifndef _FbcmFE_Complx_h
#define _FbcmFE_Complx_h

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"

namespace FbcmFE {

	class Complx {
	public:
		Complx(double r, double i)   { real_ = r; imag_ = i; }
		Complx(double r) { real_ = r; imag_ = 0.0; }
		Complx() { real_ = imag_ = 0.0; }
		~Complx();
		double Real()  {return real_;}
		double Imag() {return imag_;}
		double mag() const ;
		double magSq() const ;
		double angle() const;
		Complx Reciprocal();
		void Real(double r) {real_=r;}
		void Imag(double i) {imag_=i;}
		void Set(double Re, double Im) { Real(Re); Imag(Im); }
		void Set(Complx cplx) { Real(cplx.Real()); Imag(cplx.Imag()); }
		Complx conjugate()  { Imag(-Imag()); return *this;}
		double abs()  { return mag();}
		double abs2()  {  return magSq();}
		void SetZeroRealPart() {real_=0.0;}
		void SetZeroImagPart() {imag_=0.0;}


		Complx operator+( const double &) const;//addition
		Complx operator-( const double &) const;//subtraction
		Complx operator+( Complx &) const;//addition
		Complx operator-( Complx &) const;//subtraction
		Complx operator*( Complx &) const;//multiplication
		Complx operator*(const double &) const;//multiplication
		Complx operator/( Complx &) const;
		Complx operator/(const double &) const;
		const Complx &operator=( Complx );//assignment
		bool const   operator==( Complx &) const;//equivalent
		bool const   operator!=( Complx &) const;//not equivalent

	private:
		double real_;
		double imag_;

	};

	typedef std::vector<Complx> ComplexSignalType;
}
#endif
