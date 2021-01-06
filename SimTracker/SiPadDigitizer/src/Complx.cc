///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/Complx.h"

namespace FbcmFE {
	///-----------------------------------
	Complx::~Complx(){;}

	double Complx::magSq() const {
		return real_*real_+imag_*imag_;
	}

	double Complx::mag() const {
		return sqrt(magSq());
	}

	double Complx::angle() const{
		return atan2(imag_,real_);
	}

	Complx Complx::operator+(Complx& cplx2) const {
		return Complx(real_ + cplx2.Real(),imag_+ cplx2.Imag());
	}

	Complx Complx::operator-(Complx& cplx2) const {
		return Complx(real_ - cplx2.Real(),imag_ - cplx2.Imag());
	}

	Complx Complx::operator+(const double& dlm) const {
		return Complx(real_ + dlm, imag_);
	}

	Complx Complx::operator-(const double& dlm) const {
		return Complx(real_ - dlm, imag_);
	}

	Complx Complx::operator*(Complx& cplx2) const {
		double tmpRe=(real_ * cplx2.Real()) - (imag_ * cplx2.Imag());
		double tmpIm=(real_ * cplx2.Imag()) + (imag_ * cplx2.Real());
		return Complx(tmpRe,tmpIm);
	}

	Complx Complx::operator*(const double& dlm) const {
		double tmpRe=(real_ * dlm);
		double tmpIm=(imag_ * dlm);
		return Complx(tmpRe,tmpIm);
	}

	Complx Complx::operator/(Complx& cplx2) const {
		double tmpRe=(real_ * cplx2.Real() + imag_ * cplx2.Imag())/cplx2.magSq();
		double tmpIm=(real_ * cplx2.Imag() - imag_ * cplx2.Real())/cplx2.magSq();
		return Complx(tmpRe,tmpIm);
	}

	Complx Complx::operator/(const double& dlm) const {
		double tmpRe= real_ / dlm ;
		double tmpIm= imag_ / dlm ;
		return Complx(tmpRe,tmpIm);
	}


	Complx Complx::Reciprocal()  {
		double tmpRe=real_/magSq();
		double tmpIm=-imag_/magSq();
		real_=tmpRe;
		imag_=tmpIm;
		return *this;
	}

	const Complx& Complx::operator=(Complx right) {
		real_ = right.Real();
		imag_ = right.Imag();
		return *this;
	}


	bool const  Complx::operator==(Complx &right) const {
		if((real_ == right.Real()) && (imag_ == right.Imag()))
			return true;
		else
			return false;
	}

	bool const  Complx::operator!=( Complx &right) const {
		if((real_ != right.Real()) || (imag_ != right.Imag()))
			return true;
		else
			return false;
	}


}  

