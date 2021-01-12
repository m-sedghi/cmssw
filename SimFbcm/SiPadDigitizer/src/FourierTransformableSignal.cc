///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimFbcm/SiPadDigitizer/interface/FourierTransformableSignal.h"

namespace FbcmFE {

	FourierTransformableSignal::FourierTransformableSignal(){};

	FourierTransformableSignal::~FourierTransformableSignal(){};

	void FourierTransformableSignal::init(ComplexSignalType & CmplxSrcSig){
		CplxSignal_=CmplxSrcSig;
		CplxSignal_.shrink_to_fit();
		fft_Status_=NOTHING;
		shift_Status_=NA;
		long len=CplxSignal_.size(), m, n=1, i;
		for (m=0; len>0; m++)
			len>>=1;
		m--; /// check and correction is requared.
		for (i=0; i<m; i++)
			n <<= 1;

		nfft=n;
		pwr2=m;

	}

	void FourierTransformableSignal::init(SignalType & RealSrcSig){
		CplxSignal_.clear();
		for (auto val : RealSrcSig)
			CplxSignal_.emplace_back(Complx(val));

		CplxSignal_.shrink_to_fit();


		fft_Status_=NOTHING;
		shift_Status_=NA;
		long len=CplxSignal_.size(), m, n=1, i;
		for (m=0; len>0; m++)
			len>>=1;
		m--; /// check and correction is requared.
		for (i=0; i<m; i++)
			n <<= 1;

		nfft=n;
		pwr2=m;

	}




	void FourierTransformableSignal::fft()
	{
		if (fft_Status_!=FWD){
			MS_FFT(1,  CplxSignal_);
			fft_Status_=FWD;
			shift_Status_=DC_Side;
		}

	}

	void FourierTransformableSignal::ifft()
	{
		if (fft_Status_==FWD){
			MS_FFT(-1,  CplxSignal_);
			fft_Status_=INV;
		}
	}


	void FourierTransformableSignal::MS_FFT(int dir,  std::vector<Complx> & x)
	{

		//    This computes an in-place complex-to-complex FFT
		//    dir =  1 gives forward transform
		//    dir = -1 gives reverse transform
		long n,m,i,i1,j,k,i2,l,l1,l2;
		Complx TempCplx, c,u,t;

		// Calculate the number of points
		//m=PWR_of_Nfft;
		//Nfft=2^PWR_of_Nfft (Matlab Syntax)
		//n=Nfft;

		m=pwr2;
		n=nfft;

		// Do the bit reversal
		i2 = n >> 1;
		j = 0;

		for (i=0;  i<n-1;  i++)		{
			if (i < j) {
				// swap--> Maybe there is an easer way with vectors!
				TempCplx=x[i];
				x[i] = x[j];
				x[j] = TempCplx;
			}

			k = i2;
			while (k <= j) 	{
				j -= k;
				k >>= 1;
			}
			j += k;
		}

		// Compute the FFT

		c.Real(-1.0);
		c.Imag(0.0);



		l2 = 1;
		for (l=0 ; l < m; l++) 	{
			l1 = l2;
			l2 <<= 1;
			u.Real(1.0);
			u.Imag(0.0);

			for (j=0; j<l1 ; j++) {
				for (i=j;i<n;i+=l2) {
					i1 = i + l1;

					t=u*x[i1];
					x[i1]=x[i]-t;
					x[i]=x[i]+t;
				}
				u=u*c;
			}

			c.Imag(sqrt((1.0-c.Real())/2.0));
			if (dir == 1)
				c.Imag(-c.Imag());

			c.Real(sqrt((1.0+c.Real())/2.0));

		}


		// Scaling for Backward transform

		if (dir == -1) 	{
			for (i=0;i<n;i++) 	{
				x[i]=x[i]/((double)n);
			}
		}
	}


	void FourierTransformableSignal::ApplyFilter(Filter & fl) {
		this->fft();
		ComplexSignalType & FilterTF=fl.GetTransferFunction();
		ComplexSignalType TmpCplxSig(nfft);

		for (unsigned int i=0; i< nfft ; i++)
			TmpCplxSig[i]=CplxSignal_[i]*FilterTF[i];

		CplxSignal_.swap(TmpCplxSig);
		// let's remain in fourier domain. 
		//this->ifft(); // this will be done when want to get real signal
		//ClearImagParts(); //
	}

	void FourierTransformableSignal::ClearImagParts() {
		for (ComplexSignalType::iterator it=CplxSignal_.begin(); it < CplxSignal_.end() ; it++)
			it->SetZeroImagPart();
	}

	void FourierTransformableSignal::GetRealSignal (SignalType & SrcRealSig){
		ifft();
		ClearImagParts();
		SrcRealSig.clear();
		for (auto v:CplxSignal_)
			SrcRealSig.emplace_back(v.Real());
		SrcRealSig.shrink_to_fit();
	}


}  



