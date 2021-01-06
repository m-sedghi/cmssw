///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_HitAnalyzer_h
#define _FbcmFE_HitAnalyzer_h

#include "SimTracker/SiPadDigitizer/interface/GeneralUtilities.h"
#include "SimTracker/SiPadDigitizer/interface/FftPreparation.h"
//#include "SimTracker/SiPadDigitizer/interface/HitAnalysisInfo.h"
#include "DataFormats/FbcmDigi/interface/HitAnalysisInfo.h"

namespace FbcmFE {

	class HitAnalyzer {

    public:
		HitAnalyzer(FftPreparation & FFtPrep,LogicSignalType & CFD_LogicalSignal);
		~HitAnalyzer();
		void RunHitAnalyzer(int BXC_SlotNo);
		SignalType & GetAlignedTimeVect(){return timeVectAligned_;}

		HitAnalysisInfo GetHitAnalysisInfo();
		void SetParameters( const edm::ParameterSet * FEParamPtr );

	private:
	    unsigned int CheckHitInBx();
	    ToAStatus GetToAStatus(float ToA, float ToT, float HalfBxLen);

	    FftPreparation & FFtPrep_;

        SignalType timeVectAligned_;
		LogicSignalType & CFD_LogicalSignal_;
		std::vector<ToaTotPair> TotToaVect;
		int BXC_SlotNo_;
		double BX_Duration_;
		double AlignerDelay_;
		unsigned int FS_;
		HitStatus BxHitStatus;
		unsigned int NumOfRecognizedHits ;
        double ToAUpperCut_ ;  // for BIB study, more than one BX should be investigated
        double ToALowerCut_;   // for BIB study, more than one BX should be investigated
        double BinLen_ ;
        double BinShift_;

	};

}
#endif
