///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_HitAnalysisInfo_h
#define _FbcmFE_HitAnalysisInfo_h

#include "DataFormats/FbcmDigi/interface/ToaTotPair.h"
#include <ostream>

namespace FbcmFE {

#define HitStatus_Uncertain  -1
#define HitStatus_Zero 		  0
#define HitStatus_NonZero     1

	enum  HitStatus : char
	{
		Uncertain      = HitStatus_Uncertain,
		Zero           = HitStatus_Zero,
		NonZero        = HitStatus_NonZero
	};


	class HitAnalysisInfo {
	    public:
        
		explicit HitAnalysisInfo(int BXC_SlotNo,
						//std::vector<ToaTotPair> & TotToaVect,
                        std::vector<ToaTotPair> TotToaVect,
                        float AlignerDelay,
                        HitStatus BxHitStatus,
                        unsigned int NumOfRecognizedHits
						//,
                        //float _ToAUpperCut_,
                        //float _ToALowerCut_,
                        //float _BinLen_,
                        //float _BinShift_
						);
		HitAnalysisInfo();
		~HitAnalysisInfo(){};
		void clear();
        int BxSlotNo() const {return BXC_SlotNo_;}
		std::vector<ToaTotPair> TotToaVectort() const {return TotToaVect_;}
		float lpGBTAlignerDelay() const {return AlignerDelay_;}
		HitStatus Bx_HitStatus() const {return (HitStatus)BxHitStatus_;}
		int Bx_HitStatusInt() const {return (int)BxHitStatus_;}
		unsigned int nbrOfRecognizedHitsInBx() const { return NumOfRecognizedHitsInBx_; }
		//float ToAUpperCut() const {return ToAUpperCut_;}
		//float ToALowerCut() const {return ToALowerCut_;}
		//float BinLen() const {return BinLen_;}
		//float BinShift() const {return BinShift_;}

		inline bool operator<(const HitAnalysisInfo& other) const { return BxSlotNo() < other.BxSlotNo(); }

	private:
		int BXC_SlotNo_;
		std::vector<ToaTotPair> TotToaVect_;
		float AlignerDelay_;
		//HitStatus BxHitStatus_;
		char BxHitStatus_;
		unsigned int NumOfRecognizedHitsInBx_;
        //float ToAUpperCut_;
        //float ToALowerCut_;
        //float BinLen_;
        //float BinShift_;

	};

	std::ostream& operator<<(std::ostream& o, HitAnalysisInfo& Info);

}
#endif
