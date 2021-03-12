///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "DataFormats/FbcmDigi/interface/HitAnalysisInfo.h"

namespace FbcmFE {
    HitAnalysisInfo::HitAnalysisInfo(int BXC_SlotNo,
									//std::vector<ToaTotPair> & TotToaVect,
                                    std::vector<ToaTotPair>  TotToaVect,
                                     float AlignerDelay,
                                    HitStatus BxHitStatus,
                                    unsigned int NumOfRecognizedHits
									//,
                                    //float _ToAUpperCut_,
                                    //float _ToALowerCut_,
                                    //float _BinLen_,
                                    //float _BinShift_
									):
        BXC_SlotNo_(BXC_SlotNo),
        TotToaVect_(TotToaVect),
        AlignerDelay_(AlignerDelay),
        BxHitStatus_((char)BxHitStatus),
        NumOfRecognizedHitsInBx_(NumOfRecognizedHits)
		//,
        //ToAUpperCut_(_ToAUpperCut_),
        //ToALowerCut_(_ToALowerCut_),
        //BinLen_(_BinLen_),
        //BinShift_(_BinShift_) 
		{};


	HitAnalysisInfo::HitAnalysisInfo():
	    BXC_SlotNo_(0),
        TotToaVect_(),
        AlignerDelay_(0.0),
        //BxHitStatus_(HitStatus::Zero),
		BxHitStatus_(HitStatus_Zero),
        NumOfRecognizedHitsInBx_(0)
		//,
        //ToAUpperCut_(0.0),
        //ToALowerCut_(0.0),
        //BinLen_(0.0),
        //BinShift_(0.0)
		{};
	
	void HitAnalysisInfo::clear(){
		BXC_SlotNo_ = 0;
		TotToaVect_.clear();
		AlignerDelay_ = 0.0 ;
		//BxHitStatus_ = HitStatus::Zero;
		BxHitStatus_ = HitStatus_Zero;
		NumOfRecognizedHitsInBx_ = 0;
        //ToAUpperCut_ = 0.0;
        //ToALowerCut_ = 0.0;
        //BinLen_ = 0.0;
        //BinShift_ = 0.0;
	}

std::ostream& operator<<(std::ostream& o, HitAnalysisInfo& Info) {
   o << "BxSlotNo: " << Info.BxSlotNo() << ", "
            << "AlignerDelay: " << Info.lpGBTAlignerDelay() << ", "
            << "TotToaVect_size: " << Info.TotToaVectort().size() << ", "
            << "HitStatus: " << Info.Bx_HitStatus() << ", "
            << "nbrOfRecognizedHits: " << Info.nbrOfRecognizedHitsInBx() << "\n" ;
            //<< "ToAUpperCut: " << Info.ToAUpperCut() << ", "
            //<< "ToALowerCut: " << Info.ToALowerCut() << ", "
            //<< "BinLen: " << Info.BinLen() << ", "
            //<< "BinShift: " << Info.BinShift() << "\n";
	for (auto v:Info.TotToaVectort()){
	o << v;
	}
			
		return o;
}


}


