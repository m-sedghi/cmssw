///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#include "SimFbcm/SiPadDigitizer/interface/HitAnalyzer.h"

namespace FbcmFE {

    HitAnalyzer::HitAnalyzer(FftPreparation & FFtPrep,LogicSignalType & CFD_LogicalSignal):
    FFtPrep_(FFtPrep),
    CFD_LogicalSignal_(CFD_LogicalSignal),
    BX_Duration_(25.0),
    FS_(FFtPrep_.SamplingRepetition()),
    BxHitStatus(HitStatus::Zero),
    NumOfRecognizedHits(0)	{
		timeVectAligned_=FFtPrep_.TimeVectZeroCenter()-AlignerDelay_;
		};

	HitAnalyzer::~HitAnalyzer(){};

    void HitAnalyzer::RunHitAnalyzer(int BXC_SlotNo){

        TotToaVect.clear();
		const float Eps=0.0001;
		
        float HalfBxLen=BX_Duration_/2.0;
        int Len=timeVectAligned_.size();
        BXC_SlotNo_=BXC_SlotNo;
        unsigned int AlignerDelayIndexShift=(unsigned int)(AlignerDelay_ * FS_);
        unsigned int ZeroTimeIndex=Len/2 + AlignerDelayIndexShift;
        unsigned int BX_IndShift=(unsigned int)(BX_Duration_ * FS_);
        int BX_CenterIndex=ZeroTimeIndex + BX_IndShift * BXC_SlotNo_ ;
		
        int lowerTimeIndexCut=BX_CenterIndex + (int)(ToALowerCut_*FS_) ; // BX_CenterIndex-BX_IndShift/2;
        int UpperTimeIndexCut=BX_CenterIndex + (int)(ToAUpperCut_*FS_) ; // BX_CenterIndex+BX_IndShift/2;

        bool HitDetected=false;
        unsigned int start_ = 0, end_;
        int BckwardCNT,i;
        float ToA,ToT;
        ToaTotPair TotToa;
        ToAStatus ToAState=ToAStatus::UnKnown;

        float binshift=BinShift_-Eps;
        int16_t SubBxBinNo;
        bool LastState=CFD_LogicalSignal_ [lowerTimeIndexCut-1];
        for (i=lowerTimeIndexCut; ((i <= UpperTimeIndexCut) || HitDetected==true ) && (i < Len) ; i++)
        {

            if ( i==lowerTimeIndexCut && CFD_LogicalSignal_[lowerTimeIndexCut]) {
                for (BckwardCNT=lowerTimeIndexCut-1; (CFD_LogicalSignal_ [BckwardCNT]) && BckwardCNT >= 0  ; BckwardCNT--) ;
                BckwardCNT++;

                start_=BckwardCNT;
                HitDetected=true;
            }
            else if (CFD_LogicalSignal_[i] != LastState)  // rising-edge or falling-edge occurred
            {
                if (CFD_LogicalSignal_[i]) /// rising-edge
                {
                    start_=i;
                    HitDetected=true;
                }
                else if (HitDetected) /// falling-edge & ValidPulse
                {
                    end_=i;
                    HitDetected=false;
                    ToA=timeVectAligned_[start_]-(BXC_SlotNo_*BX_Duration_);
                    ToT=timeVectAligned_[end_]-timeVectAligned_[start_];
                    SubBxBinNo= (int16_t)(round((ToA-binshift)/BinLen_));

                    ToAState=GetToAStatus(ToA,ToT,HalfBxLen);

                    TotToa.SetPairInfo(ToA,ToT, true, ToAState ,SubBxBinNo );
                    TotToaVect.emplace_back(TotToa);

//                    std::cout << "A hit detected, ToA: " << TotToa.ToA() << ", "
//                            << "ToT: " << TotToa.ToT() << "\n\n";

                }

            }
            LastState=CFD_LogicalSignal_[i];

        }

        if (HitDetected) // if HitDetected is still true, this means that a too long hit observed reaching to the end of time span. this is not a valid hit.
        {

                    HitDetected = false;
                    ToA=timeVectAligned_[start_]-(BXC_SlotNo_*BX_Duration_);
                    SubBxBinNo= (int16_t)(round((ToA-binshift)/BinLen_));
                    //end_=i-1;
                    //ToT=timeVectAligned_[end_]-timeVectAligned_[start_];
                    ToT=1000.0;

                    ToAState=GetToAStatus(ToA,ToT,HalfBxLen);

                    TotToa.SetPairInfo(ToA,ToT,false, ToAState , SubBxBinNo );
                    TotToaVect.emplace_back(TotToa);
        }
            NumOfRecognizedHits=CheckHitInBx();
    }

    unsigned int HitAnalyzer::CheckHitInBx(){
        NumOfRecognizedHits=0;
        //float HalfBx=-BX_Duration_/2.0;
         BxHitStatus=HitStatus::Zero;
           for (ToaTotPair v:TotToaVect)
             {
                 if (v.ToAToAStatus()==ToAStatus::FullyWithinBx || v.ToAToAStatus()==ToAStatus::WithinBx_LastsAfter ){
                        NumOfRecognizedHits++;
                        if (BxHitStatus==HitStatus::Zero)
                            BxHitStatus=HitStatus::NonZero;
                     }
                else if (v.ToAToAStatus()==ToAStatus::BeforeBx_Overlapping )
                        BxHitStatus=HitStatus::Uncertain;
             }

        return NumOfRecognizedHits;

    }



    HitAnalysisInfo HitAnalyzer::GetHitAnalysisInfo(){
    HitAnalysisInfo tmp(BXC_SlotNo_,
                        TotToaVect,
                        AlignerDelay_,
                        BxHitStatus,
                        NumOfRecognizedHits );
    return tmp;
    }


    ToAStatus HitAnalyzer::GetToAStatus(float ToA, float ToT, float HalfBxLen){
        ToAStatus Status=ToAStatus::UnKnown;
        float Ending=ToA+ToT;
        if (ToA > HalfBxLen)
            Status=ToAStatus::AfterBx;
        else if ( (ToA >= -HalfBxLen) && (Ending > HalfBxLen) )
            Status=ToAStatus::WithinBx_LastsAfter;
        else if ( (ToA >= -HalfBxLen) && (Ending <= HalfBxLen) )
            Status=ToAStatus::FullyWithinBx;
        else if ( (ToA < -HalfBxLen) && (Ending > -HalfBxLen) )
            Status=ToAStatus::BeforeBx_Overlapping;
        else //if ( (ToA < -HalfBxLen) && (Ending < -HalfBxLen) )
            Status=ToAStatus::BeforeBx_NonOverlapping;

        return Status;

    }
	
	
	void HitAnalyzer::SetParameters( const edm::ParameterSet * FEParamPtr ){
		
		
		BX_Duration_ = FEParamPtr->getParameter< double >("Bx_Duration");
		AlignerDelay_ = FEParamPtr->getParameter< double >("lpGBT_AlignerDelay");
		ToAUpperCut_  = FEParamPtr->getParameter< double >("ToAUpperCut");
        ToALowerCut_  = FEParamPtr->getParameter< double >("ToALowerCut");
        BinLen_  = FEParamPtr->getParameter< double >("BinLength");
		BinShift_ = FEParamPtr->getParameter< double >("BinOffset");
		
		timeVectAligned_=FFtPrep_.TimeVectZeroCenter()-AlignerDelay_;
		
		// std::cout << "BX_Duration_: " << BX_Duration_ << "\n"
					// << "AlignerDelay_: " << AlignerDelay_ << "\n"
					// << "ToAUpperCut_: " << ToAUpperCut_ << "\n"
					// << "ToALowerCut_: " << ToALowerCut_ << "\n"
					// << "BinLen_: " << BinLen_ << "\n"
					// << "BinShift_: " << BinShift_ << "\n" ;
					
	}

}


