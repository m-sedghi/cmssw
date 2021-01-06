///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------


#include "SimTracker/SiPadDigitizer/interface/FeConfigSelector.h"

namespace FbcmFE {
    FeConfigSelector::FeConfigSelector(const std::vector< edm::ParameterSet > & SiPadFrontEndParamVect_):
	SiPadFrontEndParamVect(SiPadFrontEndParamVect_)
	{ 
		std::vector < double > Range;
		for (unsigned int j=0 ; j < SiPadFrontEndParamVect.size() ; j++ ){
			Range = SiPadFrontEndParamVect[j].getParameter< std::vector<double> >("GoodForSizeRange");
			table.emplace_back( std::make_pair(j,Range) );
		}
    };

	FeConfigSelector::~FeConfigSelector() {};

	std::pair<float, const edm::ParameterSet * > FeConfigSelector::SelectFrontEndConfig(float SensorsSize){
		int SelectedIndex=0;
		for (auto row:table){
			if (( SensorsSize >= row.second[0] ) && ( SensorsSize < row.second[1] ) ){
				SelectedIndex = row.first;
				break;
			}
		}
		
		ActiveFrontEndParamPtr = SiPadFrontEndParamVect.data() + SelectedIndex;
		//ActiveFrontEndParamPtr = & SiPadFrontEndParamVect[SelectedIndex];
		
		return std::make_pair(SensorsSize, ActiveFrontEndParamPtr);
	}


}


