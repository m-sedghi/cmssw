///--------------------------------------------------------------------
//  A set of classes for simulating the hit response of FBCM Front-end
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: November 2020
///--------------------------------------------------------------------

#ifndef _FbcmFE_FeConfigSelector_h
#define _FbcmFE_FeConfigSelector_h

#include "SimFbcm/SiPadDigitizer/interface/GeneralUtilities.h"


namespace FbcmFE {

class FeConfigSelector {
	
	public:
		FeConfigSelector(const std::vector< edm::ParameterSet > & SiPadFrontEndParamVect_);
		~FeConfigSelector();
		std::pair<float, const edm::ParameterSet * > SelectFrontEndConfig(float SensorsSize);
		const edm::ParameterSet * GetSelectedFrontEndConfic() {return ActiveFrontEndParamPtr; }

	private:
		//const std::vector< edm::ParameterSet > & SiPadFrontEndParamVect; // even though the Reference would also work in a standalone compilation,
																		   // but the reference does not work in finalizing event function in cmssw. 
																		   // So, it is essential to keep an exact copy of the ParemeterVect. 
		const std::vector< edm::ParameterSet > SiPadFrontEndParamVect; // keeping an exact copy the parameter Set, not a Reference.
		std::vector < std::pair< int ,  std::vector<double>   > > table;
		const edm::ParameterSet * ActiveFrontEndParamPtr;

	};

}
#endif
