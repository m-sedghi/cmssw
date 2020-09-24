#ifndef Geometry_Records_FbcmGeometryRecord_H
#define Geometry_Records_FbcmGeometryRecord_H

///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "boost/mpl/vector.hpp"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/GeometryFileRcd.h"
#include "CondFormats/AlignmentRecord/interface/GlobalPositionRcd.h"


class FbcmGeometryRecord 
: public edm::eventsetup::DependentRecordImplementation<
				FbcmGeometryRecord,
				boost::mpl::vector<	IdealGeometryRecord,
									GeometryFileRcd, 
									GlobalPositionRcd> > {};

#endif
