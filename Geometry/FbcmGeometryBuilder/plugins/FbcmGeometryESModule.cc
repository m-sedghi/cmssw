///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------
#include "Geometry/FbcmGeometryBuilder/src/FbcmGeometryBuilderFromDDD.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"

//#include "Geometry/Records/interface/FbcmRecoGeometryRcd.h"
#include "CondFormats/GeometryObjects/interface/RecoIdealGeometry.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"

// Alignments
#include "CondFormats/Alignment/interface/DetectorGlobalPosition.h"
#include "CondFormats/Alignment/interface/AlignmentErrorsExtended.h"
#include "CondFormats/AlignmentRecord/interface/GlobalPositionRcd.h"

//#include "CondFormats/AlignmentRecord/interface/FbcmAlignmentRcd.h"
//#include "CondFormats/AlignmentRecord/interface/FbcmAlignmentErrorExtendedRcd.h"
#include "Geometry/CommonTopologies/interface/GeometryAligner.h"

#include <memory>

using namespace edm;

class FbcmGeometryESModule : public edm::ESProducer {
public:
  /// Constructor
  FbcmGeometryESModule(const edm::ParameterSet& p);

  /// Destructor
  ~FbcmGeometryESModule() override;

  /// Produce FbcmGeometry.
  std::unique_ptr<FbcmGeometry> produce(const FbcmGeometryRecord& record);

private:
  // use the DDD as Geometry source
  bool useDDD_;
  bool applyAlignment_;
  const std::string alignmentsLabel_;
  
  edm::ESGetToken<DDCompactView, IdealGeometryRecord> cpvToken_;
  //edm::ESGetToken<RecoIdealGeometry, FbcmRecoGeometryRcd> rig_FbcmToken;
  
  //edm::ESGetToken<Alignments, GlobalPositionRcd> globalPositionToken_;
  //edm::ESGetToken<Alignments, FbcmAlignmentRcd> alignmentsToken_;
  //edm::ESGetToken<AlignmentErrorsExtended, FbcmAlignmentErrorExtendedRcd> alignmentErrorsToken_;
};

////-----------------------------------------------------------------------------------------------------------

FbcmGeometryESModule::FbcmGeometryESModule(const edm::ParameterSet& p)
    : useDDD_(p.getParameter<bool>("useDDD")),
      applyAlignment_(p.getParameter<bool>("applyAlignment")),
      alignmentsLabel_(p.getParameter<std::string>("alignmentsLabel")) {
  auto cc = setWhatProduced(this);
  
  //std::cout << "Hello From Constructor\n"; 
  
  if (useDDD_) {
    //cc.setConsumes(cpvToken_); // for CMSSW_11_0_X
	cpvToken_ = cc.consumes(); // for CMSSW_11_2_X
  } else {
    //cc.setConsumes(rig_FbcmToken);
	//rig_FbcmToken = cc.consumes();
	;
  }
  if (applyAlignment_) {
		/* globalPositionToken_ = cc.consumes(edm::ESInputTag{"", alignmentsLabel_});
		alignmentsToken_ = cc.consumes(edm::ESInputTag{"", alignmentsLabel_});
		alignmentErrorsToken_ = cc.consumes(edm::ESInputTag{"", alignmentsLabel_}); */
	edm::LogError("FBCM Geometry Produceer") << "Alignment has not been implemented \n"; 
  }
}

FbcmGeometryESModule::~FbcmGeometryESModule() {}

std::unique_ptr<FbcmGeometry> FbcmGeometryESModule::produce(const FbcmGeometryRecord& record) {
  
  // std::cout << "Hello From FbcmGeometryESModule produce Method\n"; 
  // std::cout << "useDDD_:" << useDDD_ << "\n"; 
  // std::cout << "alignmentsLabel_:" << alignmentsLabel_ << "\n"; 
  
  auto FBCM_Geometry = std::make_unique<FbcmGeometry>();
//rig : RecoIdealGeometry
  if (useDDD_) {
    auto cpv = record.getTransientHandle(cpvToken_);
    FbcmGeometryBuilderFromDDD builder;
    builder.build(*FBCM_Geometry, cpv.product());
	
	// std::vector<FbcmSiPadGeom const*> allSiPadGeoms=FBCM_Geometry->SiPads();
	// std::cout << "FBCM_Geometry has just been built!:\n";
	
		// for(const auto& SiPadGeomDet: allSiPadGeoms) {
			// std::cout << SiPadGeomDet->id();
			 ///Cool! it is OK!
		// }
	
	//return std::unique_ptr<FbcmGeometry>(builder.build(cpv.product()));
  } else {
	edm::LogError("FBCM Geometry Produceer") << "FbcmGeometryBuilderFromCondDB from CondDB has not been implemented \n"; 
    //const auto& rig_Fbcm = record.get(rig_FbcmToken);
    //FbcmGeometryBuilderFromCondDB builder;
    //builder.build(*FBCM_Geometry, rig_Fbcm);
	////return std::unique_ptr<FbcmGeometry>(builder.build(rig_Fbcm));
  }

  if (applyAlignment_) {
    edm::LogError("FBCM Geometry Produceer") << "Alignment has not been implemented \n"; 
    	
  }

  return FBCM_Geometry;

}

//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//using edm::FbcmGeometryESModule;
//DEFINE_FWK_MODULE(FbcmGeometryESModule);
DEFINE_FWK_EVENTSETUP_MODULE(FbcmGeometryESModule);
