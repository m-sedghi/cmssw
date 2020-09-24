///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/Records/interface/FbcmGeometryRecord.h"

class FbcmGeometryTest : public edm::one::EDAnalyzer<> {
public:
  explicit FbcmGeometryTest(const edm::ParameterSet&) {}

  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}
};

void FbcmGeometryTest::analyze(const edm::Event&, const edm::EventSetup& iEventSetup) {
  edm::LogVerbatim("Geometry") << "FbcmGeometryTest::analyze";
  std::cout << "Hello Mohammad from analyzer before calling ESProducer\n"; 
    
  edm::ESTransientHandle<FbcmGeometry> FbcmGeom;
  iEventSetup.get<FbcmGeometryRecord>().get(FbcmGeom);
  
   edm::ESHandle<FbcmGeometry> FbcmGeom2;
  iEventSetup.get<FbcmGeometryRecord>().get(FbcmGeom2);
  
  std::cout << "Check the delivery of FbcmGeom from ESProducer:\n"; 
  
  
  std::vector<FbcmSiPadGeom const*> allSiPadGeoms=FbcmGeom.product()->SiPads();
  
  std::vector<FbcmSiPadGeom const*> allSiPadGeoms2=FbcmGeom2->SiPads();
  
  //FbcmGeom.product();
  
  	for(const auto& SiPadGeomDet: allSiPadGeoms2) {
			std::cout << SiPadGeomDet->id();
			 ////Cool! it is OK!
		}
  
 
}

DEFINE_FWK_MODULE(FbcmGeometryTest);