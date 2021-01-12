#ifndef Fbcm_SiPadDigitizer_h
#define Fbcm_SiPadDigitizer_h
///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixMod.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ProducesCollector.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiDataCollection.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "SimFbcm/SiPadDigitizer/plugins/SiPadDigitizerAlgorithm.h"
//#include "DataFormats/Math/interface/angle_units.h"
//using angle_units::operators::convertRadToDeg;


namespace cms {
  class SiPadDigitizer : public DigiAccumulatorMixMod {
  public:
    explicit SiPadDigitizer(  const edm::ParameterSet& iConfig, 
									edm::ProducesCollector,
									edm::ConsumesCollector& iC);

    ~SiPadDigitizer() override;

    void initializeEvent(edm::Event const& e, edm::EventSetup const& c) override;
    void accumulate(edm::Event const& e, edm::EventSetup const& c) override;
    void accumulate(PileUpEventPrincipal const& e, edm::EventSetup const& c, edm::StreamID const&) override;
    void finalizeEvent(edm::Event& e, edm::EventSetup const& c) override;
	void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) override;
    virtual void beginJob() {}


    void StorePileupInformation(std::vector<int>& numInteractionList,
                                std::vector<int>& bunchCrossingList,
                                std::vector<float>& TrueInteractionList,
                                std::vector<edm::EventID>& eventInfoList,
                                int bunchSpacing) override {
			PileupInfo_ = std::make_unique<PileupMixingContent>( numInteractionList, bunchCrossingList,
																 TrueInteractionList, eventInfoList, bunchSpacing);
    }
		PileupMixingContent* getEventPileupInfo() override { return PileupInfo_.get(); }

  private:
    void accumulateSiPadHits(edm::Handle<std::vector<PSimHit> >, size_t globalSimHitIndex);
	
    bool firstInitializeEvent_;
    std::unique_ptr<SiPadDigitizerAlgorithm> SiPadDigiAlgo;
	std::map<std::string, size_t> crossingSimHitIndexOffset_;
    const std::string hitsProducer_;
	const std::string SubdetName_;
    const std::string geometryType_;
	const std::string InstanceName_;
    edm::ESHandle<FbcmGeometry> theFbcmGeom;
    edm::ESHandle<MagneticField> pSetup; 
    std::map<unsigned int, FbcmSiPadGeom const*> SiPadsIdGeomMap;
	edm::ESWatcher<FbcmGeometryRecord> theGeomWatcher;
    std::unique_ptr<PileupMixingContent> PileupInfo_;
	const std::vector<edm::ParameterSet>& conf_FE;
	edm::EDGetTokenT< std::vector<PSimHit> > m_token;
   
  };
}  // namespace cms

#endif
