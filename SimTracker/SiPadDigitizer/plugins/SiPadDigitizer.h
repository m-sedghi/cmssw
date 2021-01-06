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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ProducesCollector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Framework/interface/ESWatcher.h"
//#include "SimTracker/SiPhase2Digitizer/plugins/DigitizerUtility.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

//#include "Geometry/Records/interface/TrackerTopologyRcd.h"
//#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"



namespace CLHEP {
  class HepRandomEngine;
}


namespace edm {
  class ConsumesCollector;
  class Event;
  class EventSetup;
  class ParameterSet;
  template <typename T>
  class Handle;
  class StreamID;
}  // namespace edm

class MagneticField;
class PileUpEventPrincipal;
class PSimHit;

class SiPadDigitizerAlgorithm;
class FbcmGeometryRecord;

class FbcmSiPadGeom;
class FbcmGeometry;
class FbcmSiPadTopology;



namespace cms {
  class SiPadDigitizer : public DigiAccumulatorMixMod {
  public:
    //typedef std::unordered_map<unsigned, TrackerGeometry::ModuleType> ModuleTypeCache;  
	
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

/*
    void StorePileupInformation(std::vector<int>& numInteractionList,
                                std::vector<int>& bunchCrossingList,
                                std::vector<float>& TrueInteractionList,
                                std::vector<edm::EventID>& eventInfoList,
                                int bunchSpacing) override 
								{
			PileupInfo_ = std::make_unique<PileupMixingContent>( numInteractionList, bunchCrossingList,
																 TrueInteractionList, eventInfoList, bunchSpacing);
    }
		PileupMixingContent* getEventPileupInfo() override { return PileupInfo_.get(); }
*/
  private:
    void accumulateSiPadHits(edm::Handle<std::vector<PSimHit> >, size_t globalSimHitIndex);
	//Local3DPoint CMSUnits(Local3DPoint lp);
	//void ConvertPSimHitsToCMSUnits(std::vector<PSimHit> const& simHits_org, std::vector<PSimHit> &simHitsCMSUnit);
	
    bool firstInitializeEvent_;
    //bool firstFinalizeEvent_;
    std::unique_ptr<SiPadDigitizerAlgorithm> SiPadDigiAlgo;


    std::map<std::string, size_t> crossingSimHitIndexOffset_;

    typedef std::vector<std::string> vstring;
    const std::string hitsProducer_;
	const std::string SubdetName_;
    const vstring trackerContainers_;
    const std::string geometryType_;
	const std::string InstanceName_;
	
    edm::ESHandle<FbcmGeometry> theFbcmGeom;
    edm::ESHandle<MagneticField> pSetup; // OK
    std::map<unsigned int, FbcmSiPadGeom const*> SiPadsIdGeomMap;
	
	//edm::ESHandle<TrackerTopology> tTopoHand;
	//edm::ESWatcher<TrackerDigiGeometryRecord> theDigiGeomWatcher;
	edm::ESWatcher<FbcmGeometryRecord> theGeomWatcher;
	
    //CLHEP::HepRandomEngine* randomEngine_ = nullptr;
    std::unique_ptr<PileupMixingContent> PileupInfo_;

    //const bool pilotBlades;         // Default = false
    //const int NumberOfEndcapDisks;  // Default = 2
	//const bool IsReadoutAnalog_; // No Need, I will alway put analog data
	//const bool makeDigiSimLinks_;
	//ModuleTypeCache moduleTypeCache_;
	const std::vector<edm::ParameterSet>& conf_FE;
	unsigned int DigiEventCnt; 
	edm::EDGetTokenT< std::vector<PSimHit> > m_token;
    // infrastructure to reject dead pixels as defined in db (added by F.Blekman)
	//edm::InputTag simSrc_;
	//const edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_;
  };
}  // namespace cms

#endif
