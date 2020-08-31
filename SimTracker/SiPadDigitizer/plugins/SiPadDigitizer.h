#ifndef SiPadDigitizer_h
#define SiPadDigitizer_h

/** \class SiPadDigitizer

 ************************************************************/

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
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimTracker/SiPhase2Digitizer/plugins/DigitizerUtility.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

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
class TrackerDigiGeometryRecord;

class PixelGeomDetUnit;
class TrackerGeometry;

class PixelTopology;
using SiPadTopology = PixelTopology; // Typedef SiPad

namespace cms {
  class SiPadDigitizer : public DigiAccumulatorMixMod {
  public:
    typedef std::unordered_map<unsigned, TrackerGeometry::ModuleType> ModuleTypeCache;  
	
    explicit SiPadDigitizer(const edm::ParameterSet& iConfig, 
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
      PileupInfo_ = std::make_unique<PileupMixingContent>(
          numInteractionList, bunchCrossingList, TrueInteractionList, eventInfoList, bunchSpacing);
    }
*/

//    PileupMixingContent* getEventPileupInfo() override { return PileupInfo_.get(); }

  private:
    void accumulateSiPadHits(edm::Handle<std::vector<PSimHit> >,
                             size_t globalSimHitIndex,
                             const unsigned int tofBin);

    bool firstInitializeEvent_;
    bool firstFinalizeEvent_;
    std::unique_ptr<SiPadDigitizerAlgorithm> SiPadDigiAlgo;
	
    /** @brief Offset to add to the index of each sim hit to account for which crossing it's in.
*
* I need to know what each sim hit index will be when the hits from all crossing frames are merged into
* one collection (assuming the MixingModule is configured to create the crossing frame for all sim hits).
* To do this I'll record how many hits were in each crossing, and then add that on to the index for a given
* hit in a given crossing. This assumes that the crossings are processed in the same order here as they are
* put into the crossing frame, which I'm pretty sure is true.<br/>
* The key is the name of the sim hit collection. */

    std::map<std::string, size_t> crossingSimHitIndexOffset_;

    typedef std::vector<std::string> vstring;
    const std::string hitsProducer_;
	const std::string SubdetName_;
    const vstring trackerContainers_;
    const std::string geometryType_;
    edm::ESHandle<TrackerGeometry> pDD_;
    edm::ESHandle<MagneticField> pSetup;
    std::map<unsigned int, PixelGeomDetUnit const*> detectorUnits;
	edm::ESHandle<TrackerTopology> tTopoHand;
	edm::ESWatcher<TrackerDigiGeometryRecord> theDigiGeomWatcher;
    //CLHEP::HepRandomEngine* randomEngine_ = nullptr;
    std::unique_ptr<PileupMixingContent> PileupInfo_;

    //const bool pilotBlades;         // Default = false
    //const int NumberOfEndcapDisks;  // Default = 2
	const bool IsReadoutAnalog_; // No Need, I will alway put analog data
	const bool makeDigiSimLinks_;
	ModuleTypeCache moduleTypeCache_;
    // infrastructure to reject dead pixels as defined in db (added by F.Blekman)
	//edm::InputTag simSrc_;
	//const edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_;
  };
}  // namespace cms

#endif
