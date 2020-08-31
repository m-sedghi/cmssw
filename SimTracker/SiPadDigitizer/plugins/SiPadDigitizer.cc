//
// system include files
#include <memory>
#include <set>

// user include files
#include "SiPadDigitizer.h"
#include "SiPadDigitizerAlgorithm.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"


#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"
#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"
//Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Exception.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//using namespace std;
namespace cms {
  SiPadDigitizer::SiPadDigitizer(const edm::ParameterSet& iConfig,
                                     edm::ProducesCollector producesCollector,
                                     edm::ConsumesCollector& iC)
      : firstInitializeEvent_(true),
        firstFinalizeEvent_(true), // No need ?!
		SiPadDigiAlgo(),
        hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")), //g4SimHits
		SubdetName_(iConfig.getParameter<std::string>("SubdetName")), // FBCMHits
        geometryType_(iConfig.getParameter<std::string>("GeometryType")),
		IsReadoutAnalog_(iConfig.getParameter<bool>("isReadoutAnalog")),
		makeDigiSimLinks_(iConfig.getUntrackedParameter<bool>("makeDigiSimLinks")) 	// ?
		//simSrc_(iConfig.getParameter<edm::InputTag>("SimFBCM_Tag")),
		//simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(simSrc_))
		 //_trackerContainers_(iConfig.getParameter<std::vector<std::string> >("RoutList")),
		 //pilotBlades(iConfig.exists("enablePilotBlades") ? iConfig.getParameter<bool>("enablePilotBlades") : false),
        //NumberOfEndcapDisks(iConfig.exists("NumPixelEndcap") ? iConfig.getParameter<int>("NumPixelEndcap") : 2) 
		{
    edm::LogInfo("SiPadDigitizer ") << "Enter the SiPad Digitizer";
	
    const std::string alias("simSiPadDigis");
	
    producesCollector.produces<edm::DetSetVector<PixelDigi> >("SiPad").setBranchAlias(alias);
	producesCollector.produces<edm::DetSetVector<PixelDigi> >("SiPadOverT").setBranchAlias(alias);
    if (makeDigiSimLinks_) {
      producesCollector.produces<edm::DetSetVector<PixelDigiSimLink> >("SiPadSimLink").setBranchAlias(alias);
    }
	
	
	
	edm::InputTag tag(hitsProducer_, SubdetName_);
	iC.consumes<std::vector<PSimHit> >(tag);
	//iC.consumes<edm::DetSetVector<PixelDigiSimLink> >(edm::InputTag("simSiPadDigis", "SiPad"));
	
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "SiPadDigitizer requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }
    SiPadDigiAlgo.reset(new SiPadDigitizerAlgorithm(iConfig));
	//(std::unique_ptr<SiPadDigitizerAlgorithm>) (const edm::ParameterSet&)
	//SiPadDigiAlgo(iConfig);
  }

  SiPadDigitizer::~SiPadDigitizer() { edm::LogInfo("SiPadDigitizer ") << "Destruct the SiPad Digitizer"; }

 void SiPadDigitizer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) 
  {    
  //std::cout << "SiPadDigitizer--LumiBlock " << "\n";
    iSetup.get<IdealMagneticFieldRecord>().get(pSetup);
    iSetup.get<TrackerTopologyRcd>().get(tTopoHand);

    if (theDigiGeomWatcher.check(iSetup)) {
      iSetup.get<TrackerDigiGeometryRecord>().get(geometryType_, pDD_);
      //reset cache
      ModuleTypeCache().swap(moduleTypeCache_);
      detectorUnits.clear();
      for (auto const& det_u : pDD_->detUnits()) {
        unsigned int detId_raw = det_u->geographicalId().rawId();
        DetId detId = DetId(detId_raw);
        if (DetId(detId).det() == DetId::Detector::Tracker) {
          const Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u);
          assert(pixdet);
          detectorUnits.insert(std::make_pair(detId_raw, pixdet));
        }
      }
    }
	
  }


  
  void SiPadDigitizer::accumulateSiPadHits(edm::Handle<std::vector<PSimHit> > hSimHits,
                                                   size_t globalSimHitIndex,
                                                   const unsigned int tofBin)
  {
     if (hSimHits.isValid()) {
      std::set<unsigned int> detIds;
      std::vector<PSimHit> const& simHits = *hSimHits.product();
      int indx = 0;
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex)
	  {
        unsigned int detId_raw = (*it).detUnitId();
        if (detectorUnits.find(detId_raw) == detectorUnits.end())
          continue;
        if (detIds.insert(detId_raw).second) 
		{
          // The insert succeeded, so this detector element has not yet been processed.
          const Phase2TrackerGeomDetUnit* phase2det = detectorUnits[detId_raw];
          
		  // access to magnetic field in global coordinates
          GlobalVector bfield = pSetup->inTesla(phase2det->surface().position());
          LogDebug("PixelDigitizer") << "B-field(T) at " << phase2det->surface().position()
                                     << "(cm): " << pSetup->inTesla(phase2det->surface().position());
          //next line temporarily was disabled for compilation
		  SiPadDigiAlgo->accumulateSimHits(it, itEnd, globalSimHitIndex, tofBin, phase2det, bfield);
		  
        }
        indx++;
      }
    }

	
  }



  void SiPadDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup) {
	  
	  // Cache random number engine
    edm::Service<edm::RandomNumberGenerator> rng;
	if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "SiPad requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }
	
	 if (firstInitializeEvent_) 
	 {
      SiPadDigiAlgo->init(iSetup);
      firstInitializeEvent_ = false;
	 }
	
	
    //randomEngine_ = &rng->getEngine(e.streamID());

    SiPadDigiAlgo->initializeEvent(rng->getEngine(e.streamID()));
	
   

    // Make sure that the first crossing processed starts indexing the sim hits from zero.
    // This variable is used so that the sim hits from all crossing frames have sequential
    // indices used to create the digi-sim link (if configured to do so) rather than starting
    // from zero for each crossing.
    crossingSimHitIndexOffset_.clear();

  }

  void SiPadDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
    // Step A: Get Inputs

	  edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer_, SubdetName_);
	  
	  //edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(tag));
      //iEvent.getByToken(simHitToken_, simHits);
	  
      iEvent.getByLabel(tag, simHits);
	  
      unsigned int tofBin = PixelDigiSimLink::LowTof;

      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
      //       std::cout << "index offset, current hit count = " << crossingSimHitIndexOffset_[tag.encode()] << ", " << simHits->size() << std::endl;
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
    
  }

  void SiPadDigitizer::accumulate(PileUpEventPrincipal const& iEvent,
                                    edm::EventSetup const& iSetup,
                                    edm::StreamID const& streamID) 
	{
       edm::Handle<std::vector<PSimHit> > simHits;
       edm::InputTag tag(hitsProducer_, SubdetName_);
		
		//edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(tag));
      //iEvent.getByToken(simHitToken_, simHits);
      iEvent.getByLabel(tag, simHits);
	  
	  
      unsigned int tofBin = PixelDigiSimLink::LowTof;

      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
      //       std::cout << "index offset, current hit count = " << crossingSimHitIndexOffset_[tag.encode()] << ", " << simHits->size() << std::endl;
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
	}

  void SiPadDigitizer::finalizeEvent(edm::Event& iEvent, const edm::EventSetup& iSetup) 
  {
    ///------------------------
    //addSiPadCollection(iEvent, iSetup); // becuase in this fuction, I will call only one function, 
	// I will ghather the code of addSiPadCollection directly inside this function. 
    //-------- the function body of addSiPadCollection ----------
		
	////-------
    const TrackerTopology* tTopo = tTopoHand.product();
    std::vector<edm::DetSet<PixelDigi> > digiVector;
	std::vector<edm::DetSet<PixelDigi> > digiVector_OverT;
    std::vector<edm::DetSet<PixelDigiSimLink> > digiLinkVector;
	
	
    for (auto const& det_u : pDD_->detUnits()) 
	{
      //DetId detId_raw = DetId(det_u->geographicalId().rawId()); // was needed for AlgoType

      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
	  //next line temporarily was disabled for compilation
	  SiPadDigiAlgo->digitize(dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u), digi_map, tTopo);
      	  
      edm::DetSet<PixelDigi> collector(det_u->geographicalId().rawId());
	  edm::DetSet<PixelDigi> collector_OverT(det_u->geographicalId().rawId());
      edm::DetSet<PixelDigiSimLink> linkcollector(det_u->geographicalId().rawId());
      for (auto const& digi_p : digi_map) 
	  {
        DigitizerUtility::DigiSimInfo info = digi_p.second;
		
		// the following two lines are equivalent to the function: addToCollector
		// but the diffrence is info.sig_tot or info.ot_bit
		// edm::DetSet<Phase2TrackerDigi>& collector ----> info.ot_bit
		// edm::DetSet<PixelDigi>& collector ------------> info.sig_tot
		// the addToCollector function uses as tempale: DigiType for DetSet<Phase2TrackerDigi> or DetSet<PixelDigi>
		
        std::pair<int, int> ip = PixelDigi::channelToPixel(digi_p.first);
        collector.data.emplace_back(ip.first, ip.second, info.sig_tot);
		collector_OverT.data.emplace_back(ip.first, ip.second, info.ot_bit);
		
        for (auto const& sim_p : info.simInfoList)
		{
          linkcollector.data.emplace_back(digi_p.first,
                                          sim_p.second->trackId(),
                                          sim_p.second->hitIndex(),
                                          sim_p.second->tofBin(),
                                          sim_p.second->eventId(),
                                          sim_p.first);
        }
      }
      if (!collector.data.empty())
        digiVector.push_back(std::move(collector));
      if (!collector_OverT.data.empty())
        digiVector_OverT.push_back(std::move(collector_OverT));
      if (!linkcollector.data.empty())
        digiLinkVector.push_back(std::move(linkcollector));
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<PixelDigi> > output(new edm::DetSetVector<PixelDigi>(digiVector));
	std::unique_ptr<edm::DetSetVector<PixelDigi> > outputOverT(new edm::DetSetVector<PixelDigi>(digiVector_OverT));
	std::unique_ptr<edm::DetSetVector<PixelDigiSimLink> > outputlink(new edm::DetSetVector<PixelDigiSimLink>(digiLinkVector));

    // Step D: write output to file
    iEvent.put(std::move(output), "SiPad");
	iEvent.put(std::move(outputOverT), "SiPadOverT");
    if (makeDigiSimLinks_) {
      iEvent.put(std::move(outputlink), "SiPadSimLink");
    }
	 
  }


}  // namespace cms

// Check if it poosible to do in other file
//using cms::SiPadDigitizer;
//DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);

#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"

using cms::SiPadDigitizer;
DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);