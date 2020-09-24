///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

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
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"


#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"

#include "DataFormats/FbcmDigi/interface/SiPadAmplitude.h"
#include "DataFormats/FbcmDigi/interface/SiPadAmplitudeCollection.h"

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

#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//



//using SiPadDigi = PixelDigi;
//using SiPadDigiSimLink = PixelDigiSimLink;
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
	//std::cout << "SiPadDigitizer-Constructor " << "\n"; // call 0
	
    const std::string alias("simSiPadAmpl");
	
    producesCollector.produces<edm::DetSetVector<SiPadAmplitude> >("SiPad").setBranchAlias(alias);
	producesCollector.produces<edm::DetSetVector<PixelDigi> >("SiPad2").setBranchAlias(alias);
	//producesCollector.produces<edm::DetSetVector<SiPadAmplitude> >("SiPadOverT").setBranchAlias(alias);
    if (makeDigiSimLinks_) {
      //producesCollector.produces<edm::DetSetVector<SiPadDigiSimLink> >("SiPadSimLink").setBranchAlias(alias);
    }
	
	
	
	edm::InputTag tag(hitsProducer_, SubdetName_);
	iC.consumes<std::vector<PSimHit> >(tag);
	//iC.consumes<edm::DetSetVector<SiPadDigiSimLink> >(edm::InputTag("simSiPadDigis", "SiPad"));
	
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
    //std::cout << "SiPadDigitizer-LumiBlock " << "\n"; //call 1
	iSetup.get<IdealMagneticFieldRecord>().get(pSetup);
    ///edm::ESTransientHandle<FbcmGeometry> FbcmGeom;
		
	// for this version, Ideal fro Digi and Alignment have not been implemented. 
	if (theGeomWatcher.check(iSetup)) {
      //iSetup.get<FbcmGeometryRecord>().get(geometryType_, theFbcmGeom);
	  iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom);
	  SiPadsIdGeomMap.clear();
//	  std::cout <<"Indide the theGeomWatcher\n"; 
	  
      for (auto const& SiPadUnit : theFbcmGeom->SiPads()) {
	    unsigned int detId_raw = SiPadUnit->geographicalId().rawId(); 
		FbcmDetId DetID_Fbcm(detId_raw); // FbcmDetId class itself checks if it is a valid ID for FBMC or not!
		//The other way is as the following:
		//FbcmDetId DetID_Fbcm = SiPadUnit->id();		
		
		//const FbcmSiPadGeom* SiPadDetptr = dynamic_cast<const FbcmSiPadGeom*>(SiPadUnit);
		  const FbcmSiPadGeom* SiPadDetptr = SiPadUnit;
          assert(SiPadDetptr);
          SiPadsIdGeomMap.insert(std::make_pair(DetID_Fbcm.rawId(), SiPadDetptr));
        
      }
	  
	  // for (auto const& it : SiPadsIdGeomMap)
		// {
			// std::cout << it.first << ":  "; 
			// std::cout << it.second->id();
		// }
	  
	  // throw cms::Exception("Manual Stop im LuminosityBlock: ") << " After the for-loop";
	  
    }
	
  }


  
  void SiPadDigitizer::accumulateSiPadHits(edm::Handle<std::vector<PSimHit> > hSimHits,
                                                   size_t globalSimHitIndex)
  {
	 //std::cout << "SiPadDigitizer-accumulateSiPadHits " << "\n"; // call 4
     if (hSimHits.isValid()) {
		 
		 
		 
      std::set<unsigned int> detIds;
      std::vector<PSimHit> const& simHits = *hSimHits.product();
	  
	  //std::cout << "A hSimHits found in accumulateSiPadHits. Size= " << simHits.size()  << "\n";
	  
      int indx = 0;
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex)
	  {
        unsigned int detId_raw = (*it).detUnitId();
			
			
			//FbcmDetId fbdetId(detId_raw);
			std::cout << "--New SimHits found: " << detId_raw << ", i.e.: ";
		
        if (SiPadsIdGeomMap.find(detId_raw) == SiPadsIdGeomMap.end())
          continue;
        if (detIds.insert(detId_raw).second) // if it is a new detId_raw
		{
			FbcmDetId dc(detId_raw);
			std::cout << dc; 
          // The insert succeeded, so this detector element has not yet been processed.
          const FbcmSiPadGeom* SiPadSensorGeom = SiPadsIdGeomMap[detId_raw];
          
		  std::cout <<"Read From retured FbcmSiPadGeom: " << SiPadSensorGeom->id(); 
		  std::cout << SiPadSensorGeom->surface().position(); 
		  std::cout << "B-Field(T)" << pSetup->inTesla(SiPadSensorGeom->surface().position()) << "\n"; 
		  
		  std::cout << "tof:" << it->timeOfFlight() << ", "
					<< "pabs:" << it->pabs() << ", "
					<< "energyLoss:" << it->energyLoss() << ", "
					<< "trackId:" << it->trackId() << ", "
					<< "exitPoint:" << it->exitPoint() << ", "
					<< "localPosition:" << it->localPosition() << "\n";
					
		  std::cout << (*it) << "\n";
		  
		  
		  // access to magnetic field in global coordinates
		  // next line temporarily disabled
          GlobalVector bfield = pSetup->inTesla(SiPadSensorGeom->surface().position());
		  
          LogDebug("PixelDigitizer") << "B-field(T) at " << SiPadSensorGeom->surface().position()
                                     << "(cm): " << pSetup->inTesla(SiPadSensorGeom->surface().position());
          
		  	enum { LowTof, HighTof };
			unsigned int tofBin = HighTof;
		  //next line temporarily was disabled for compilation
		 SiPadDigiAlgo->accumulateSimHits(it, itEnd, globalSimHitIndex, tofBin, SiPadSensorGeom, bfield);
		  
        }
        indx++;
      }
    }

	
  }



  void SiPadDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup) 
  {
	  //std::cout << "SiPadDigitizer-initializeEvent " << "\n"; // call 2
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

  void SiPadDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
  {
  //std::cout << "SiPadDigitizer-accumulate-Event " << "\n"; // call 3-1

  // Step A: Get Inputs

	  edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer_, SubdetName_);
	  
	  //edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(tag));
      //iEvent.getByToken(simHitToken_, simHits);
	  
      iEvent.getByLabel(tag, simHits);
	  
	  //enum { LowTof, HighTof };
	  
      //unsigned int tofBin = PixelDigiSimLink::LowTof;
	  //tofBin= PixelDigiSimLink:HighTof; 
	  //std::cout << "LowTof=" << tofBin << ", " << "HighTof=" << PixelDigiSimLink::HighTof << "\n"; 
	  

	// from other refrecens: 36 ns ~ 3*12.06 ns where, 12.06 ns is sigam of gausian spread of electrons
	//https://indico.cern.ch/event/7522/contributions/1248050/subcontributions/112540/attachments/1048712/1494938/IEEE2006CMSTkSimulation.pdf
	// ToF < 36 ns ==> LowToF
	// ToF > 36 ns ==> HighToF
	// for Fbcm, we should consider HighToF, becuase Tof~59.8 ns
      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()]);
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
		//std::cout << "SiPadDigitizer-accumulate-PileUpEvent " << "\n"; // call 3-2
		
       edm::Handle<std::vector<PSimHit> > simHits;
       edm::InputTag tag(hitsProducer_, SubdetName_);
		
		//edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(tag));
      //iEvent.getByToken(simHitToken_, simHits);
      iEvent.getByLabel(tag, simHits);
	  
	  
      //unsigned int tofBin = PixelDigiSimLink::LowTof;

      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()]);
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
    //std::cout << "SiPadDigitizer-finalizeEvent " << "\n"; // call 5, at the end of the envent

    //const TrackerTopology* tTopo = tTopoHand.product();
	const TrackerTopology* tTopo = NULL; // I don' need this 
    std::vector<edm::DetSet<SiPadAmplitude> > AmplVector;
	
	//PixelDigi
	
	std::vector<edm::DetSet<PixelDigi> > digiVector_OverT;
    //std::vector<edm::DetSet<SiPadDigiSimLink> > digiLinkVector;
	
	
    for (auto const& SiPadUnit : theFbcmGeom->SiPads()) 
	{
      //DetId detId_raw = DetId(SiPadUnit->geographicalId().rawId()); // was needed for AlgoType
		//std::cout << "Mohammad1304 first Loop\n";
		
      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
	  //next line temporarily was disabled for compilation
	  SiPadDigiAlgo->digitize(SiPadUnit, digi_map, tTopo); // this fill out the digi_map
	  //SiPadDigiAlgo->GetAmplitude(SiPadUnit, digi_map); // this fill out the digi_map
      	  FbcmDetId SiPdetId(SiPadUnit->geographicalId().rawId());
		  
		  //std::cout << SiPdetId; 
		  //std::cout << SiPadUnit->geographicalId().rawId() << "\n"; 
		  
      edm::DetSet<SiPadAmplitude> collector(SiPadUnit->geographicalId().rawId());
	  edm::DetSet<PixelDigi> collector_OverT(SiPadUnit->geographicalId().rawId());
      //edm::DetSet<SiPadDigiSimLink> linkcollector(SiPadUnit->geographicalId().rawId());
      for (auto const& digi_p : digi_map) 
	  {
        DigitizerUtility::DigiSimInfo info = digi_p.second;
		
		// the following two lines are equivalent to the function: addToCollector
		// but the diffrence is info.sig_tot or info.ot_bit
		// edm::DetSet<Phase2TrackerDigi>& collector ----> info.ot_bit
		// edm::DetSet<SiPadAmplitude>& collector ------------> info.sig_tot
		// the addToCollector function uses as tempale: DigiType for DetSet<Phase2TrackerDigi> or DetSet<SiPadAmplitude>
		std::cout << "Hello Mohammad\n";
		std::cout << SiPdetId; 
		// Not completly implemented !!!
		collector.data.emplace_back(
		SiPdetId.Side(),
		SiPdetId.Station(),
		SiPdetId.SiliconDie(),
		SiPdetId.SiPad(),
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		info.sig_tot,
		0.0,
		0.0,
		0.0,
		SiPdetId.rawId(),
		0);
		
		
        std::pair<int, int> ip = PixelDigi::channelToPixel(digi_p.first);
        //collector.data.emplace_back(ip.first, ip.second, info.sig_tot);
		collector_OverT.data.emplace_back(ip.first, ip.second, info.ot_bit);
		
        // for (auto const& sim_p : info.simInfoList)
		// {
          // linkcollector.data.emplace_back(digi_p.first,
                                          // sim_p.second->trackId(),
                                          // sim_p.second->hitIndex(),
                                          // sim_p.second->tofBin(),
                                          // sim_p.second->eventId(),
                                          // sim_p.first);
        // }
		
		
      }
      if (!collector.data.empty())
	  {
        AmplVector.push_back(std::move(collector));
		std::cout << "AmplVector is not empty\n"; 
	  }
		
      if (!collector_OverT.data.empty())
	  {
		  digiVector_OverT.push_back(std::move(collector_OverT));
		  std::cout << "collector_OverT is not empty\n"; 
	  }
      // if (!linkcollector.data.empty())
        // digiLinkVector.push_back(std::move(linkcollector));
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<SiPadAmplitude> > output(new edm::DetSetVector<SiPadAmplitude>(AmplVector));
	std::unique_ptr<edm::DetSetVector<PixelDigi> > outputOverT(new edm::DetSetVector<PixelDigi>(digiVector_OverT));
	//std::unique_ptr<edm::DetSetVector<SiPadDigiSimLink> > outputlink(new edm::DetSetVector<SiPadDigiSimLink>(digiLinkVector));

    // Step D: write output to file
    iEvent.put(std::move(output), "SiPad");
	iEvent.put(std::move(outputOverT), "SiPad2");
    if (makeDigiSimLinks_) {
      //iEvent.put(std::move(outputlink), "SiPadSimLink");
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