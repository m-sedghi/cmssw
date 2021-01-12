///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <memory>
#include <set>

#include "SiPadDigitizer.h"

namespace cms {
  SiPadDigitizer::SiPadDigitizer(const edm::ParameterSet& iConfig,
                                       edm::ProducesCollector producesCollector,
                                       edm::ConsumesCollector& iC)
      : firstInitializeEvent_(true),
        SiPadDigiAlgo(),
        hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")), //g4SimHits
		SubdetName_(iConfig.getParameter<std::string>("SubdetName")), // FBCMHits
        geometryType_(iConfig.getParameter<std::string>("GeometryType")),
		InstanceName_(iConfig.getParameter<std::string>("InstanceName")),
		conf_FE(iConfig.getParameter<std::vector<edm::ParameterSet>>("SiPadFrontEndParam")) {
    
	edm::LogInfo("SiPadDigitizer ") << "Entered the SiPadDigitizer";
	const std::string alias("simSiPadDigis");
	
    producesCollector.produces<edm::DetSetVector<SiPadDigiData> >(InstanceName_).setBranchAlias(alias);
	
	edm::InputTag tag(hitsProducer_, SubdetName_);
	iC.consumes<std::vector<PSimHit> >(tag);
	m_token= iC.consumes<std::vector<PSimHit> >(tag);
	
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "SiPadDigitizer requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }
    SiPadDigiAlgo.reset(new SiPadDigitizerAlgorithm(iConfig));
  }

  SiPadDigitizer::~SiPadDigitizer() { edm::LogInfo("SiPadDigitizer ") << "Destruct the SiPad Digitizer"; }

 void SiPadDigitizer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) 
  {    
    //std::cout << "SiPadDigitizer-LumiBlock " << "\n"; //call 1
	iSetup.get<IdealMagneticFieldRecord>().get(pSetup);
    		
	// for this version, assuming Ideal-geometry for Digi. Alignment has not been implemented. 
	if (theGeomWatcher.check(iSetup)) {
      //iSetup.get<FbcmGeometryRecord>().get(geometryType_, theFbcmGeom);
	  iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom); // for the moment, FbcmGeometryRecord has no label for idealGeomentry
	  SiPadsIdGeomMap.clear();
	  
      for (auto const& SiPadUnit : theFbcmGeom->SiPads()) {
	    unsigned int detId_raw = SiPadUnit->geographicalId().rawId(); 
		FbcmDetId DetID_Fbcm(detId_raw); // FbcmDetId class itself checks if it is a valid ID for FBMC or not!
		//Another way is as the following:
		//FbcmDetId DetID_Fbcm = SiPadUnit->id();		
		
		  const FbcmSiPadGeom* SiPadDetptr = SiPadUnit;
          assert(SiPadDetptr);
          SiPadsIdGeomMap.insert(std::make_pair(DetID_Fbcm.rawId(), SiPadDetptr));
      }
    }
  }

void SiPadDigitizer::accumulateSiPadHits(edm::Handle<std::vector<PSimHit> > hSimHits,
                                         size_t globalSimHitIndex)
  {
	 //std::cout << "SiPadDigitizer-accumulateSiPadHits " << "\n"; // call 4
     if (hSimHits.isValid()) {
	
	std::vector<PSimHit> const& simHits = *hSimHits.product();
	  std::set<unsigned int> detIds;
	  
      int indx = 0;
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex)
	  {
        unsigned int detId_raw = (*it).detUnitId();
		FbcmDetId fbdetId(detId_raw);
			
        if (SiPadsIdGeomMap.find(detId_raw) == SiPadsIdGeomMap.end())
          continue;
        if (detIds.insert(detId_raw).second) // if it is a new detId_raw
		{ // The insert succeeded, so this detector element has not yet been processed.
          const FbcmSiPadGeom* SiPadSensorGeom = SiPadsIdGeomMap[detId_raw];
          // access to magnetic field in global coordinates
		  GlobalVector bfield = pSetup->inTesla(SiPadSensorGeom->surface().position());
		  
          LogDebug("PixelDigitizer") << "B-field(T) at " << SiPadSensorGeom->surface().position()
                                     << "(cm): " << pSetup->inTesla(SiPadSensorGeom->surface().position());
          
		  	enum { LowTof, HighTof };
			unsigned int tofBin = LowTof; 
		 SiPadDigiAlgo->accumulateSimHits(it, itEnd, globalSimHitIndex, tofBin, SiPadSensorGeom, bfield);
        }
        indx++;
      }
    }
 }

  void SiPadDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup)   {
	  //std::cout << "SiPadDigitizer-initializeEvent " << "\n"; // call 2
	  // Cache random number engine
    edm::Service<edm::RandomNumberGenerator> rng;
	if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "SiPad requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }
	
	 if (firstInitializeEvent_) {
      SiPadDigiAlgo->init(iSetup);
      firstInitializeEvent_ = false;
	 }

    SiPadDigiAlgo->initializeEvent(rng->getEngine(e.streamID()));

    // Make sure that the first crossing processed starts indexing the sim hits from zero.
    // This variable is used so that the sim hits from all crossing frames have sequential
    // indices used to create the digi-sim link (if configured to do so) rather than starting
    // from zero for each crossing.
    crossingSimHitIndexOffset_.clear();
  }

  void SiPadDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  //std::cout << "SiPadDigitizer-accumulate-Event " << "\n"; // call 3-1

  // Step A: Get Inputs

	  edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer_, SubdetName_);
	  
	  //edm::EDGetTokenT< std::vector<PSimHit> > simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(tag));
      //iEvent.getByToken(simHitToken_, simHits);
      iEvent.getByLabel(tag, simHits); 
	// from other refrecens: 36 ns ~ 3*12.06 ns where, 12.06 ns is sigam of gausian spread of electrons
	//https://indico.cern.ch/event/7522/contributions/1248050/subcontributions/112540/attachments/1048712/1494938/IEEE2006CMSTkSimulation.pdf
	// ToF < 36 ns ==> LowToF
	// ToF > 36 ns ==> HighToF
	// for Fbcm, we should consider LowToF, becuase Tof~9.4 ns for the zero Bx
      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()]);
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
	  
      accumulateSiPadHits(simHits, crossingSimHitIndexOffset_[tag.encode()]);
     
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
	}

  void SiPadDigitizer::finalizeEvent(edm::Event& iEvent, const edm::EventSetup& iSetup) 
  {
    std::vector<edm::DetSet<SiPadDigiData> > DigiDataVector;
			
    for (auto const& SiPadUnit : theFbcmGeom->SiPads()) {
      std::map<int, SiPadDigiData> SiPadDigilMap; // first: channelNo(pixel) per SiPad --> by default each SiPad comprises one pixel
												  // Second: DigiDataResult per SiPad (channel)
	  
	  SiPadDigiAlgo->GetDigiResults(SiPadUnit, SiPadDigilMap); // this fills out the SiPadDigilMap
	  
      //FbcmDetId SiPdetId(SiPadUnit->geographicalId().rawId());
      edm::DetSet<SiPadDigiData> collector(SiPadUnit->geographicalId().rawId()); 
      for (auto const& Ch_DigiRes : SiPadDigilMap) // for each pixel(channel) in SiPad
	  {
        SiPadDigiData DigiResult = Ch_DigiRes.second;
		collector.data.emplace_back(
					  DigiResult.SideIndex(),
					  DigiResult.StationIndex(),
					  DigiResult.SiliconDieIndex(),
					  DigiResult.SiPadIndex(),
					  DigiResult.Radius(),
					  DigiResult.Phi_Degrees(),
					  DigiResult.Area(),
					  DigiResult.ChargeSum()  ,
					  DigiResult.CahrgePsimVector()  ,
					  DigiResult.BxSlotHitAnalysisVector() 
					);
      }
      if (!collector.data.empty()) {
        DigiDataVector.push_back(std::move(collector));
	  }
		
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<SiPadDigiData> > output(new edm::DetSetVector<SiPadDigiData>(DigiDataVector));

    // Step D: write output to file
    iEvent.put(std::move(output), InstanceName_);
	
  }
  
}  // namespace cms


#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"
using cms::SiPadDigitizer;
DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);
