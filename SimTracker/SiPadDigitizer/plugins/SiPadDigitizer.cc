///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <memory>
#include <set>

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

//#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
//#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
//#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
//#include "Geometry/CommonTopologies/interface/PixelTopology.h"
//#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
//#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"



#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiDataCollection.h"



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"

//#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"

//Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "DataFormats/Math/interface/angle_units.h"
using angle_units::operators::convertRadToDeg;


namespace cms {
  SiPadDigitizer::SiPadDigitizer(const edm::ParameterSet& iConfig,
                                     edm::ProducesCollector producesCollector,
                                     edm::ConsumesCollector& iC)
      : firstInitializeEvent_(true),
        //firstFinalizeEvent_(true), // No need ?!
		SiPadDigiAlgo(),
        hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")), //g4SimHits
		SubdetName_(iConfig.getParameter<std::string>("SubdetName")), // FBCMHits
        geometryType_(iConfig.getParameter<std::string>("GeometryType")),
		InstanceName_(iConfig.getParameter<std::string>("InstanceName")),
		conf_FE(iConfig.getParameter<std::vector<edm::ParameterSet>>("SiPadFrontEndParam")),
		DigiEventCnt(0)
		//IsReadoutAnalog_(iConfig.getParameter<bool>("isReadoutAnalog")),
		//makeDigiSimLinks_(iConfig.getUntrackedParameter<bool>("makeDigiSimLinks")) 	// ?
		//simSrc_(iConfig.getParameter<edm::InputTag>("SimFBCM_Tag")),
		//simHitToken_(edm::ConsumesCollector::consumes< std::vector<PSimHit> >(simSrc_))
		 //_trackerContainers_(iConfig.getParameter<std::vector<std::string> >("RoutList")),
		 //pilotBlades(iConfig.exists("enablePilotBlades") ? iConfig.getParameter<bool>("enablePilotBlades") : false),
        //NumberOfEndcapDisks(iConfig.exists("NumPixelEndcap") ? iConfig.getParameter<int>("NumPixelEndcap") : 2) 
		{
    edm::LogInfo("SiPadDigitizer ") << "Enter the SiPad Digitizer";
	//std::cout << "SiPadDigitizer-Constructor " << "\n"; // call 0
	//std::cout << "InstanceName: " << InstanceName_ << std::endl; 
	//double ddddd=conf_FE[0].getParameter<double>("MaxFEOutputVoltage");
	//std::cout << "\n\n\n" <<  ddddd << "\n\n\n" ;
    const std::string alias("simSiPadDigis");
	
    producesCollector.produces<edm::DetSetVector<SiPadDigiData> >(InstanceName_).setBranchAlias(alias);
	//producesCollector.produces<edm::DetSetVector<PixelDigi> >("SiPad2").setBranchAlias(alias);
	
	//producesCollector.produces<edm::DetSetVector<SiPadDigiData> >("SiPadOverT").setBranchAlias(alias);
    //if (makeDigiSimLinks_) {
      //producesCollector.produces<edm::DetSetVector<SiPadDigiSimLink> >("SiPadSimLink").setBranchAlias(alias);
    //}
	
	
	
	edm::InputTag tag(hitsProducer_, SubdetName_);
	iC.consumes<std::vector<PSimHit> >(tag);
	//iC.consumes<edm::DetSetVector<SiPadDigiSimLink> >(edm::InputTag("simSiPadDigis", "SiPad"));
	m_token= iC.consumes<std::vector<PSimHit> >(tag);
	
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
	

	/*
	std::vector<PSimHit> const& simHits_org= *hSimHits.product();	
	std::vector<PSimHit> simHitsCMSUnit ;	 
	
	
	/// Important note regarding ConvertPSimHitsToCMSUnits: 
	/// This function is needed only when the “TimingSD” class is used as the parent of FbcmSD for the SensitiveDetector.
	/// Because “TimingSD” stores the PSimhits with the unit of “mm” similar to G4; however, “cm” is used for CMS. 
	/// The “TkAccumulatingSensitiveDetector” or “MuonSensitiveDetector”, converts and saves the PsimHits with “cm” units for the distances. 
	/// Thus, if sometime another SensetiveDetector class should be used, then you should double check the dimensions of lengths.
	/// Maybe the following function would not be required. 
	ConvertPSimHitsToCMSUnits(simHits_org, simHitsCMSUnit);
	std::vector<PSimHit> const& simHits = simHitsCMSUnit;
    */
	std::vector<PSimHit> const& simHits = *hSimHits.product();
	  
	  
	  std::set<unsigned int> detIds;
	  
	  //std::cout << "A hSimHits found in accumulateSiPadHits. Size= " << simHits.size()  << "\n";
	  
      int indx = 0;
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex)
	  {
        unsigned int detId_raw = (*it).detUnitId();
			
			
			FbcmDetId fbdetId(detId_raw);
			//std::cout << "--New SimHits found: " << detId_raw << ", i.e.: ";
		
        if (SiPadsIdGeomMap.find(detId_raw) == SiPadsIdGeomMap.end())
          continue;
        if (detIds.insert(detId_raw).second) // if it is a new detId_raw
		{
			//FbcmDetId dc(detId_raw);
			//std::cout << dc; 
          // The insert succeeded, so this detector element has not yet been processed.
          const FbcmSiPadGeom* SiPadSensorGeom = SiPadsIdGeomMap[detId_raw];
          
		  
		  //std::cout <<"Entry and Exitp Points cm " << (*it).entryPoint() << ", " << (*it).exitPoint() << "\n"; 
		  
		  //std::cout <<"NewSipad " << SiPadSensorGeom->id(); 
		  //std::cout <<"x(): " << SiPadSensorGeom->surface().position().x() << std::endl; 
		  //std::cout <<"y(): " << SiPadSensorGeom->surface().position().y() << std::endl; 
		  //Local3DPoint lp1(-1.,-1.,0.);
		  //float x_=SiPadSensorGeom->surface().position().x();
		  //float y_=SiPadSensorGeom->surface().position().y();
		  //float RRadius=TMath::Sqrt(x_*x_+y_*y_); 
		  //float RRadius=SiPadSensorGeom->surface().position().perp(); //OK
		  //float pppi=SiPadSensorGeom->surface().position().barePhi();  // OK // convertRadToDeg(barePhi())
		  //std::cout <<"R, BarePhi: " << lp1.perp() << ", "  << lp1.barePhi() << ", " << convertRadToDeg(lp1.barePhi()) << "\n"; 
		  
		  //std:: cout << "SiPad Center:" << SiPadSensorGeom->surface().toGlobal( LocalPoint(0.,0.,0.)) << std::endl;
		  //std:: cout << "hit poiint:" << SiPadSensorGeom->surface().toGlobal( (*it).localPosition() ) << std::endl;
		  
		  //std::cout << "B-Field(T)" << pSetup->inTesla(SiPadSensorGeom->surface().position()) << "\n"; 
		  
		   //std::cout << "tof:" << it->timeOfFlight() << ", "
			//		  << "pabs:" << it->pabs() << ", "
				//	  << "energyLoss:" << it->energyLoss() << ", "
					//  << "trackId:" << it->trackId() << ", "
					  //<< "exitPoint:" << it->exitPoint() << ", "
					  //<< "localPosition:" << it->localPosition() << "\n";
					
		   //std::cout << (*it) << "\n";
		  
		  
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
	
	
	
	/* 
		////---------------------------------
		edm::Handle<std::vector<PSimHit> > simHitsTest;
		edm::InputTag tag(hitsProducer_, SubdetName_);
		std::cout << hitsProducer_ << " " << SubdetName_  << " ";
			
		e.getByToken(m_token, simHitsTest);
		//e.getByLabel(tag, simHitsTest);
		DigiEventCnt++;
		std::cout << " Finnailzed EventCnt: " << DigiEventCnt ;
			
		std::cout << "Finalzed PSimHits per event. size: " << simHitsTest->size() << " \n";
		std::cout << "New event, BCX:: " << e.bunchCrossing()
					 << " iEvent Princ ID: " << e.id()
					 << "::\n" ;
		std::vector<PSimHit> const& simHits_Test= *simHitsTest.product();	
		for (auto ps : simHits_Test)
			std::cout << ps.detUnitId()  << " Tof: "<< ps.tof() << "\n" ;
		
		std::cout << "----------------\n" ;	  
	*/

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
	// for Fbcm, we should consider LowToF, becuase Tof~9.4 ns
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

		//////---------------------------------
		//std::cout << "\nall PSimHits per event. size: " << simHits->size() << "\n";
		//std::cout << "New event, BCX:: " << iEvent.bunchCrossing()
			//		 << " iEvent Princ ID: " << iEvent.principal().id()
				//	 << "::\n" ;
		//std::vector<PSimHit> const& simHits_Test= *simHits.product();	
	//	for (auto ps : simHits_Test)
//			std::cout << ps.detUnitId() <<"\n" ;// << " Tof: "<< ps.tof() << "\n" ;
		
		//std::cout << "----------------\n" ;	  
	  
	  
	  
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
    std::vector<edm::DetSet<SiPadDigiData> > DigiDataVector;
			
    for (auto const& SiPadUnit : theFbcmGeom->SiPads()) 
	{
      std::map<int, SiPadDigiData> SiPadDigilMap; // first: channelNo(pixel) per SiPad --> by default each SiPad comprises one pixel
												  // Second: DigiDataResult per SiPad (channel)
	  
		SiPadDigiAlgo->GetDigiResults(SiPadUnit, SiPadDigilMap); // this fills out the SiPadDigilMap
	  
      	FbcmDetId SiPdetId(SiPadUnit->geographicalId().rawId());
		  
		  //std::cout << SiPdetId; 
		  //std::cout << SiPadUnit->geographicalId().rawId() << "\n"; 
		  
      edm::DetSet<SiPadDigiData> collector(SiPadUnit->geographicalId().rawId());
	    
/* 	  edm::DetSet<unsigned int> UN_Collector(SiPadUnit->geographicalId().rawId());
	  UN_Collector.data.emplace_back(SiPadUnit->geographicalId().rawId());
	  if (!UN_Collector.data.empty())
	  {
        UI_Vector.push_back(std::move(UN_Collector));
		std::cout << "UN_Collector is not empty\n"; 
	  } */
	  
	  
      //edm::DetSet<SiPadDigiSimLink> linkcollector(SiPadUnit->geographicalId().rawId());
      for (auto const& Ch_DigiRes : SiPadDigilMap) // for each pixel(channel) in SiPad
	  {
        SiPadDigiData DigiResult = Ch_DigiRes.second;
	
		 //std::cout << "\n---------------- InSide SiPad Digitizer ----------------\n";
		// std::cout << SiPdetId; 
		//std:: cout << DigiResult; 
		


  
		
		
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
      if (!collector.data.empty())
	  {
        DigiDataVector.push_back(std::move(collector));
		//std::cout << "DigiDataVector is not empty\n"; 
	  }
		
    }

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<SiPadDigiData> > output(new edm::DetSetVector<SiPadDigiData>(DigiDataVector));

    // Step D: write output to file
    iEvent.put(std::move(output), InstanceName_);
	//iEvent.put(std::move(outputOverT), "SiPad2");
	
    //if (makeDigiSimLinks_) {
      //iEvent.put(std::move(outputlink), "SiPadSimLink");
    //}
	 
  }
  
  /*
  Local3DPoint SiPadDigitizer::CMSUnits(Local3DPoint lp)
  {
	  return Local3DPoint(lp.x()*0.1,lp.y()*0.1,lp.z()*0.1); // conversion from mm to cm
	  
  }

	void SiPadDigitizer::ConvertPSimHitsToCMSUnits(std::vector<PSimHit> const& simHits_org, std::vector<PSimHit> &simHitsCMSUnit) {
		
/// Important note regarding ConvertPSimHitsToCMSUnits: 
/// This function is needed only when the “TimingSD” class is used as the parent of FbcmSD for the SensitiveDetector.
/// Because “TimingSD” stores the PSimhits with the unit of “mm” similar to G4; however, “cm” is used for CMS. 
/// The “TkAccumulatingSensitiveDetector” or “MuonSensitiveDetector”, converts and saves the PsimHits with “cm” units for the distances. 
/// Thus, if sometime another SensetiveDetector class should be used, then you shoul double check the dimensions of lengths.
/// Maybe this function would not be required. 
			
		 for(const auto& PSimIt: simHits_org) {
			 
		 //std::cout <<"entryP and ExitP Orginal: " << PSimIt.entryPoint() << PSimIt.exitPoint() << " \n"; 
			 
			 simHitsCMSUnit.push_back (PSimHit(
									CMSUnits(PSimIt.entryPoint()),
									CMSUnits(PSimIt.exitPoint()),
									PSimIt.pabs(),
									PSimIt.tof(),
									PSimIt.energyLoss(),
									PSimIt.particleType(),
									PSimIt.detUnitId(),
									PSimIt.trackId(),
									PSimIt.thetaAtEntry(),
									PSimIt.phiAtEntry(),
									PSimIt.processType()
									)
			 );
		}
		
	}


*/
}  // namespace cms

// Check if it poosible to do in other file
//using cms::SiPadDigitizer;
//DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);

#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"

using cms::SiPadDigitizer;
DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);
