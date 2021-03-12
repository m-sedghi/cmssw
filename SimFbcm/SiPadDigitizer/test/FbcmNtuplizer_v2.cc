
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/...
#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include <string> 

#define MAX_SIMHIT_SIZE 100
#define MAX_DIGHIT_SIZE 100
#define AREA_GROUP_LEN 8

#define TEST_IND 7

using namespace std;

class FbcmNtuplizer_v2 : public edm::EDAnalyzer {
public:
  explicit FbcmNtuplizer_v2(const edm::ParameterSet&);
  ~FbcmNtuplizer_v2();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


  // ----------member data ---------------------------
  edm::EDGetTokenT< edm::DetSetVector<SiPadDigiData> > TokenTag_;
  edm::ESHandle<FbcmGeometry> theFbcmGeom;
  
  int nbrOfDiesPerRing;

  TH1* hNEvents;
  
  TTree* theTree;
  TTree* FixedValuesTree;
  TTree* SizeTest_tree[AREA_GROUP_LEN];
  //  TTree*  Flux_Tree;


  float SensorArea;
  float SensorRho;
  int NumOfSensors;
  //float SumChargePerBx;

  int nSimParticles;
  int BxSlotCnt; 
  int nValidDigiToAs;
  int nValidDigiToTs;
  
  int sum_nSimHitsPerEvnt[AREA_GROUP_LEN];
  int sum_nDigiHitsPerEvnt[AREA_GROUP_LEN];
  
  float SimLambda[AREA_GROUP_LEN];
  float DigiLambda[AREA_GROUP_LEN];
  float SimFlux[AREA_GROUP_LEN];
  float DigiFlux[AREA_GROUP_LEN];
  
  
  
  int SimPdgIds[MAX_SIMHIT_SIZE];
  float SimPts[MAX_SIMHIT_SIZE];
  float SimCharges[MAX_SIMHIT_SIZE];
  float SimTofs[MAX_SIMHIT_SIZE];
  float SimTof_perBx[MAX_SIMHIT_SIZE];

  int DigiHitStatus[3];
  int nDigiHits[3];
  
  float DigiToAs[MAX_DIGHIT_SIZE] ;
  float DigiToTs[MAX_DIGHIT_SIZE];
  
  /*
    void resetValues(){
    sensorRho = sensorArea = SIMTotalCharge = -1.0;
    // SIMTOA->clear();
    // SIMPt->clear();
    // SIMPdgId->clear();
    // SIMCharge->clear();
    nSIMParticles=0;

    for(int i=0; i < 100 ; i++){
    SIMPdgIds[i] = 0;
    SIMPts[i] = -100;
    SIMCharges[i] = -100;
    SIMTOAs[i] = -10000;
    }

    for(int i=0; i<3; i++){
    DIGIHitStatus[i] = -100;
    DIGInHits[i] = -100;
    }
    };
  */
  void ZeroInitArrays();
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FbcmNtuplizer_v2::FbcmNtuplizer_v2(const edm::ParameterSet& iConfig) :
  TokenTag_(consumes< edm::DetSetVector<SiPadDigiData> >(iConfig.getParameter<edm::InputTag>("FbcmDigiTag")))
{
  edm::Service<TFileService> fs;
  //theTree = fs->make<TTree>( iConfig.getParameter< string >("TreeName").c_str() , "all hits" );
  //theTree = fs->make<TTree>( "Pu" , "all hits" );
  FixedValuesTree = fs->make<TTree>( "GeometryInfo" , "FixedValues" );
  hNEvents = fs->make<TH1I>("hNEvents" , "" , 1 , 0 , 1 );

  Int_t split = 99;
  for (int i = 0 ; i < AREA_GROUP_LEN ; i++)
    SizeTest_tree[i] = fs->make<TTree>( ("SensorSize_"+std::to_string(i)).c_str() , ("Group"+std::to_string(i)).c_str() , split );
  
  //Flux_Tree = fs->make<TTree>( "FluxInfo" , "FluxInfo" );
  
  //  resetValues();
  ZeroInitArrays();
  Int_t bsize = 32000; //default

  for (int i = 0 ; i < 8 ; i++) {
    SizeTest_tree[i]->Branch("SensorRho" , &SensorRho , "SensorRho/f[8,22]", bsize );
    SizeTest_tree[i]->Branch("SensorArea" , &SensorArea , "SensorArea/f[0.01,1.0]" , bsize  );
    //SizeTest_tree[i]->Branch("SimSumChargePerBx" , &SumChargePerBx , "SimSumChargePerBx/F");
	  
    SizeTest_tree[i]->Branch("nSimParticles" , &nSimParticles);
    SizeTest_tree[i]->Branch("SimPdgId" , SimPdgIds , "SimPdgId[nSimParticles]/I" , bsize );
    SizeTest_tree[i]->Branch("SimPt" , SimPts , "SimPt[nSimParticles]/F" , bsize );
    SizeTest_tree[i]->Branch("SimCharge" , SimCharges , "SimCharge[nSimParticles]/F", bsize );
    SizeTest_tree[i]->Branch("SimTof" , SimTofs , "SimTof[nSimParticles]/F", bsize ); 
	  
	  
    SizeTest_tree[i]->Branch("SimTof_perBx" , SimTof_perBx , "SimTof_perBx[nSimParticles]/F", bsize ); 
	  
    SizeTest_tree[i]->Branch("BxSlotCnt" , &BxSlotCnt);
    SizeTest_tree[i]->Branch("DigiHitStatus" , DigiHitStatus , "DigiHitStatus[BxSlotCnt]/I", bsize );
    SizeTest_tree[i]->Branch("nDigiHits" , nDigiHits , "nDigiHits[BxSlotCnt]/I" , bsize ); 
	  
    SizeTest_tree[i]->Branch("nValidDigiToAs" , &nValidDigiToAs);
    SizeTest_tree[i]->Branch("nValidDigiToTs" , &nValidDigiToTs);
    SizeTest_tree[i]->Branch("DigiToA" , DigiToAs , "DigiToAs[nValidDigiToAs]/F" , bsize );
    SizeTest_tree[i]->Branch("DigiToT" , DigiToTs , "DigiToTs[nValidDigiToTs]/F", bsize );
	  
  }

  // Flux_Tree->Branch("sum_nDigiHitsTest1" , &SimLambda[TEST_IND] , "sum_nDigiHitsTest1/F" );
  // Flux_Tree->Branch("sum_nDigiHitsTest2" , &DigiLambda[TEST_IND] , "sum_nDigiHitsTest2/F" );

  // Flux_Tree->Branch("sum_nSimHits" , sum_nSimHitsPerEvnt , "sum_nSimHits[8]/I" );
  // Flux_Tree->Branch("sum_nDigiHits" , sum_nDigiHitsPerEvnt , "sum_nDigiHits[8]/I" );
  // Flux_Tree->Branch("SimLambda" , SimLambda , "SimLambda[8]/F" );
  // Flux_Tree->Branch("DigiLambda" , DigiLambda , "DigiLambda[8]/F" );
  // Flux_Tree->Branch("SimFlux" , SimFlux , "SimFlux[8]/F" );
  // Flux_Tree->Branch("DigiFlux" , DigiFlux , "DigiFlux[8]/F" ); 



  /*
    theTree->Branch("SensorRho" , &sensorRho , "SensorRho/f[8,22]" );
    theTree->Branch("SensorArea" , &sensorArea , "SensorArea/f[0.01,1.0]" );
    theTree->Branch("SIMTotalCharge" , &SIMTotalCharge , "SIMTotalCharge/f");
    theTree->Branch("nSIMParticles" , &nSIMParticles);
    theTree->Branch("SIMPdgId" , SIMPdgIds , "SIMPdgId[nSIMParticles]/I" );
    theTree->Branch("SIMPt" , SIMPts , "SIMPt[nSIMParticles]/F" );
    theTree->Branch("SIMCharge" , SIMCharges , "SIMCharge[nSIMParticles]/F");
    theTree->Branch("SIMTOA" , SIMTOAs , "SIMTOA[nSIMParticles]/F"); 
    theTree->Branch("DIGInHits" , DIGInHits , "DIGInHits[3]/I" ); 
    theTree->Branch("DIGIHitStatus" , DIGIHitStatus , "DIGIHitStatus[3]/I");
  
  */
	  
}


FbcmNtuplizer_v2::~FbcmNtuplizer_v2()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FbcmNtuplizer_v2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  hNEvents->Fill( 0.5 );

  edm::Handle< edm::DetSetVector<SiPadDigiData> > handle;
  iEvent.getByToken(TokenTag_,handle);


  
  for (int i=0 ; i< AREA_GROUP_LEN ; i++) {
    sum_nSimHitsPerEvnt[i] = 0;
    sum_nDigiHitsPerEvnt[i] = 0;
  }
  
  int SiDieId;
  int SensorGroupIndex;
  for (edm::DetSetVector<SiPadDigiData>::const_iterator itDetSet = handle->begin() ;  itDetSet < handle->end() ; ++itDetSet) {
    // each Det in fbcmDets:
    //std::cout<< "\nRaw DetID:  " << itDetSet->detId() << "\n"; // each detId --> multiple Data, but we have only one data
    //char ss;
    //std::cin >> ss;
	  
    if( itDetSet->size() != 1 )
      cout << "size == " << itDetSet->size() << endl; 
    for (std::vector<SiPadDigiData>::const_iterator Digidata_it = itDetSet->begin() ; Digidata_it < itDetSet->end() ; ++Digidata_it) {
      //const SiPadDigiData Digidata= *(Digidata_it);
      //std::cout << *(Digidata_it) ;
	  
      SiDieId=Digidata_it->SiliconDieIndex() ;
      SensorGroupIndex = SiDieId % nbrOfDiesPerRing;
      //std::cout << SiDieId << " --> " << SensorGroupIndex << "\n";
	  
	
	  
      SensorRho=Digidata_it->Radius() ;
      SensorArea=Digidata_it->Area() ;
	  
      //SimBxChargeSum
	  	  
      //resetValues();

      //this->SIMTotalCharge = Digidata_it->ChargeSum();
	 
      //	  float SimHitTime;
      nSimParticles = 0 ;
      //SumChargePerBx = 0.0;
      for(auto sim : Digidata_it->CahrgePsimVector()){
	SimPdgIds[nSimParticles] = sim.second.ParticleType() ;
	SimPts[nSimParticles] = sim.second.Pabs() ;
	SimTofs[nSimParticles] = sim.second.Tof() ;
	SimCharges[nSimParticles] = sim.first;
			
	SimTof_perBx[nSimParticles] = sim.second.Tof() - 25.0 * sim.second.BunchCrossing() ;
	//SimHitTime = sim.second.time() - 25.0 * sim.second.BunchCrossing() ;

	nSimParticles++;
      }
	  
      sum_nSimHitsPerEvnt[SensorGroupIndex] += nSimParticles;
	  
      nValidDigiToAs=0;
      nValidDigiToTs=0;
      BxSlotCnt=0; // just to make sure starting from non-zero index if in the case of using BxSlotNo < 0. 
      for(auto hit : Digidata_it->BxSlotHitAnalysisVector() ) {	  
	DigiHitStatus[BxSlotCnt] = hit.Bx_HitStatusInt();
	nDigiHits[BxSlotCnt] = hit.nbrOfRecognizedHitsInBx(); 
	BxSlotCnt++;
		
	for (auto ToaTot : hit.TotToaVectort()) {
	  if ( (ToaTot.ToAToAStatusInt() == ToAStatus::FullyWithinBx || ToaTot.ToAToAStatusInt() == ToAStatus::WithinBx_LastsAfter ) ) {
				
	    DigiToAs[nValidDigiToAs++] = ToaTot.ToA();
	    if (ToaTot.IsToTValid())
	      DigiToTs[nValidDigiToTs++] = ToaTot.ToT();
	  }
	}
		
      }
      sum_nDigiHitsPerEvnt[SensorGroupIndex] += nValidDigiToAs;


      //theTree->Fill();
      SizeTest_tree[SensorGroupIndex]->Fill();
      //std::cout << "detected raw ID: " << fbdetId.rawId() << "  : ";
      //std::cout << fbdetId;
					
    }
    //std::cout << "----------------------------\n";
		
		
    //histo->Fill(itDetSet->detId());
  }
  
  // for (int i =0 ; i < AREA_GROUP_LEN ; i++) {
  //   SimLambda[i] = (float)sum_nSimHitsPerEvnt[i]/nbrOfSenosrsVect[i]/7.0;
  //   DigiLambda[i] = (float)sum_nDigiHitsPerEvnt[i]/nbrOfSenosrsVect[i]/3.0;
  //   SimFlux[i] = (float) sum_nSimHitsPerEvnt[i]/TotalSensAreaPerGroup[i]/7.0; // devided by 7 due to 7 bunch crossings
  //   DigiFlux[i] = (float) sum_nDigiHitsPerEvnt[i]/TotalSensAreaPerGroup[i]/3.0; // devided by 3 because we had consdered only 3 BXs, i.e. {0,1,2},  in the Digi
	
  // }
  
  // Flux_Tree->Fill(); 

}


// ------------ method called once each job just before starting event loop  ------------
void 
FbcmNtuplizer_v2::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FbcmNtuplizer_v2::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void FbcmNtuplizer_v2::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
  iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom); 
  const std::vector<const FbcmStationGeom*> AllStatitons = theFbcmGeom->Stations();
  const std::vector<FbcmSiPadGeom const*> allSiPadGeoms = theFbcmGeom->SiPads();
  nbrOfDiesPerRing = AllStatitons[0]->NumOfDiesPerRing(); 


  //std::cout << "Total Num of Stations: " << AllStatitons.size() << "\n";

  int SiDieId=0;
  int SensorGroupIndex=0;
  float padX, padY , padRho;

  FixedValuesTree->Branch("SensorGroupIndex" , &SensorGroupIndex );
  FixedValuesTree->Branch("SensorX" , &padX );
  FixedValuesTree->Branch("SensorY" , &padY );
  FixedValuesTree->Branch("SensorRho" , &padRho );

//  for(int side=0 ; side < 0 ; side++) // since "AllStatitons" includes all stations in both ends (+z/-z), so this loop make a wrong number. 
//    for(auto station: AllStatitons){
//      for(auto die: station->SiliconDies()){
//	SiDieId=die->id().SiliconDie();
//	SensorGroupIndex = SiDieId % nbrOfDiesPerRing;

//	for (auto siPad : die->SiPads() )
	for(const auto& siPad : allSiPadGeoms)
	{
	  std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();

	  padX = SiPadDimension.first;
	  padY = SiPadDimension.second;
	  // we could also store only the Area instead of padX and padY
	  padRho = siPad->surface().position().perp(); 
	  
	  SiDieId=siPad->id().SiliconDie();
      SensorGroupIndex = SiDieId % nbrOfDiesPerRing;
	  
	  FixedValuesTree->Fill();
	}
//      }
//    }
  
}

void FbcmNtuplizer_v2::ZeroInitArrays(){

  // for (int i=0 ; i< AREA_GROUP_LEN ; i++) {
  //   nbrOfSenosrsVect[i]=0;
  //   EachSensorAreaVect[i]=0.0;
  //   TotalSensAreaPerGroup[i]=0.0;
  // }
}

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  FbcmNtuplizer_v2::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  FbcmNtuplizer_v2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  FbcmNtuplizer_v2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FbcmNtuplizer_v2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FbcmNtuplizer_v2);
