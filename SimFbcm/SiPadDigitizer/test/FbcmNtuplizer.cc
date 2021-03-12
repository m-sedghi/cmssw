
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

using namespace std;

class FbcmNtuplizer : public edm::EDAnalyzer {
public:
  explicit FbcmNtuplizer(const edm::ParameterSet&);
  ~FbcmNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


  // ----------member data ---------------------------
  edm::EDGetTokenT< edm::DetSetVector<SiPadDigiData> > TokenTag_;
  
  TTree* theTree;

  float sensorArea;
  float sensorRho;

  float SIMTotalCharge;

  int nSIMParticles;
  int SIMPdgIds[100];
  float SIMPts[100];
  float SIMCharges[100];
  float SIMTOAs[100];

  int DIGIHitStatus[3];
  int DIGInHits[3];

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
FbcmNtuplizer::FbcmNtuplizer(const edm::ParameterSet& iConfig) :
  TokenTag_(consumes< edm::DetSetVector<SiPadDigiData> >(iConfig.getParameter<edm::InputTag>("FbcmDigiTag")))
{
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>( iConfig.getParameter< string >("TreeName").c_str() , "all hits" );

  resetValues();

  theTree->Branch("SensorRho" , &sensorRho , "SensorRho/f[8,22]" );
  theTree->Branch("SensorArea" , &sensorArea , "SensorArea/f[0.01,1.0]" );
  theTree->Branch("SIMTotalCharge" , &SIMTotalCharge , "SIMTotalCharge/f");
  theTree->Branch("nSIMParticles" , &nSIMParticles );
  theTree->Branch("SIMPdgId" , SIMPdgIds , "SIMPdgId[nSIMParticles]/I" );
  theTree->Branch("SIMPt" , SIMPts , "SIMPt[nSIMParticles]/F" );
  theTree->Branch("SIMCharge" , SIMCharges , "SIMCharge[nSIMParticles]/F");
  theTree->Branch("SIMTOA" , SIMTOAs , "SIMTOA[nSIMParticles]/F");

  theTree->Branch("DIGInHits" , DIGInHits , "DIGInHits[3]/I" );
  theTree->Branch("DIGIHitStatus" , DIGIHitStatus , "DIGIHitStatus[3]/I");
}


FbcmNtuplizer::~FbcmNtuplizer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FbcmNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  edm::Handle< edm::DetSetVector<SiPadDigiData> > handle;
  iEvent.getByToken(TokenTag_,handle);
   
  for (edm::DetSetVector<SiPadDigiData>::const_iterator itDetSet = handle->begin() ;  itDetSet < handle->end() ; ++itDetSet) {
    //std::cout<< "\nRaw DetID:  " << itDetSet->detId() << "\n"; // each detId --> multiple Data
    if( itDetSet->size() != 1 )
      cout << "size == " << itDetSet->size() << endl; 
    for (std::vector<SiPadDigiData>::const_iterator Det_It = itDetSet->begin() ; Det_It < itDetSet->end() ; ++Det_It) {
      //const SiPadDigiData tmp= *(Det_It);
      //std::cout << *(Det_It) ;
      resetValues();

      this->sensorRho = Det_It->Radius();
      this->sensorArea = Det_It->Area();
      this->SIMTotalCharge = Det_It->ChargeSum();

      for(auto sim : Det_It->CahrgePsimVector()){
	SIMPdgIds[nSIMParticles] = sim.second.ParticleType() ;
	SIMPts[nSIMParticles] = sim.second.Pabs() ;
	SIMTOAs[nSIMParticles] = sim.second.Tof() ;

	SIMCharges[nSIMParticles] = sim.first;
	nSIMParticles++;
      }

      for(auto hit : Det_It->BxSlotHitAnalysisVector() ){
	DIGIHitStatus[hit.BxSlotNo()] = hit.Bx_HitStatusInt();
	DIGInHits[hit.BxSlotNo()] = hit.nbrOfRecognizedHitsInBx();
      }


      theTree->Fill();
      //std::cout << "detected raw ID: " << fbdetId.rawId() << "  : ";
      //std::cout << fbdetId;
					
    }
    //std::cout << "----------------------------\n";
		
		
    //histo->Fill(itDetSet->detId());
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
FbcmNtuplizer::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FbcmNtuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  FbcmNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  FbcmNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  FbcmNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  FbcmNtuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FbcmNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FbcmNtuplizer);
