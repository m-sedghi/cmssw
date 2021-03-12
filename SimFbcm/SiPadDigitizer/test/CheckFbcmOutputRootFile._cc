
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


class FbcmOutputRootFileTester : public edm::EDAnalyzer {
   public:
      explicit FbcmOutputRootFileTester(const edm::ParameterSet&);
      ~FbcmOutputRootFileTester();

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
      const edm::InputTag FbcmDigiTag_;
	  edm::EDGetTokenT< edm::DetSetVector<SiPadDigiData> > TokenTag_;
	  edm::ESHandle<FbcmGeometry> theFbcmGeom;
	  
      TH1D *histo;

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
FbcmOutputRootFileTester::FbcmOutputRootFileTester(const edm::ParameterSet& iConfig) :
		FbcmDigiTag_(iConfig.getParameter<edm::InputTag>("FbcmDigiTag")),
		TokenTag_(consumes< edm::DetSetVector<SiPadDigiData> >(iConfig.getParameter<edm::InputTag>("FbcmDigiTag")))
{
   //now do what ever initialization is needed
  //edm::Service<TFileService> fs;
  //histo = fs->make<TH1D>("pt" , "Pt" , 50 , 0 , 50 );
}


FbcmOutputRootFileTester::~FbcmOutputRootFileTester()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FbcmOutputRootFileTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace std;

	//iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom);

   edm::Handle< edm::DetSetVector<SiPadDigiData> > handle;
   ///both "getByToken" and "getByLabel" work well
   iEvent.getByToken(TokenTag_,handle);
   //iEvent.getByLabel(FbcmDigiTag_,handle);
   //iEvent.getByLabel(edm::InputTag("simFbcmDigis","SiPad"),handle);
   
	for (edm::DetSetVector<SiPadDigiData>::const_iterator itDetSet = handle->begin() ;  itDetSet < handle->end() ; ++itDetSet) {
		std::cout<< "\nRaw DetID:  " << itDetSet->detId() << "\n"; // each detId --> multiple Data
		for (std::vector<SiPadDigiData>::const_iterator Det_It = itDetSet->begin() ; Det_It < itDetSet->end() ; ++Det_It) {
			//const SiPadDigiData tmp= *(Det_It);
			std::cout << *(Det_It) ;
			FbcmDetId fbdetId( (int)(Det_It->SideIndex()),
					(int)(Det_It->StationIndex()),
					(int)(Det_It->SiliconDieIndex()),
					(int)(Det_It->SiPadIndex()) );
					std::cout << "detected raw ID: " << fbdetId.rawId() << "  : ";
					std::cout << fbdetId;
					
		}
		std::cout << "----------------------------\n";
		
		
		//histo->Fill(itDetSet->detId());
	}


/* #ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
    */
}


// ------------ method called once each job just before starting event loop  ------------
void 
FbcmOutputRootFileTester::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FbcmOutputRootFileTester::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
FbcmOutputRootFileTester::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
	iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom);
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
FbcmOutputRootFileTester::endRun(edm::Run const&, edm::EventSetup const&)
{
		std::cout << "Run finished\n";
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
FbcmOutputRootFileTester::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
FbcmOutputRootFileTester::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FbcmOutputRootFileTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FbcmOutputRootFileTester);
