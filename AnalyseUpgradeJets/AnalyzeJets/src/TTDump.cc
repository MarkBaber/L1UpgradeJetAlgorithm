#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TLorentzVector.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h" 
// #include "DataFormats/L1Trigger/interface/L1JetParticle.h"
// #include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

// #include "DataFormats/JetReco/interface/JetID.h" 
// #include "RecoJets/JetProducers/interface/JetIDHelper.h"
// #include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/SLHC/interface/EtaPhiContainer.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"

#include <algorithm>  // for sorting

using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;

const Double_t PI = 3.141592654;


// Switch for TT-level info
#define TT_INFO


// ========================================
// class declaration
// ========================================

class TTDump : public edm::EDAnalyzer {

public:
  explicit TTDump(const edm::ParameterSet&);
  ~TTDump();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);




  // ----------member data ---------------------------
  ParameterSet conf_;
  edm::InputTag m_jetCollection;
  edm::InputTag m_jetIDMap;


  // TT variables
  int TT_Et, TT_Eta, TT_Phi;
  std::vector <int> TT_Et_Arr, TT_Eta_Arr, TT_Phi_Arr;
  int TT_index; // TT index in the 1D vector (start top left tower (eta,phi) = (-28,1) scanning lines of constant phi left to right
  //                            to bottom right TT   (eta,phi) = (28,72)

  // Should reserve space in the CONSTRUCTOR!!!!
  std::vector<int> TT_ArrayNew,muonPx;


  // ************************************************************************************
  // *                                       ROOT                                       *
  // ************************************************************************************

  
  // TT
  TTree *ttTree;



  // TT tree branches
  // ****************************************

  std::vector<int>      ttE;
  std::vector<int>      ttH;
  std::vector<int>      ttiEta;
  std::vector<int>      ttiPhi;
  std::vector<int>      ttEcalFG;
  std::vector<int>      ttHcalFG;


  


};





TTDump::TTDump(const edm::ParameterSet& iConfig): conf_(iConfig){




  // Allocate enough space to store all the TT in the barrel and endcaps
  TT_ArrayNew.resize(4032);


  //now do what ever initialization is needed
  edm::Service<TFileService> tfs;




  // TT tree
  // ****************************************

  ttTree = tfs->make<TTree>("TT","TT");

  ttTree->Branch( "E",      &ttE);
  ttTree->Branch( "H",      &ttH);
  ttTree->Branch( "iEta",   &ttiEta);
  ttTree->Branch( "iPhi",   &ttiPhi);
  ttTree->Branch( "EcalFG", &ttEcalFG);
  ttTree->Branch( "HcalFG", &ttHcalFG);


}



TTDump::~TTDump()
{
}



// ------------ method called once each job just before starting event loop  ------------
void TTDump::beginJob(){
}








// **********************************************************************
// *                              analyze()                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************


void TTDump::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool evValid = true;



  
  // TT collection
  edm::Handle<l1slhc::L1CaloTowerCollection> caloTowers;
  iEvent.getByLabel(conf_.getParameter<edm::InputTag>("CalorimeterTowers"), caloTowers);
  if(!caloTowers.isValid()){
    edm::LogWarning("MissingProduct") << conf_.getParameter<edm::InputTag>("CalorimeterTowers") << std::endl;
    evValid = false;
  }






  // ------------------------------------------------------------------------------------------------------------------------------------------------------
  // -                                                              Fill Histograms                                                                       -
  // ------------------------------------------------------------------------------------------------------------------------------------------------------


 
  if(!evValid){
    std::cout << "Invalid event\n";	
  }
  else{   // valid event




    
    // ****************************************************************************************************
    // *                                        TT distributions                                          *
    // ****************************************************************************************************

    
#ifdef TT_INFO


    ttE.clear();
    ttH.clear();
    ttiEta.clear();
    ttiPhi.clear();
    ttEcalFG.clear();
    ttHcalFG.clear();
    
    



    for( l1slhc::L1CaloTowerCollection::const_iterator lTT_It = caloTowers->begin() ; 
	 lTT_It != caloTowers->end() ; ++lTT_It ){

      // ****************************************
      // *    Load the calorimeter tower data   *
      // ****************************************
      int E      = lTT_It->E();
      int H      = lTT_It->H();
      //      int EplusH = E + H;
      int iEta   = lTT_It->iEta();
      int iPhi   = lTT_It->iPhi();
      int EcalFG = lTT_It->EcalFG();
      int HcalFG = lTT_It->HcalFG();
    

      // Restrict to central TTs
      if (abs(iEta) > 28)
	continue;


      ttE.push_back(E);
      ttH.push_back(H);
      ttiEta.push_back(iEta);
      ttiPhi.push_back(iPhi);
      ttEcalFG.push_back(EcalFG);
      ttHcalFG.push_back(HcalFG);


    }


    ttTree->Fill();

#endif






  } // end valid event

} // end analyser














// ------------ method called once each job just after ending the event loop  ------------
void TTDump::endJob(){
}





// ------------ method called when starting to processes a run  ------------
void 
TTDump::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTDump::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TTDump::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TTDump::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTDump::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}









// ------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------




//define this as a plug-in
DEFINE_FWK_MODULE(TTDump);
