// Original Author:  Mark Baber Imperial College, London
//
// Greater than seed (GTS) filtering


// system include files
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "SimDataFormats/SLHC/interface/L1TowerJet.h"
#include "SimDataFormats/SLHC/interface/L1TowerJetFwd.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include <iostream>
#include <fstream>

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"



#include <algorithm>

//
// class declaration
//


using namespace l1slhc;
using namespace edm;
using namespace std;
using namespace reco;
using namespace l1extra;

bool TowerJetRankDescending ( l1slhc::L1TowerJet jet1, l1slhc::L1TowerJet jet2 ){      return ( jet1.p4().Pt() > jet2.p4().Pt() ); }



int getDeltaiEta( int iEta1, int iEta2 ){
  int iEta1Index = iEta1;
  int iEta2Index = iEta2;
  if (iEta1 > 0){ iEta1Index--; }
  if (iEta2 > 0){ iEta2Index--; }

  return iEta1Index - iEta2Index;
}

int getDeltaiPhi( int iPhi1, int iPhi2 ){

  int deltaiPhi = iPhi1 - iPhi2;
  if ( deltaiPhi > 36){
    deltaiPhi -= 72;
  }
  else if ( deltaiPhi < -36){
    deltaiPhi += 72;
  }

  return deltaiPhi;
}






class L1TowerJetGTSFilter : public edm::EDProducer {
   public:
      explicit L1TowerJetGTSFilter(const edm::ParameterSet&);
      ~L1TowerJetGTSFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      virtual double Median( std::vector< double > aVec );
    
       // ----------member data ---------------------------
      ParameterSet conf_;
    
      // Limit of the number of jets to be retained in filtering
      int mNumOfOutputJets;  

      std::auto_ptr < l1slhc::L1TowerJetCollection > mOutputCollection;


      // Tower geometry converter 
      TriggerTowerGeometry mTowerGeo;


  int seedThreshold;

  
  double protoGlobalRho;


};

//
// constructors and destructor
//

L1TowerJetGTSFilter::L1TowerJetGTSFilter(const edm::ParameterSet& iConfig):
  conf_(iConfig),
  mNumOfOutputJets( iConfig.getParameter<uint32_t>("NumOfOutputJets") ),
  seedThreshold( iConfig.getParameter<double>("SeedEnergyThreshold") )
{
  produces< double >("ProtoGlobalRho");
  produces< double >("SeedTTE");
  produces< double >("SeedTTRegion");
  produces< l1slhc::L1TowerJetCollection >("GTSFilteredTowerJets");
  seedThreshold *= 2; // Convert to 2 GeV units
}

L1TowerJetGTSFilter::~L1TowerJetGTSFilter()
{
}




// ------------ method called to produce the data  ------------
void
L1TowerJetGTSFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    bool evValid = true;

    edm::Handle< L1TowerJetCollection > PreFilteredJets;
    iEvent.getByLabel(conf_.getParameter<edm::InputTag>("PreFilteredJets"), PreFilteredJets);
    if(!PreFilteredJets.isValid()){
        evValid=false;
    }

    if( !evValid ) {
        throw cms::Exception("MissingProduct") << conf_.getParameter<edm::InputTag>("preFiltJets")      
	  				       << std::endl; 
    }
    else{

      std::auto_ptr< double > outputProtoGlobalRho(new double());
      auto_ptr< L1TowerJetCollection > filteredJets(new L1TowerJetCollection());
      l1slhc::L1TowerJetCollection AllPreFilteredJets;
      l1slhc::L1TowerJetCollection unsortedJets;


      // Number of jets currently retained
      int lJetCounter(0);
      int iEtaMin = -28;

      
      // Protojet rho container
      std::vector<double> protoGlobalRhos;
      

      // Make jets to the edge of the calorimeter
      for (l1slhc::L1TowerJetCollection::const_iterator PreFilt_It = PreFilteredJets->begin(); PreFilt_It != PreFilteredJets->end(); ++PreFilt_It ){

	// Calculate protojet rho
	double protoRho = PreFilt_It->E()/ PreFilt_It->JetRealArea();
	protoGlobalRhos.push_back( protoRho );



	L1TowerJet tempJet = *PreFilt_It;

	int iEtaMax =  28 - (PreFilt_It->JetSize() -1);
	
	int iEta = PreFilt_It->iEta();
	//	int iPhi = PreFilt_It->iPhi();
	int halfJetSize = PreFilt_It->JetSize()/2;


	AllPreFilteredJets.insert( tempJet.iEta(), tempJet.iPhi(), tempJet );

	// ********************************************************************************
	// *                               Make small jets                                *
	// ********************************************************************************
	if ( (iEta == iEtaMin) || (iEta == iEtaMax) ){
	 	  
	  int sign(1);
	  if (iEta == iEtaMin){
	    sign = -1;
	  }

	  for ( int deltaIEta = 1; deltaIEta <= halfJetSize; deltaIEta++){

	    const L1CaloTowerRefVector jetConstituents = tempJet.getConstituents();

	    int delIEta     = sign*deltaIEta;
	    int newIEta     = iEta + delIEta;
	    int newEtaEdge(newIEta);
	    if (iEta == iEtaMin){
	      newEtaEdge = newIEta + PreFilt_It->JetSize() - 1; 
	    }

	    // Correct the jet area
	    double areaCorrection = ( PreFilt_It->JetSize() * deltaIEta );
	    double corrJetArea    = PreFilt_It->JetArea()  - areaCorrection;
	    tempJet.SetJetArea( corrJetArea );	    

	    double areaRealCorrection(0);
	    //	    std::cout << "iEta = " << newIEta << "\n";
	    for ( int delIEtaOut = 1; delIEtaOut <= abs(delIEta); ++delIEtaOut ){
	      
	      int iEtaOutside = newEtaEdge - sign*delIEtaOut;

	      //	      std::cout << "\tiEtaOut = " << iEtaOutside << "\n";
	      areaRealCorrection += ( PreFilt_It->JetSize() * mTowerGeo.towerArea( iEtaOutside, 1 ) );

	    }
	    double corrJetRealArea    = PreFilt_It->JetRealArea() - areaRealCorrection;
	    tempJet.SetJetRealArea( corrJetRealArea );

// 	    std::cout << corrJetRealArea << "\t" << tempJet.JetRealArea() << "\t" << PreFilt_It->JetRealArea() << "\t" 
// 		      << corrJetArea     << "\t" << tempJet.JetArea()     << "\t" << PreFilt_It->JetArea() << "\n";

// 	    std::cout << "Jet area      = "   << PreFilt_It->JetArea()     << "\tCorrected jetArea     = " << corrJetArea
// 		      << "\nJet real area = " << PreFilt_It->JetRealArea() << "\tCorrected jetRealArea = " << corrJetRealArea << "\n\n";



	    for ( L1CaloTowerRefVector::const_iterator lConstituentIt = jetConstituents.begin() ; 
		  lConstituentIt != jetConstituents.end(); ++lConstituentIt ){
	      
	      int TTiEta = (*lConstituentIt)->iEta();
	      int TTiPhi = (*lConstituentIt)->iPhi();
	      
	      // Remove TTs from the tower jet 
	      if (      (iEta == iEtaMin) && ( TTiEta > newEtaEdge ) ){ 
		tempJet.removeConstituent( TTiEta, TTiPhi );
	      }
	      else if ( (iEta == iEtaMax) && ( TTiEta < newEtaEdge ) ){

		tempJet.removeConstituent( TTiEta, TTiPhi );

	      }

	    }

	    // Calculate the energy weighted eta and phi and the centrality of the jet 
	    tempJet.calculateWeightedJetCenter();
	    tempJet.calculateCentrality();

	    //	    if ( originalE != tempJet.E() ){
// 	      std::cout << iEta << "\tnewIEta = " << newIEta << "\toldE = " << originalE << "\tE = " << tempJet.E() 
// 			<< "\tEta = " << tempJet.WeightedEta() << "\tPhi = " << tempJet.WeightedPhi() 
// 			<< "\tConst = " << nConst << "\n";
	      //	    }	    

	    tempJet.setP4( math::PtEtaPhiMLorentzVector( tempJet.E(), tempJet.WeightedEta(), tempJet.WeightedPhi(), 0. ) );
	    tempJet.setiEta( newIEta );
	    
	    AllPreFilteredJets.insert( newIEta, tempJet.iPhi(), tempJet );
	  }
	}

      }


      //      for (l1slhc::L1TowerJetCollection::const_iterator PreFilt_It = PreFilteredJets->begin(); PreFilt_It != PreFilteredJets->end(); ++PreFilt_It ){
      for (l1slhc::L1TowerJetCollection::const_iterator PreFilt_It = AllPreFilteredJets.begin(); PreFilt_It != AllPreFilteredJets.end(); ++PreFilt_It ){



	const L1CaloTowerRefVector jetConstituents = PreFilt_It->getConstituents();
	
	int iEta = PreFilt_It->iEta();
	int iPhi = PreFilt_It->iPhi();
	int halfJetSize = PreFilt_It->JetSize()/2;

	int seediEta = iEta + halfJetSize;
	// Check if jumped the iEta = 1 transition
	if ( (seediEta > -1 ) && (iEta < 0) ){
	  seediEta++;
	}

	int seediPhi = iPhi + halfJetSize;
	// Check if jumped the iPhi = 1 transition	
	if ( seediPhi > 72 ){
	  seediPhi -= 72;
	}


	// Check if there is a jet seed
	L1TowerJet tempJet = *PreFilt_It;
	int seedE(0), seedH(0), seedEplusH(0);
	bool seedExists(false);
	for ( L1CaloTowerRefVector::const_iterator lConstituentIt = jetConstituents.begin() ; lConstituentIt != jetConstituents.end(); ++lConstituentIt ){

	  int TTiEta = (*lConstituentIt)->iEta();
          int TTiPhi = (*lConstituentIt)->iPhi();

	  // Seed TT exists
	  if ( ( TTiEta == seediEta ) && ( TTiPhi == seediPhi ) ){
	    seedE = (*lConstituentIt)->E();
	    seedH = (*lConstituentIt)->H();
	    seedEplusH = seedE + seedH;
	    if ( seedEplusH < seedThreshold ){ break; }
	    seedExists = true;
	    break;
	  }

	}
	if (!seedExists){ continue; } // No seed was found




	//	std::cout << "\n\nJet:   iEta = " << iEta   << "\tiPhi = " << iPhi << "\n";

	// Check if jet fires the 'greater than seed' veto
 	bool GTSVeto = false;
  	for ( L1CaloTowerRefVector::const_iterator lConstituentIt = jetConstituents.begin() ; lConstituentIt != jetConstituents.end(); ++lConstituentIt ){
	  
	  
	  int TTiEta = (*lConstituentIt)->iEta();
	  int TTiPhi = (*lConstituentIt)->iPhi();

	  // Exclude seed TT from check
 	  if ( (iEta == seediEta) && (iPhi == seediPhi) ){ continue;}
	  
	  int E      = (*lConstituentIt)->E();
	  int H      = (*lConstituentIt)->H();
	  int EplusH = E + H;

// 	  std::cout << "Const: iEta = " << TTiEta << "\tiPhi = " << TTiPhi
// 		    << "\tE = "         << E      << "\tH = "    << H 
// 		    << "\tEplusH = "    << EplusH << "\n";
	      

	  // Check TT against veto mask
	  if ( EplusH > seedEplusH ){
		GTSVeto = true;
	  }
	  else if( EplusH == seedEplusH ){

	    // Calculate deltaiEta, deltaiPhi
	    int deltaiEta        = getDeltaiEta( seediEta, TTiEta );
	    int deltaiPhi        = getDeltaiPhi( seediPhi, TTiPhi );
	    int deltaiEtaiPhiSum = deltaiEta + deltaiPhi;

	    if ( (deltaiEta < 0) && (deltaiEtaiPhiSum <= 0) ){
		GTSVeto = true;
		break;
	    }
	    else if ( (deltaiEta >= 0) && (deltaiEtaiPhiSum < 0) ){
		GTSVeto = true;
		break;
	    }

 	  }
	  
  	} // End jet constituent loop


// 	if ( GTSVeto ){
// 	  std::cout << "Seed: iEta = " << seediEta   << "\tiPhi = " << seediPhi 
// 		    << "\tE = "        << seedE      << "\tH = "    << seedH
// 		    << "\tEplusH = "   << seedEplusH << "\n";
// 	  std::cout << "Veto = " <<GTSVeto << "\n";	  
// 	}

 	// Store non-vetoed jets
	if ( !GTSVeto ){
	  unsortedJets.insert( tempJet.iEta(), tempJet.iPhi(), tempJet );
	  
 	  lJetCounter++;
  	  if( lJetCounter >= mNumOfOutputJets ){ break; }  // Restrict number of jets retained in filtering  
	  
	}
    
      } // End jet loop         
      

      // Sort and store jets
      std::sort( unsortedJets.begin(),         unsortedJets.end(),         TowerJetRankDescending );
      for ( L1TowerJetCollection::const_iterator sorted_It = unsortedJets.begin(); sorted_It != unsortedJets.end(); ++sorted_It ){
	//	int iEta = sorted_It.iEta()
	  filteredJets->insert( sorted_It->iEta(), sorted_It->iPhi(), *sorted_It );
      }



      protoGlobalRho = Median( protoGlobalRhos );
      *outputProtoGlobalRho = protoGlobalRho;
      iEvent.put( outputProtoGlobalRho, "ProtoGlobalRho");


      // Store the GTS filtered jets
      iEvent.put( filteredJets, "GTSFilteredTowerJets" );




    } // End valid event



    
    
}





double 
L1TowerJetGTSFilter::Median( std::vector< double > aVec ){ 
 
  if ( aVec.size() == 0 ){ 
    return -1; 
  } 
 
  // Order vector collection                                                                                                                                
  sort( aVec.begin(), aVec.end() ); 
 
  double median(0); 
  int size = aVec.size();                                                                                                                                   
  int halfSize = size/2;                                                                                                                                    
  if( size == 0 ){                                                                                                                                          
    median = 0;                                                                                                                                             
  }                                                                                                                                                         
  else if( size == 1 ){                                                                                                                                     
    median = aVec[0];                                                                                                                                       
  }                                                                                                                                                         
  else if( size%2 == 0 ){ 
    // Even number of entries, take average of the values around center 
    median = ( aVec[ halfSize - 1 ] + aVec[ halfSize ] ) * 0.5;                                                                                             
  }                                                                                                                                                         
  else{                                                                                                                                                     
    // Odd number of entries, halfSize is central element                                                             
                                                                                                                                                           
    median = aVec[ halfSize ]; 
  } 
 
  return median; 
}                     







// ------------ method called once each job just before starting event loop  ------------
void 
L1TowerJetGTSFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TowerJetGTSFilter::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
L1TowerJetGTSFilter::beginRun(edm::Run&, edm::EventSetup const&)
{    
}

// ------------ method called when ending the processing of a run  ------------
void 
L1TowerJetGTSFilter::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
L1TowerJetGTSFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
L1TowerJetGTSFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TowerJetGTSFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TowerJetGTSFilter);
