#ifndef HIST_BOOKER
#define HIST_BOOKER




#include <map>
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"





// Object for storing histogram binning parameters
struct histParam{

  int    bins;  // Number of bins
  double low;   // Range low
  double high;  // Range high

  // Initialise object
  histParam( int _bins, double _low, double _high ):bins(_bins),low(_low),high(_high){};

};

// ********************************************************************************
// *                           Histogram binning parameters                       *
// ********************************************************************************
// prefix 'p' denotes these are booking parameters

   const histParam pDeltaPT(121,   -60.5, 60.5);

   const histParam pResponse(200,   -0.01, 3.99);
   const histParam pDeltaPTRel(121,   -1.0125, 2.0125);
   const histParam pDeltaEta(121,  -0.605, 0.605);
   const histParam pDeltaPhi(121,  -0.605, 0.605);
   const histParam pDeltaR(100,  -0.005, 0.995);

/*    const histParam pL1PT(600,   -0.5, 599.5); */
/*    const histParam pOffPT(440,  -0.5, 439.5); */
   const histParam pL1PT(200,   -0.5, 199.5);
   const histParam pOffPT(200,  -0.5, 199.5);
   const histParam pL1PTLarge(600,   -0.5, 599.5);
   const histParam pOffPTLarge(440,  -0.5, 439.5);
   const histParam pOffEta(61, -3.05, 3.05);
   const histParam pOffPhi(63, -3.15, 3.15);
  
   const histParam pL1HT( 151, -5, 1505);
   const histParam pOffHT(151, -5, 1505);
   const histParam pL1MHT( 51, -5, 505);
   const histParam pOffMHT(51, -5, 505);

//  const double PUbin(81),      PUlow(-0.5),      PUhigh(80.5);
//  const double JetBin(50),     JetLow(-0.5),     JetHigh(49.5);
const histParam pPU(81, -0.5, 80.5);
// Is this even needed??? Duplication of PU?????????????????????????????????????????
const histParam pJetMultiplicity(81, -0.5, 80.5);

// TT multiplicity
const histParam pEFired(101,-0.5,100.5);
const histParam pHFired(101,-0.5,100.5);
const histParam pEorHFired(201,-0.5,200.5);
const histParam pEandHFired(101,-0.5,100.5);

// TT Energy
const histParam pETotal(1501,-0.5,1500.5);
const histParam pHTotal(1501,-0.5,1500.5);
const histParam pEplusHTotal(1501,-0.5,1500.5);
const histParam pLETotal(201,-0.5,200.5);
const histParam pLHTotal(201,-0.5,200.5);
const histParam pLEplusHTotal(201,-0.5,200.5);


const histParam pEMax(201,-0.5,200.5);
const histParam pHMax(201,-0.5,200.5);




const histParam pMedianE(51, -0.5, 50.5);
const histParam pMedianH(51, -0.5, 50.5);
const histParam pMedianEplusH(51, -0.5, 50.5);


const histParam piEta(57,-28.5,28.5);




class histBooker{

  public:

    // constructor
    histBooker(std::map< TString, TH1*> *h1D, std::map< TString, TH1*> *h2D):histogram1D(h1D), histogram2D(h2D){};

    // Member functions
    void book1D(         TString histName, TFileDirectory f, TString histTitle, histParam xParam );
    void book1DPtEtaPhi( TString histName, TFileDirectory f, TString histTitle );
    void bookTProfile(   TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam );
    void book2D(         TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam );
    void book2DTProf(    TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam );


  private:
    // List of sub directories (key) and their currently booked histograms for input validation
    std::map< TString, std::list < TString> > hist1DList;
    std::map< TString, std::list < TString> > hist2DList;
  //    std::list<TString> histSubDir2DList;

    // Pointers to histogram containers
    std::map< TString, TH1*> *histogram1D;
    std::map< TString, TH1*> *histogram2D;



};



#endif
