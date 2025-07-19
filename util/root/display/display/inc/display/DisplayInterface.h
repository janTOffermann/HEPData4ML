#ifndef ROOT_EVENT_DISPLAY
#define ROOT_EVENT_DISPLAY

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
// #include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif


// standard lib includes
#include <vector>

using namespace std;
namespace Display{

  class EventDisplay;

  class DisplayInterface{
    public:
      DisplayInterface(TString delphesCardFilepath);
      ~DisplayInterface();

      void DisplayEvent(TString delphesCardFilepath);
      EventDisplay* GetEventDisplay(){return eventDisplay_;};

      // // Functions for adding data. These are just basically pass-throughs to eventDisplay_.
      // // TODO: Do we want to restructure this? Maybe just make
      // void AddCaloData(TString name, vector<Double_t> etaMin, vector<Double_t> etaMax, vector<Double_t> phiMin, vector<Double_t> phiMax, vector<Double_t> Eem, vector<Double_t> Ehad){
      //   eventDisplay_->AddCaloData(name,etaMin,etaMax,phiMin,phiMax,Eem,Ehad);
      // };
      // void AddJetData(TString name, vector<Double_t> a0, vector<Double_t> a1, vector<Double_t> a2, vector<Double_t> a3, Double_t radius, Bool_t useCartesian=kTRUE){
      //   eventDisplay_->AddJetData(name,a0,a1,a2,a3,radius,useCartesian);
      // };
      // void AddJetData_EPxPyPz(TString name, vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz, Double_t radius){
      //   eventDisplay_->AddJetData_EPxPyPz(name,E,px,py,pz,radius);
      // };
      // void AddJetData_PtEtaPhiM(TString name, vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, Double_t radius){
      //   eventDisplay_->AddJetData_PtEtaPhiM(name,pt,eta,phi,m,radius);
      // };

      // void AddTrackTypeData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Int_t> pdgId, Bool_t toMuon, const enum EColor color){
      //   eventDisplay_->AddTrackTypeData(name,pt,eta,phi,m,x,y,z,pdgId,toMuon,color);
      // }
      // void AddTrackData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Int_t> pdgId);
      // void AddElectronData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Double_t> charge);
      // void AddMuonData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Double_t> charge);


      void Test(TString filename);


    private:
      EventDisplay* eventDisplay_;
      TString delphesCardFilepath_;

      void TestData();



  };
}

#endif