/*
  *
  * This code is largely copied from the Delphes framework's "display" library,
  * as due to some technicalities it was difficult for me to simply write a
  * pure extension of the class (due to how loading headers via PyROOT works,
  * and my not knowing the best way to reliably point ROOT's gInterpreter to
  * the relevant Delphes include path). It has some modifications and additions
  * of features to fit my particular use case, but the core functionality
  * is retained from the Delphes code.
  * For completeness of documentation, I am including the original top comment
  * below. - Jan
  *
  * Delphes: a framework for fast simulation of a generic collider experiment
  * Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
 *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
 *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EventDisplay_h
#define EventDisplay_h

#include <vector>
#include <map>

#include "RQ_OBJECT.h"
#include "Rtypes.h"
#include "TEveElement.h"
#include "TEveCalo.h"
#include "TEveTrackPropagator.h"
#include "TColor.h"

class TAxis;
class TChain;
class TGHtml;
class TGStatusBar;
// class DelphesDisplay;
class Delphes3DGeometry;
class DelphesBranchBase;
class DelphesHtmlSummary;
class DelphesPlotSummary;
class ExDelphesDisplay;

using namespace std;

namespace Display{

  class BranchBase;

  class EventDisplay
  {
    RQ_OBJECT("EventDisplay")
  public:
    EventDisplay();
    EventDisplay(Delphes3DGeometry &det3D);
    ~EventDisplay();
    void EventChanged(Int_t); // *SIGNAL*

    // Functions for adding data.
    // Note that some (electron, muon) are variations of the AddTrackData() function.
    void AddCaloData(TString name, vector<Double_t> etaMin, vector<Double_t> etaMax, vector<Double_t> phiMin, vector<Double_t> phiMax, vector<Double_t> Eem, vector<Double_t> Ehad);
    void AddJetData(TString name, vector<Double_t> a0, vector<Double_t> a1, vector<Double_t> a2, vector<Double_t> a3, Double_t radius, Bool_t useCartesian=kTRUE, const enum EColor color = kYellow);
    void AddJetData_EPxPyPz(TString name, vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz, Double_t radius, const enum EColor color = kYellow);
    void AddJetData_PtEtaPhiM(TString name, vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, Double_t radius, const enum EColor color = kYellow);

    // We add a few different kinds of "Track-type" data: electons, muons, and actual tracks.
    // To manage this, we have an "AddTrackTypeData()" function, that will then be used by
    // specialized functions for the different data types.
    void AddTrackTypeData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Int_t> pdgId, Bool_t toMuon, const enum EColor color);
    void AddTrackData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Int_t> pdgId);
    void AddElectronData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Double_t> charge);
    void AddMuonData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, vector<Double_t> charge);
    void AddPhotonData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m);
    void AddMETData(TString name,vector<Double_t> pt, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m);
    void AddGenParticleData(TString name,vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz, vector<Double_t> xProd, vector<Double_t> yProd, vector<Double_t> zProd, vector<Double_t> xDecay, vector<Double_t> yDecay, vector<Double_t> zDecay, vector<Bool_t> stable, vector<Int_t> pdgId); // TODO: Eventually support displaced particle production vertices!
    void AddVertexData(TString name, vector<Double_t> x, vector<Double_t> y, vector<Double_t> z, const enum EColor color = kBlue);

    // Add jet constituent data, in the format of CaloData.
    // This puts things in a special collection!
    void AddJetConstituentCaloData(TString name, vector<Double_t> etaMin, vector<Double_t> etaMax, vector<Double_t> phiMin, vector<Double_t> phiMax, vector<Double_t> Eem, vector<Double_t> Ehad);

    // Functions for creating data containers.
    // These can be explicitly called, or they will be invoked by the
    // functions for adding data, as needed.
    // NOTE: These used to create TClonesArrays, which are now unused.
    //       Can use those by modifying the "Add___Data()" functions to
    //       populate those arrays instead, and then use BranchBase::ReadBranch(),
    //       which would have to be reimplemented.
    void AddCaloContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kBlack);
    void AddJetContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kYellow, Bool_t truthLevel = kFALSE);
    void AddTrackContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kCyan, Bool_t truthLevel = kFALSE, Bool_t toMuon = kFALSE);
    void AddMuonContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kRed, Bool_t truthLevel = kFALSE){
      AddTrackContainer(name,/*Int_t expectedSize,*/ color, truthLevel, kTRUE);
    }
    void AddMETContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kViolet);
    void AddVertexContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kBlue);

    void AddJetConstituentCaloContainer(TString name, /*Int_t expectedSize=10,*/ const enum EColor color = kBlack);


    void Display(Int_t eventNumber=-1);
    void SetMode(Int_t mode){mode_ = mode;};
    void SetColorMode(Int_t mode){colorMode_ = mode;};

    void SetJetPtMin(Double_t val){jetPtMin_ = val;};
    void SetTrackPtMin(Double_t val){trackPtMin_ = val;};
    void SetTruthParticlePtMin(Double_t val){truthParticlePtMin_ = val;};

  protected:
    void update_html_summary();
    void make_gui();
    void load_event();

    void BuildDetector(Delphes3DGeometry &det3D);
    void ResetContainers();

    void AddContainersToScene();
    // void CreateUnifiedCaloDataContainer();
    void InitScene();

    void InitCalorimeters();
    void ScaleLego();
    void AddJetsToCaloDisplay();
    void AddRefsToCaloDisplay();

    enum EColor FetchJetColor();
    void ResetJetColor(){jetColorIndex_ = 0;};

    // Need this unfortunate function due to some pretty irritating
    // TEveCaloLego behaviour... - Jan
    vector<Double_t> ConvertEtaPhiToCaloLego(Double_t eta, Double_t phi);
    Double_t ConvertRadiusToCaloLego(Double_t radius);

    // Configuration and global variables.
    Int_t mode_ = 0;
    Int_t event_id_;
    Int_t event_id_tmp_;
    Double_t tkRadius_, totRadius_, tkHalfLength_, muHalfLength_, bz_;
    TAxis *etaAxis_, *phiAxis_;
    TChain *chain_;
    map<TString, BranchBase *> elements_;
    map<TString, BranchBase *> jetMapElements_;
    ExDelphesDisplay *delphesDisplay_;
    DelphesHtmlSummary *htmlSummary_;
    TGHtml *gHtml_;
    DelphesPlotSummary *plotSummary_;
    TGStatusBar *fStatusBar_;
    TEveElementList *geometry_ = 0;
    vector<TString> addedContainers_ = {};

    vector<TEveTrackPropagator *> trkProp_ = {};

    TEveCaloLego *lego_ = 0;
    vector<TString> jetContainers_ = {}; // keep track of which containers represent jets
    map<TString, vector<Float_t>> jetEta_ = {};
    map<TString, vector<Float_t>> jetPhi_ = {};
    map<TString, Float_t> jetRadius_ = {};
    Double_t legoEtaScale_ = 1.;
    Double_t legoPhiScale_ = 1.;

    Int_t colorMode_ = 0;
    vector<enum EColor> jetColors_ = {};
    Int_t jetColorIndex_ = 0;

    Double_t jetPtMin_ = 0.;
    Double_t trackPtMin_ = 0.;
    Double_t truthParticlePtMin_ = 0.;


    // gui controls
  public:
    void Fwd();

    void Bck();

    void PreSetEv(char *ev);

    void GoTo();

    void InitSummaryPlots();

    void DisplayProgress(Int_t p);
  };
}
#endif //EventDisplay_h

