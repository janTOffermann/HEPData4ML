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

#include <display/EventDisplay.h>
#include <display/BranchElement.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>

#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TEveArrow.h"
#include "TEveBrowser.h"
// #include "TEveElement.h"
#include "TEveEventManager.h"
#include "TEveGeoNode.h"
#include "TEveJetCone.h"
#include "TEveManager.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveTrans.h"
#include "TEveViewer.h"
#include "TGButton.h"
#include "TGHtml.h"
#include "TGNumberEntry.h"
#include "TGProgressBar.h"
#include "TGStatusBar.h"
#include "TGTextEntry.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRootBrowser.h"
#include "TRootEmbeddedCanvas.h"
#include "TSystem.h"

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include "display/Delphes3DGeometry.h"
#include "display/DelphesBranchElement.h"
#include "display/CaloData.h"
#include "display/DelphesDisplay.h"
#include "display/DelphesHtmlSummary.h"
#include "display/DelphesPlotSummary.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootConfReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

namespace Display{
  EventDisplay::EventDisplay(){
    event_id_ = 0;
    tkRadius_ = 1.29;
    totRadius_ = 2.0;
    tkHalfLength_ = 3.0;
    muHalfLength_ = 6.0;
    bz_ = 3.8;
    // chain_ = new TChain("Delphes");
    treeReader_ = 0;
    delphesDisplay_ = 0;
    etaAxis_ = 0;
    phiAxis_ = 0;
  }

  EventDisplay::~EventDisplay(){
    // delete chain_;
    // for(auto entry : calo3d_) delete entry;
    // for(auto entry : calo3d_lego_) delete entry;
  }

  void EventDisplay::EventChanged(Int_t e){
    if(e != event_id_){
      event_id_ = e;
      Emit("EventChanged(Int_t)", e);
      load_event();
    }
  }

  void EventDisplay::BuildDetector(Delphes3DGeometry &det3D){
    TGeoManager *geom = gGeoManager;
    tkRadius_ = det3D.getTrackerRadius();
    totRadius_ = det3D.getDetectorRadius();
    tkHalfLength_ = det3D.getTrackerHalfLength();
    muHalfLength_ = det3D.getDetectorHalfLength();
    bz_ = det3D.getBField();
    etaAxis_ = det3D.getCaloAxes().first;
    phiAxis_ = det3D.getCaloAxes().second;
    TGeoVolume *top = det3D.getDetector(false);
    geom->SetTopVolume(top);
    geometry_ = new TEveElementList("Geometry");
    TObjArray *nodes = top->GetNodes();
    TIter itNodes(nodes);
    TGeoNode *nodeobj;
    TEveGeoTopNode *node;
    while((nodeobj = (TGeoNode *)itNodes.Next()))
    {
      node = new TEveGeoTopNode(gGeoManager, nodeobj);
      node->UseNodeTrans();
      geometry_->AddElement(node);
    }
  }

  EventDisplay::EventDisplay(Delphes3DGeometry &det3D){
    event_id_ = 0;
    tkRadius_ = 1.29;
    totRadius_ = 2.0;
    tkHalfLength_ = 3.0;
    bz_ = 3.8;
    // chain_ = new TChain("Delphes");
    treeReader_ = 0;
    delphesDisplay_ = 0;

    // initialize the application
    TEveManager::Create(kTRUE, "IV");
    fStatusBar_ = gEve->GetBrowser()->GetStatusBar();

    // build the detector
    BuildDetector(det3D);
  }

  void EventDisplay::Display(Int_t eventNumber){

    if(eventNumber > -1) event_id_ = eventNumber;

    // viewers and scenes, plus adding the data containers to the scene.
    InitScene();

    // the GUI: control panel, summary tab
    make_gui();
    gSystem->ProcessEvents();

    // currently gives a crash
    load_event();


    gEve->Redraw3D(kTRUE);
  }

  void EventDisplay::AddContainersToScene(){
    for(map<TString, BranchBase *>::iterator element = elements_.begin(); element != elements_.end(); ++element)
    {
      TString name = element->first;
      auto it = find(addedContainers_.begin(), addedContainers_.end(), name);
      if (it!=addedContainers_.end()) continue;

      BranchElement<TEveTrackList> *item_v1 = dynamic_cast<BranchElement<TEveTrackList> *>(element->second);
      BranchElement<TEveElementList> *item_v2 = dynamic_cast<BranchElement<TEveElementList> *>(element->second);
      if(item_v1) gEve->AddElement(item_v1->GetContainer());
      if(item_v2) gEve->AddElement(item_v2->GetContainer());
      addedContainers_.push_back(name);
    }
  }

  CaloData* CombineCaloData(const std::vector<CaloData*>& caloObjects) {
      if (caloObjects.empty()) return nullptr;

      // Create new combined CaloData with same structure
      CaloData* combined = new CaloData(2);
      combined->RefSliceInfo(0).Setup("ECAL", 0.1, kRed);
      combined->RefSliceInfo(1).Setup("HCAL", 0.1, kBlue);
      combined->IncDenyDestroy();
      Int_t combinedTowerCounter = 0;

      // Merge data from all CaloData objects
      for (const auto& caloObj : caloObjects) {
          if (!caloObj) continue;

          if (auto caloDataVec = dynamic_cast<TEveCaloDataVec*>(caloObj)) {

            // Get the geometry
            vector<Float_t> etaMin, etaMax, phiMin, phiMax;
            for (auto entry : caloDataVec->GetCellGeom()){
              etaMin.push_back(entry.EtaMin());
              etaMax.push_back(entry.EtaMax());
              phiMin.push_back(entry.PhiMin());
              phiMax.push_back(entry.PhiMax());
            }

            // loop over slices
            Int_t nSlices = caloDataVec->GetNSlices();
            Int_t nCells = caloDataVec->GetNCells();

            // loop over towers/cells
            for(Int_t i = 0; i < nCells; i++){
                combined->AddTower(etaMin.at(i),etaMax.at(i),phiMin.at(i),phiMax.at(i));
              // loop over slices
              for(Int_t j = 0; j < nSlices; j++){
                vector<Float_t> sliceVals = caloDataVec->GetSliceVals(j); // across all towers (unfortunate given how the loop is constructed!)
                combined->FillSlice(j,sliceVals.at(i));
              }
                combinedTowerCounter++;
            }
          }
      }
      combined->DataChanged();
      return combined;
  }

  // void EventDisplay::CreateUnifiedCaloDataContainer(){
  //   vector<CaloData*> caloDataVec = {};
  //   for(map<TString, BranchBase *>::iterator data = elements_.begin(); data != elements_.end(); ++data){
  //     if(TString(data->second->GetClassName()).Contains("CaloData")){
  //       caloDataVec.push_back(dynamic_cast<BranchElement<CaloData> *>((data->second))->GetContainer());
  //     }
  //   }
  //   caloDataVec.clear();
  //   CaloData* combinedCaloData = CombineCaloData(caloDataVec);

  //   TString name = "UnifiedCaloData";
  //   BranchElement<CaloData> *clist = new BranchElement<CaloData>(combinedCaloData, name, kBlack);
  //   // clist->GetContainer()->SetEtaBins(etaAxis_);
  //   // clist->GetContainer()->SetPhiBins(phiAxis_);
  //   elements_[name] = clist;

  //   // BranchElement<CaloData>* br = dynamic_cast<BranchElement<CaloData> *>(elements_[name]);
  //   // br->SetContainer(combinedCaloData);
  //   return;
  // }

  void EventDisplay::InitScene(){

    // Add the actual data containers to the scene.
    AddContainersToScene();

    // TODO: Create a unified CaloData container, to use for the viewing.
    //       The current function is broken!
    // CreateUnifiedCaloDataContainer();

    delphesDisplay_ = new DelphesDisplay;
    gEve->AddGlobalElement(geometry_);
    delphesDisplay_->ImportGeomRPhi(geometry_);
    delphesDisplay_->ImportGeomRhoZ(geometry_);
    // find the first calo data and use that to initialize the calo display
    // cout << "Starting check for CaloData containers, for lego view." << endl;
    for(map<TString, BranchBase *>::iterator data = elements_.begin(); data != elements_.end(); ++data)
    {

      // TODO: Find a way to deal with *multiple* CaloData containers.
      if(TString(data->second->GetClassName()).Contains("CaloData")){
        // // TODO: Temporary disable of EFlowTracks
        // if(data->first.Contains("Track")) continue;
        // if(data->first.Contains("Photon")) continue;
        // cout << "Adding to Lego display: " << data->first << endl;

        CaloData *container = dynamic_cast<BranchElement<CaloData> *>((data->second))->GetContainer();
        assert(container);
        TEveCalo3D *calo3d = new TEveCalo3D(container);
        calo3d->SetBarrelRadius(tkRadius_);
        calo3d->SetEndCapPos(tkHalfLength_);
        gEve->AddGlobalElement(calo3d);
        delphesDisplay_->ImportCaloRPhi(calo3d);
        delphesDisplay_->ImportCaloRhoZ(calo3d);
        TEveCaloLego *lego = new TEveCaloLego(container);
        lego->InitMainTrans();
        //    lego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
        lego->RefMainTrans().SetScale(100, 100, TMath::Pi());
        lego->SetAutoRebin(kFALSE);
        lego->Set2DMode(TEveCaloLego::kValSizeOutline);
        delphesDisplay_->ImportCaloLego(lego);
        // calo3d_.push_back(calo3d);
        // calo3d_lego_.push_back(lego);
        break; // would like to draw multiple ones, but they don't line up nicely
      }
    }
  }

  void EventDisplay::AddCaloData(TString name,
    vector<Double_t> etaMin, vector<Double_t> etaMax,
    vector<Double_t> phiMin, vector<Double_t> phiMax,
    vector<Double_t> Eem, vector<Double_t> Ehad){

    if(elements_.find(name) == elements_.end()) AddCaloContainer(name);
    Size_t nElements = etaMin.size();

    BranchElement<CaloData>* br = dynamic_cast<BranchElement<CaloData> *>(elements_[name]);

    for(Size_t i = 0; i < nElements; i++){
      br->GetContainer()->AddTower(etaMin.at(i),etaMax.at(i),phiMin.at(i),phiMax.at(i));
      br->GetContainer()->FillSlice(0, Eem.at(i));
      br->GetContainer()->FillSlice(1, Ehad.at(i));
      // cout << "Filled eta=(" << etaMin.at(i) << ", " << etaMax.at(i) << "), phi=(" << phiMin.at(i) << ", " << phiMax.at(i) << ") with (" << Eem.at(i) << ", " << Ehad.at(i) << ")" << endl;
    }
    br->GetContainer()->DataChanged();
  }

  void EventDisplay::AddJetData(TString name,
    vector<Double_t> a0, vector<Double_t> a1,
    vector<Double_t> a2, vector<Double_t> a3,
    Double_t radius, Bool_t useCartesian,
    const enum EColor color){

    if(useCartesian) AddJetData_EPxPyPz(name,a0,a1,a2,a3,radius,color);
    else AddJetData_PtEtaPhiM(name,a0,a1,a2,a3,radius,color);
  }

  void EventDisplay::AddJetData_EPxPyPz(TString name,
    vector<Double_t> E, vector<Double_t> px,
    vector<Double_t> py, vector<Double_t> pz,
    Double_t radius,
    const enum EColor color){

    if(elements_.find(name) == elements_.end()) AddJetContainer(name,color);
    Size_t nElements = E.size();

    BranchElement<TEveElementList>* br = dynamic_cast<BranchElement<TEveElementList> *>(elements_[name]);
    TEveJetCone *eveJetCone;
    ROOT::Math::PxPyPzEVector vec;

    Int_t counter = 0;
    for(Size_t i = 0; i < nElements; i++){
      vec.SetCoordinates(px.at(i),py.at(i),pz.at(i),E.at(i));
      eveJetCone = new TEveJetCone();
      eveJetCone->SetTitle(Form("jet [%d]: Pt=%f, Eta=%f, \nPhi=%f, M=%f", counter, vec.Pt(), vec.Eta(), vec.Phi(), vec.M()));
      eveJetCone->SetName(Form("jet [%d]", counter++));
      eveJetCone->SetMainTransparency(60);
      eveJetCone->SetLineColor(br->GetColor());
      eveJetCone->SetFillColor(br->GetColor());
      eveJetCone->SetCylinder(tkRadius_ - 10, tkHalfLength_ - 10); // TODO: Why the -10?
      eveJetCone->SetPickable(kTRUE);
      eveJetCone->AddEllipticCone(vec.Eta(), vec.Phi(), 0.5 * radius, 0.5 * radius); // check if radius or 1/2 radius?
      br->GetContainer()->AddElement(eveJetCone);
    }
  }

  // TODO: Templating?
  void EventDisplay::AddJetData_PtEtaPhiM(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m,
    Double_t radius,
    const enum EColor color){

    if(elements_.find(name) == elements_.end()) AddJetContainer(name,color);
    Size_t nElements = pt.size();

    BranchElement<TEveElementList>* br = dynamic_cast<BranchElement<TEveElementList> *>(elements_[name]);
    TEveJetCone *eveJetCone;
    ROOT::Math::PtEtaPhiMVector vec;

    Int_t counter = 0;
    for(Size_t i = 0; i < nElements; i++){
      vec.SetCoordinates(pt.at(i),eta.at(i),phi.at(i),m.at(i));
      eveJetCone = new TEveJetCone();
      eveJetCone->SetTitle(Form("jet [%d]: Pt=%f, Eta=%f, \nPhi=%f, M=%f", counter, vec.Pt(), vec.Eta(), vec.Phi(), vec.M()));
      eveJetCone->SetName(Form("jet [%d]", counter++));
      eveJetCone->SetMainTransparency(60);
      eveJetCone->SetLineColor(br->GetColor());
      eveJetCone->SetFillColor(br->GetColor());
      eveJetCone->SetCylinder(tkRadius_ - 10, tkHalfLength_ - 10); // TODO: Why the -10?
      eveJetCone->SetPickable(kTRUE);
      eveJetCone->AddEllipticCone(vec.Eta(), vec.Phi(), 0.5 * radius, 0.5 * radius); // check if radius or 1/2 radius?
      br->GetContainer()->AddElement(eveJetCone);
    }
  }

  void EventDisplay::AddTrackTypeData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m,
    vector<Double_t> x, vector<Double_t> y,
    vector<Double_t> z, vector<Int_t> pdgId,
    Bool_t toMuon, const enum EColor color){
    if(elements_.find(name) == elements_.end()){
      if(toMuon) AddMuonContainer(name,color);
      else AddTrackContainer(name,color);
    }
    Size_t nElements = pt.size();

    BranchElement<TEveTrackList>* br = dynamic_cast<BranchElement<TEveTrackList> *>(elements_[name]);
    ROOT::Math::PtEtaPhiMVector vec;
    ROOT::Math::XYZTVector posVec;
    TEveTrack *eveTrack;

    TEveTrackPropagator *trkProp = br->GetContainer()->GetPropagator();
    trkProp->SetMagField(0., 0., br->GetTrackingVolume().at(2));
    trkProp->SetMaxR(br->GetTrackingVolume().at(0));
    trkProp->SetMaxZ(br->GetTrackingVolume().at(1));

    for(Int_t i = 0; i < nElements; i++){
    vec.SetCoordinates(pt.at(i),eta.at(i),phi.at(i),m.at(i));
    posVec.SetCoordinates(x.at(i)/10.,y.at(i)/10.,z.at(i)/10.,0.); // note the conversion
      TParticle pb(
        pdgId.at(i), 1, 0, 0, 0, 0,
        vec.Px(), vec.Py(), vec.Pz(), vec.E(),
        posVec.X(), posVec.Y(), posVec.Z(), posVec.T()
      );
      eveTrack = new TEveTrack(&pb, i, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), i));
      eveTrack->SetIndex(i); // TODO: Is this OK?
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(br->GetContainer());
      br->GetContainer()->AddElement(eveTrack);
      eveTrack->SetLineColor(br->GetColor());
      eveTrack->MakeTrack();
    }
  }

  void EventDisplay::AddTrackData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m,
    vector<Double_t> x, vector<Double_t> y,
    vector<Double_t> z, vector<Int_t> pdgId){

    AddTrackTypeData(name,
      pt,eta,phi,m,
      x,y,z,pdgId,kFALSE,kBlue
    );
  }

  void EventDisplay::AddElectronData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m,
    vector<Double_t> x, vector<Double_t> y,
    vector<Double_t> z, vector<Double_t> charge){

    vector<Int_t> pdgId(charge.size());
    std::transform(charge.begin(), charge.end(), pdgId.begin(),
                   [](Int_t q) { return q > 0. ? -11 : 11; });

    AddTrackTypeData(name,
      pt,eta,phi,m,
      x,y,z,pdgId,kFALSE,kGreen
  );
  }

  void EventDisplay::AddMuonData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m,
    vector<Double_t> x, vector<Double_t> y,
    vector<Double_t> z, vector<Double_t> charge){

    vector<Int_t> pdgId(charge.size());
    std::transform(charge.begin(), charge.end(), pdgId.begin(),
                   [](Int_t q) { return q > 0. ? -13 : 13; });

    AddTrackTypeData(name,
      pt,eta,phi,m,
      x,y,z,pdgId,kTRUE,kRed
    );
  }

  void EventDisplay::AddPhotonData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m){
    if(elements_.find(name) == elements_.end()){
      AddTrackContainer(name,kOrange);
    }
    Size_t nElements = pt.size();

    BranchElement<TEveTrackList>* br = dynamic_cast<BranchElement<TEveTrackList> *>(elements_[name]);
    ROOT::Math::PtEtaPhiMVector vec;
    TEveTrack *eveTrack;

    TEveTrackPropagator *trkProp = br->GetContainer()->GetPropagator();
    trkProp->SetMagField(0., 0., br->GetTrackingVolume().at(2));
    trkProp->SetMaxR(br->GetTrackingVolume().at(0));
    trkProp->SetMaxZ(br->GetTrackingVolume().at(1));

    for(Int_t i = 0; i < nElements; i++){
    vec.SetCoordinates(pt.at(i),eta.at(i),phi.at(i),m.at(i));
      TParticle pb(
       22, 1, 0, 0, 0, 0,
        vec.Px(), vec.Py(), vec.Pz(), vec.E(),
        0., 0., 0., 0. // assuming always handling prompt photons
      );
      eveTrack = new TEveTrack(&pb, i, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), i));
      eveTrack->SetIndex(i);
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(br->GetContainer());
      eveTrack->SetLineStyle(7);
      br->GetContainer()->AddElement(eveTrack);
      eveTrack->SetLineColor(br->GetColor());
      eveTrack->MakeTrack();
    }
  }

  void EventDisplay::AddMETData(TString name,
    vector<Double_t> pt, vector<Double_t> eta,
    vector<Double_t> phi, vector<Double_t> m){

    if(elements_.find(name) == elements_.end()){
      AddTrackContainer(name,kMagenta);
    }
    Size_t nElements = pt.size();

    BranchElement<TEveTrackList>* br = dynamic_cast<BranchElement<TEveTrackList> *>(elements_[name]);
    ROOT::Math::PtEtaPhiMVector vec;
    TEveTrack *eveTrack;

    TEveTrackPropagator *trkProp = br->GetContainer()->GetPropagator();
    trkProp->SetMagField(0., 0., br->GetTrackingVolume().at(2));
    trkProp->SetMaxR(br->GetTrackingVolume().at(0));
    trkProp->SetMaxZ(br->GetTrackingVolume().at(1));

    for(Int_t i = 0; i < nElements; i++){
    vec.SetCoordinates(pt.at(i),eta.at(i),phi.at(i),m.at(i));
      TParticle pb(
       13, 1, 0, 0, 0, 0,
        vec.Px(), vec.Py(), vec.Pz(), vec.E(),
        0., 0., 0., 0.
      );
      eveTrack = new TEveTrack(&pb, i, trkProp);
      eveTrack->SetName("Missing Et");
      eveTrack->SetStdTitle();
      eveTrack->SetRnrPoints(0);
      eveTrack->SetMarkerColor(kMagenta);
      eveTrack->SetMarkerStyle(4);
      eveTrack->SetMarkerSize(2.);
      eveTrack->SetLineWidth(2.);
      eveTrack->SetLineStyle(7);
      br->GetContainer()->AddElement(eveTrack);
      eveTrack->SetLineColor(br->GetColor());
      eveTrack->MakeTrack();
    }
  }

  void EventDisplay::AddGenParticleData(TString name,
    vector<Double_t> E, vector<Double_t> px,
    vector<Double_t> py, vector<Double_t> pz,
    vector<Double_t> x, vector<Double_t> y,
    vector<Double_t> z,
    vector<Int_t> pdgId){

    if(elements_.find(name) == elements_.end()){
      AddTrackContainer(name,kCyan);
    }
    Size_t nElements = E.size();

    BranchElement<TEveTrackList>* br = dynamic_cast<BranchElement<TEveTrackList> *>(elements_[name]);
    ROOT::Math::PxPyPzEVector vec;
    ROOT::Math::XYZTVector posVec;
    TEveTrack *eveTrack;

    TEveTrackPropagator *trkProp = br->GetContainer()->GetPropagator();
    trkProp->SetMagField(0., 0., br->GetTrackingVolume().at(2));
    trkProp->SetMaxR(br->GetTrackingVolume().at(0));
    trkProp->SetMaxZ(br->GetTrackingVolume().at(1));

    for(Int_t i = 0; i < nElements; i++){
    vec.SetCoordinates(px.at(i),py.at(i),pz.at(i),E.at(i));
    posVec.SetCoordinates(x.at(i)/10.,y.at(i)/10.,z.at(i)/10.,0.); // note the conversion
      TParticle pb(
        pdgId.at(i), 1, 0, 0, 0, 0,
        vec.Px(), vec.Py(), vec.Pz(), vec.E(),
        posVec.X(), posVec.Y(), posVec.Z(), posVec.T()
      );
      eveTrack = new TEveTrack(&pb, i, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), i));
      eveTrack->SetIndex(i);
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(br->GetContainer());
      if(eveTrack->GetCharge() == 0) eveTrack->SetLineStyle(7);
      br->GetContainer()->AddElement(eveTrack);
      eveTrack->SetLineColor(br->GetColor());
      eveTrack->MakeTrack();
    }
  }

  void EventDisplay::AddCaloContainer(TString name, /*Int_t expectedSize,*/ const enum EColor color){
    // TClonesArray* arr = new TClonesArray(name,expectedSize); // will tuck this away in the BranchElement.
    BranchElement<CaloData> *clist = new BranchElement<CaloData>(name, color);
    clist->GetContainer()->SetEtaBins(etaAxis_);
    clist->GetContainer()->SetPhiBins(phiAxis_);
    elements_[name] = clist;
  }

  void EventDisplay::AddJetContainer(TString name, /*Int_t expectedSize,*/ const enum EColor color, Bool_t truthLevel){
    // TClonesArray* arr = new TClonesArray(name,expectedSize); // will tuck this away in the BranchElement.
    BranchElement<TEveElementList> *elist = new BranchElement<TEveElementList>(name, color);
    if(truthLevel){
      elist->GetContainer()->SetRnrSelf(false);
      elist->GetContainer()->SetRnrChildren(false);
    }
    elist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
    elements_[name] = elist;
  }

  void EventDisplay::AddTrackContainer(TString name, /*Int_t expectedSize,*/ const enum EColor color, Bool_t truthLevel, Bool_t toMuon){
    // TClonesArray* arr = new TClonesArray(name,expectedSize); // will tuck this away in the BranchElement.
    BranchElement<TEveTrackList> *tlist = new BranchElement<TEveTrackList>(name, color);
    if(truthLevel){
        tlist->GetContainer()->SetRnrSelf(false);
        tlist->GetContainer()->SetRnrChildren(false);
    }

    if(toMuon)  tlist->SetTrackingVolume(totRadius_, muHalfLength_, bz_);
    else tlist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
    elements_[name] = tlist;
  }

  void EventDisplay::AddMETContainer(TString name, /*Int_t expectedSize,*/ const enum EColor color){
    // TClonesArray* arr = new TClonesArray(name,expectedSize); // will tuck this away in the BranchElement.
    BranchElement<TEveElementList> *elist = new BranchElement<TEveElementList>(name, color);
    elist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_); // TODO: Configurable length? Could go out to muon system.
    elements_[name] = elist;
  }

  void EventDisplay::ResetContainers(){
    for(map<TString, BranchBase *>::iterator data = elements_.begin(); data != elements_.end(); ++data){
      data->second->Reset(); // TODO: Check this?
    }
  }

  void EventDisplay::load_event()
  {
    // Load event specified in global event_id_.
    // The contents of previous event are removed.

    // NOTE: Currently disabling the ability to load different events,
    //       since with my current modifications this class is no longer
    //       "plugged into" an input file. The user just passes lists of
    //       objects via the API, but these currently all go into the same event.

    // safety
    if(event_id_ < 0) return;
    // if(event_id_ >= treeReader_->GetEntries() || event_id_ < 0) return;

    // message
    fStatusBar_->SetText(Form("Loading event %d.", event_id_), 1);
    gSystem->ProcessEvents();

    // // clear the previous event
    // TODO: This is broken with current usage, displays an event and immediately erases it
    // gEve->GetViewers()->DeleteAnnotations();
    // ResetContainers();

    // // Load selected branches with data from specified event
    // treeReader_->ReadEntry(event_id_);
    // for(map<TString, BranchBase *>::iterator data = elements_.begin(); data != elements_.end(); ++data)
    // {
    //   data->second->ReadBranch();
    // }

    // // update display
    TEveElement *top = (TEveElement *)gEve->GetCurrentEvent();
    delphesDisplay_->DestroyEventRPhi();
    delphesDisplay_->ImportEventRPhi(top);
    delphesDisplay_->DestroyEventRhoZ();
    delphesDisplay_->ImportEventRhoZ(top);

    update_html_summary();
    // plotSummary_->FillEvent();
    // plotSummary_->Draw();

    gEve->Redraw3D(kFALSE, kTRUE);
    fStatusBar_->SetText(Form("Loaded event %d.", event_id_), 1);
    gSystem->ProcessEvents();
  }

  void EventDisplay::update_html_summary()
  {
    // Update summary of current event.

    TEveElement::List_i i;
    TEveElement::List_i j;
    Int_t k;
    TEveElement *el;
    DelphesHtmlObjTable *table;
    TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
    if(mgr)
    {
      htmlSummary_->Clear("D");
      for(i = mgr->BeginChildren(); i != mgr->EndChildren(); ++i)
      {
        el = ((TEveElement *)(*i));
        if(el->IsA() == TEvePointSet::Class())
        {
          TEvePointSet *ps = (TEvePointSet *)el;
          TString ename = ps->GetElementName();
          TString etitle = ps->GetElementTitle();
          if(ename.First('\'') != kNPOS)
            ename.Remove(ename.First('\''));
          etitle.Remove(0, 2);
          Int_t nel = atoi(etitle.Data());
          table = htmlSummary_->AddTable(ename, 0, nel);
        }
        else if(el->IsA() == TEveTrackList::Class())
        {
          TEveTrackList *tracks = (TEveTrackList *)el;
          TString ename = tracks->GetElementName();
          if(ename.First('\'') != kNPOS)
            ename.Remove(ename.First('\''));
          table = htmlSummary_->AddTable(ename.Data(), 5,
            tracks->NumChildren(), kTRUE, "first");
          table->SetLabel(0, "Momentum");
          table->SetLabel(1, "P_t");
          table->SetLabel(2, "Phi");
          table->SetLabel(3, "Theta");
          table->SetLabel(4, "Eta");
          k = 0;
          for(j = tracks->BeginChildren(); j != tracks->EndChildren(); ++j)
          {
            Float_t p = ((TEveTrack *)(*j))->GetMomentum().Mag();
            table->SetValue(0, k, p);
            Float_t pt = ((TEveTrack *)(*j))->GetMomentum().Perp();
            table->SetValue(1, k, pt);
            Float_t phi = ((TEveTrack *)(*j))->GetMomentum().Phi();
            table->SetValue(2, k, phi);
            Float_t theta = ((TEveTrack *)(*j))->GetMomentum().Theta();
            table->SetValue(3, k, theta);
            Float_t eta = theta > 0.0005 && theta < 3.1413 ? ((TEveTrack *)(*j))->GetMomentum().Eta() : 1e10;
            table->SetValue(4, k, eta);
            ++k;
          }
        }
      }
      htmlSummary_->Build();
      gHtml_->Clear();
      gHtml_->ParseText((char *)htmlSummary_->Html().Data());
      gHtml_->Layout();
    }
  }

  /******************************************************************************/
  // GUI
  /******************************************************************************/

  void EventDisplay::make_gui()
  {
    // Create minimal GUI for event navigation.

    // add a tab on the left
    TEveBrowser *browser = gEve->GetBrowser();
    browser->SetWindowName("HEPData4ML Event Display");
    browser->StartEmbedding(TRootBrowser::kLeft);

    // set the main title
    Int_t xDim = 1000;
    Int_t yDim = 600;
    TGMainFrame *frmMain = new TGMainFrame(gClient->GetRoot(), xDim, yDim);
    frmMain->SetWindowName("HEPData4ML Event Display");
    frmMain->SetCleanup(kDeepCleanup);

    // build the navigation menu
    TString icondir;
    if(gSystem->Getenv("ROOTSYS"))
      icondir = Form("%s/icons/", gSystem->Getenv("ROOTSYS"));
    if(!gSystem->OpenDirectory(icondir))
      icondir = Form("%s/icons/", (const char *)gSystem->GetFromPipe("root-config --etcdir"));
    TGGroupFrame *vf = new TGGroupFrame(frmMain, "Event navigation", kVerticalFrame | kFitWidth);
    {
      TGHorizontalFrame *hf = new TGHorizontalFrame(frmMain);
      {
        TGPictureButton *b = 0;

        b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoBack.gif"));
        hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
        b->Connect("Clicked()", "EventDisplay", this, "Bck()");

        // TODO: Hack
        Int_t nentries = 1;

        TGNumberEntry *numberEntry = new TGNumberEntry(hf, 0, 9, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, nentries);
        hf->AddFrame(numberEntry, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 2, 0, 10, 10));
        this->Connect("EventChanged(Int_t)", "TGNumberEntry", numberEntry, "SetIntNumber(Long_t)");
        numberEntry->GetNumberEntry()->Connect("TextChanged(char*)", "EventDisplay", this, "PreSetEv(char*)");
        numberEntry->GetNumberEntry()->Connect("ReturnPressed()", "EventDisplay", this, "GoTo()");

        b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoForward.gif"));
        hf->AddFrame(b, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 2, 10, 10, 10));
        b->Connect("Clicked()", "EventDisplay", this, "Fwd()");
      }
      vf->AddFrame(hf, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

      TGHProgressBar *progress = new TGHProgressBar(frmMain, TGProgressBar::kFancy, 100);
      // TODO: Hack
      Int_t nentries = 1;
      progress->SetMax(nentries);
      progress->ShowPosition(kTRUE, kFALSE, "Event %.0f");
      progress->SetBarColor("green");
      vf->AddFrame(progress, new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 5));
      this->Connect("EventChanged(Int_t)", "TGHProgressBar", progress, "SetPosition(Float_t)");
    }
    frmMain->AddFrame(vf, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
    vf = new TGGroupFrame(frmMain, "Batch operations", kVerticalFrame | kFitWidth);
    {
      TGTextButton *b = new TGTextButton(vf, "Initialize Summary Plots");
      vf->AddFrame(b, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX, 10, 10, 10, 10));
      b->Connect("Clicked()", "EventDisplay", this, "InitSummaryPlots()");
    }
    frmMain->AddFrame(vf, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

    frmMain->MapSubwindows();
    frmMain->Resize();
    frmMain->MapWindow();
    browser->StopEmbedding();
    browser->SetTabTitle("Event Control", 0);

    // the summary tab
    htmlSummary_ = new DelphesHtmlSummary("HEPData4ML Event Display Summary Table");
    TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    gHtml_ = new TGHtml(0, 100, 100);
    TEveWindowFrame *wf = slot->MakeFrame(gHtml_);
    gHtml_->MapSubwindows();
    wf->SetElementName("Summary tables");

    // plot tab
    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    TEveWindowTab *tab = slot->MakeTab();
    tab->SetElementName("Summary plots");
    tab->SetShowTitleBar(kFALSE);
    // plotSummary_ = new DelphesPlotSummary(tab);
    // plotSummary_->Init(elements_);
    // plotSummary_->Connect("Progress(Int_t)", "EventDisplay", this, "DisplayProgress(Int_t)");
    fStatusBar_->SetText("Ready.", 1);
  }

  void EventDisplay::Fwd()
  {
    return;
    // if(event_id_ < treeReader_->GetEntries() - 2)
    // {
    //   EventChanged(event_id_ + 1);
    // }
    // else
    // {
    //   printf("Already at last event.\n");
    // }
  }

  void EventDisplay::Bck()
  {
    return;
    // if(event_id_ > 0)
    // {
    //   EventChanged(event_id_ - 1);
    // }
    // else
    // {
    //   printf("Already at first event.\n");
    // }
  }

  void EventDisplay::PreSetEv(char *ev)
  {
    event_id_tmp_ = Int_t(atoi(ev));
  }

  void EventDisplay::GoTo()
  {
    return;
    // if(event_id_tmp_ >= 0 && event_id_tmp_ < treeReader_->GetEntries() - 1)
    // {
    //   EventChanged(event_id_tmp_);
    // }
    // else
    // {
    //   printf("Error: no such event.\n");
    // }
  }

  void EventDisplay::InitSummaryPlots()
  {
    return;
    // plotSummary_->FillSample(treeReader_, event_id_);
    // plotSummary_->FillEvent();
    // plotSummary_->Draw();
  }

  void EventDisplay::DisplayProgress(Int_t p)
  {
    fStatusBar_->SetText(Form("Processing... %d %%", p), 1);
    gSystem->ProcessEvents();
  }

}