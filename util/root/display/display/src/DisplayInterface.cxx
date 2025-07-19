#include <display/DisplayInterface.h>
#include <display/EventDisplay.h>

// standard library includes
#include <iostream>

// ROOT includes
#include "TSystem.h"
#include "TGeoManager.h"

// Delphes includes
#include "classes/DelphesClasses.h"
#include "display/DelphesEventDisplay.h"
#include "display/Delphes3DGeometry.h"
#include "display/TestUtils.h"

using namespace std;

namespace Display{

  DisplayInterface::DisplayInterface(TString delphesCardFilepath){

    delphesCardFilepath_ = delphesCardFilepath;

    // load the libraries
    gSystem->Load("libGeom");
    gSystem->Load("libGuiHtml");
    gSystem->Load("libDelphesDisplay");

    Delphes3DGeometry det3D_geom(new TGeoManager("delphes", "Delphes geometry"), kTRUE);
    det3D_geom.readFile(delphesCardFilepath_, "ParticlePropagator", "ChargedHadronTrackingEfficiency", "MuonEfficiency", "HCal");

    eventDisplay_ = new EventDisplay(det3D_geom);

  };

  DisplayInterface::~DisplayInterface(){};

  void DisplayInterface::DisplayEvent(TString delphesCardFilepath){

    if(!delphesCardFilepath.EqualTo(delphesCardFilepath_)){
      delete eventDisplay_;
      Delphes3DGeometry det3D_geom(new TGeoManager("delphes", "Delphes geometry"), kTRUE);
      det3D_geom.readFile(delphesCardFilepath_, "ParticlePropagator", "ChargedHadronTrackingEfficiency", "MuonEfficiency", "HCal");
      eventDisplay_ = new EventDisplay(det3D_geom);
    }

    // TestData();

    eventDisplay_->Display();

    // // display
    // // det3D_geom.getDetector()->Draw();
    // det3D_geom.getDetector()->Draw("ogl");
  }

  void DisplayInterface::Test(TString filename){
    cout << "DisplayInterface::Test()" << endl;
    TGeoManager::Import(filename);
   // gGeoManager->DefaultColors();
    gGeoManager->SetMaxVisNodes(5000);
    // gGeoManager->SetVisLevel(4);
    gGeoManager->GetVolume("ATLS")->Draw("ogl");
    // new TBrowser;
  }

  void DisplayInterface::TestData(){

    cout << "Generating some test data." << endl;
    vector<ROOT::Math::PtEtaPhiMVector> jets =  generateRandomJetVecs(3);
    TString jetName = "TestJets";

    cout << "Adding some test data." << endl;

    eventDisplay_->AddJetData_PtEtaPhiM(
      jetName,
      extractValues(jets, [](const auto& j) { return j.Pt(); }),
      extractValues(jets, [](const auto& j) { return j.Eta(); }),
      extractValues(jets, [](const auto& j) { return j.Phi(); }),
      extractValues(jets, [](const auto& j) { return j.M(); }),
      0.4
    );

    vector<Track> tracks = generateRandomTracks(5);

    eventDisplay_->AddTrackData("TestTracks",
      extractValues(tracks, [](const auto& t) { return t.pt;}),
      extractValues(tracks, [](const auto& t) { return t.eta;}),
      extractValues(tracks, [](const auto& t) { return t.phi;}),
      extractValues(tracks, [](const auto& t) { return t.m;}),
      extractValues(tracks, [](const auto& t) { return t.x;}),
      extractValues(tracks, [](const auto& t) { return t.y;}),
      extractValues(tracks, [](const auto& t) { return t.z;}),
      extractValuesInt(tracks, [](const auto& t) { return t.pdgId;})
    );

    vector<Track> electrons = generateRandomElectrons(2);

    eventDisplay_->AddElectronData("TestElectrons",
      extractValues(electrons, [](const auto& t) { return t.pt;}),
      extractValues(electrons, [](const auto& t) { return t.eta;}),
      extractValues(electrons, [](const auto& t) { return t.phi;}),
      extractValues(electrons, [](const auto& t) { return t.m;}),
      extractValues(electrons, [](const auto& t) { return t.x;}),
      extractValues(electrons, [](const auto& t) { return t.y;}),
      extractValues(electrons, [](const auto& t) { return t.z;}),
      extractValues(electrons, [](const auto& t) { return t.pdgId > 0 ? -1. : 1.;})
    );

    vector<Track> muons = generateRandomElectrons(1);
    eventDisplay_->AddMuonData("TestMuons",
      extractValues(muons, [](const auto& t) { return t.pt;}),
      extractValues(muons, [](const auto& t) { return t.eta;}),
      extractValues(muons, [](const auto& t) { return t.phi;}),
      extractValues(muons, [](const auto& t) { return t.m;}),
      extractValues(muons, [](const auto& t) { return t.x;}),
      extractValues(muons, [](const auto& t) { return t.y;}),
      extractValues(muons, [](const auto& t) { return t.z;}),
      extractValues(muons, [](const auto& t) { return t.pdgId > 0 ? -1. : 1.;})
    );

    cout << "... Added." << endl;
  }

}