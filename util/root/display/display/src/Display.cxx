#include <display/Display.h>

// standard library includes
#include <iostream>

// ROOT includes
#include "TSystem.h"
#include "TGeoManager.h"

// Delphes includes
// #include "classes/DelphesClasses.h"
#include "display/DelphesEventDisplay.h"
#include "display/Delphes3DGeometry.h"

using namespace std;

namespace EventDisplay{

  Display::Display(){
    TString s = "Created Display class.";
    cout << s << endl;
  };

  Display::~Display(){};


  void Display::DisplayEvent(TString delphesCardFilepath){

    cout << "Display::DisplayEvent()" << endl;

    // load the libraries
    gSystem->Load("libGeom");
    gSystem->Load("libGuiHtml");
    gSystem->Load("libDelphesDisplay");


    Delphes3DGeometry det3D_geom(new TGeoManager("delphes", "Delphes geometry"), kFALSE);
    det3D_geom.readFile(delphesCardFilepath, "ParticlePropagator", "ChargedHadronTrackingEfficiency", "MuonEfficiency", "HCal");



  }

}