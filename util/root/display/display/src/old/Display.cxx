#include <display/Display.h>

// standard library includes
#include <iostream>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TEveArrow.h"
#include "TEveBrowser.h"
#include "TEveCalo.h"
#include "TEveElement.h"
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

// // Delphes includes
// #include "display/DelphesEventDisplay.h"
// #include "display/Delphes3DGeometry.h"
// #include "display/DelphesBranchElement.h"

using namespace std;

namespace EventDisplay{

  EventDisplay::EventDisplay(){
    event_id_ = 0;
    tkRadius_ = 1.29;
    totRadius_ = 2.0;
    tkHalfLength_ = 3.0;
    bz_ = 3.8;
    chain_ = new TChain("Delphes");
    treeReader_ = 0;
    delphesDisplay_ = 0;
  };

  EventDisplay::~EventDisplay(){

    // clean up TClonesArrays
    for (auto it = clone_arrays_.begin(); it != clone_arrays_.end(); it++)  delete it->second;
  };

  void EventDisplay::BuildDetector(Delphes3DGeometry &det3D){
    // initialize the application
    TEveManager::Create(kTRUE, "IV");
    fStatusBar_ = gEve->GetBrowser()->GetStatusBar();
    TGeoManager *geom = gGeoManager;

    // build the detector
    tkRadius_ = det3D.getTrackerRadius();
    totRadius_ = det3D.getDetectorRadius();
    tkHalfLength_ = det3D.getTrackerHalfLength();
    muHalfLength_ = det3D.getDetectorHalfLength();
    bz_ = det3D.getBField();
    etaAxis_ = det3D.getCaloAxes().first;
    phiAxis_ = det3D.getCaloAxes().second;
    TGeoVolume *top = det3D.getDetector(false);
    geom->SetTopVolume(top);
    TEveElementList *geometry = new TEveElementList("Geometry");
    TObjArray *nodes = top->GetNodes();
    TIter itNodes(nodes);
    TGeoNode *nodeobj;
    TEveGeoTopNode *node;
    while((nodeobj = (TGeoNode *)itNodes.Next()))
    {
      node = new TEveGeoTopNode(gGeoManager, nodeobj);
      node->UseNodeTrans();
      geometry->AddElement(node);
    }
    return;
  }

  // void EventDisplay::AddCaloData(TString name, std::vector<BranchBase *> elements){

  //   // Create a TClonesArray for this data.
  //   TString className = "DelphesCaloData";
  //   Int_t expectedSize = 1000;
  //   clone_arrays_[name] = new TClonesArray(className,expectedSize);

  //   // BranchElement<DelphesCaloData> *clist = new BranchElement<DelphesCaloData>(name, clone_arrays_[name], kBlack);
  //   // clist->GetContainer()->SetEtaBins(etaAxis_);
  //   // clist->GetContainer()->SetPhiBins(phiAxis_);
  //   // elements.push_back(clist);
  // }

  //   void EventDisplay::AddTrackData(TString name, std::vector<BranchBase *> elements){

  //   // Create a TClonesArray for this data.
  //   TString className = "TEveTrackList";
  //   Int_t expectedSize = 1000;
  //   clone_arrays_[name] = new TClonesArray(className,expectedSize);

  //   BranchElement<TEveTrackList> *clist = new BranchElement<TEveTrackList>(name, clone_arrays_[name], kBlack);
  //   clist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
  //   elements.push_back(clist);
  // }


}