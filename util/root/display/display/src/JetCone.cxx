#include "display/JetCone.h"
#include "TGeoManager.h"
#include "TGeoMedium.h"
#include "TEveManager.h"
#include <iostream>

// ClassImp(Display::JetCone)

namespace Display
{

  //______________________________________________________________________________
  JetCone::JetCone() :
    TEveJetCone(),
    fTube(nullptr),
    fVol(nullptr),
    fShape(nullptr)
  {

  }

  //______________________________________________________________________________
  JetCone::JetCone(const Text_t* n, const Text_t* t) :
    TEveJetCone(n, t),
    fTube(nullptr),
    fVol(nullptr),
    fShape(nullptr){}

  //______________________________________________________________________________
  JetCone::~JetCone()
  {
    if(fShape) fShape = nullptr; // TEve will handle destruction?
  }

  void JetCone::ElementChanged(Bool_t updateScenes, Bool_t redraw){
      TEveJetCone::ElementChanged(updateScenes, redraw);

      if (fShape) {
          // Bool_t shouldRender = GetRnrSelf() && GetRnrState();
          // fShape->SetRnrSelf(shouldRender);
          // fShape->SetRnrChildren(shouldRender);
          fShape->SetRnrSelf(GetRnrSelf());
          fShape->SetRnrChildren(GetRnrChildren());

          if (updateScenes && fCaloScene) {
              fCaloScene->Changed();
          }
      }
      SyncShapeVisibility();
  }

  Bool_t JetCone::SetRnrSelf(Bool_t rnr){
      Bool_t result = TEveJetCone::SetRnrSelf(rnr);
      SyncShapeVisibility();
      return result;
  }

  void JetCone::SyncShapeVisibility(){
      if (fShape) {
          // Use GetRnrState() which considers parent visibility too
          Bool_t effective_visibility = GetRnrState();
          fShape->SetRnrSelf(effective_visibility);

          if (fCaloScene) {
              fCaloScene->Changed();
          }
      }
  }

void JetCone::ParentTEvesChanged()
{
    // TEveJetCone::ParentTEvesChanged();
    SyncShapeVisibility();
}


}