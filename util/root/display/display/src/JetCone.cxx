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

  void JetCone::ElementChanged(Bool_t update_scenes, Bool_t redraw){
      TEveJetCone::ElementChanged(update_scenes, redraw);

      if (fShape) {
          fShape->SetRnrSelf(GetRnrSelf());
          fShape->SetRnrChildren(GetRnrChildren());

          if (update_scenes && fCaloScene) {
              fCaloScene->Changed();
          }
      }
  }

}