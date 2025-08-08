#ifndef JetElementList_h
#define JetElementList_h

#include "TEveElement.h"
#include "TEveScene.h"

namespace Display{
  class JetCone;
  class JetElementList : public TEveElementList{
  public:
    JetElementList(const char* name = "Jets", const char* title = "") : TEveElementList(name, title) {}

    virtual Bool_t SetRnrChildren(Bool_t rnr) override;

    // virtual Bool_t SetRnrSelf(Bool_t rnr) override;

  private:
      TEveScene* fCaloScene = nullptr;

      void SyncShapes(Bool_t visible);

  public:
      void SetCaloScene(TEveScene* scene) { fCaloScene = scene; }
  };
}

#endif /* JetElementList_h */