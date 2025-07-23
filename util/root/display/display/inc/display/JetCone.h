#ifndef JetCone_h
#define JetCone_h

#include "TEveJetCone.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TEveGeoShape.h"
#include "TEveScene.h"

namespace Display {

class JetCone : public TEveJetCone
{
public:
  JetCone();
  JetCone(const Text_t* n, const Text_t* t);
  virtual ~JetCone();

  TEveGeoShape* GetJetCircle(){return fShape;};
  void SetJetCircle(TEveGeoShape* shape){fShape = shape;};

  virtual void ElementChanged(Bool_t update_scenes = kTRUE, Bool_t redraw = kFALSE);

  // virtual Bool_t SetRnrSelf(Bool_t rnr);

  void SetCaloScene(TEveScene* scene){fCaloScene = scene;};


  ClassDef(JetCone, 1); // ROOT dictionary generation

protected:
  TGeoTube*     fTube;   // ROOT naming convention with 'f' prefix
  TGeoVolume*   fVol;
  TEveGeoShape* fShape;

private:
  TEveScene* fCaloScene;


};

} // namespace Display

#endif /* JetCone_h */