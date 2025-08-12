#ifndef Track_h
#define Track_h

#include "TEveTrack.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TEveGeoShape.h"
#include "TEveScene.h"

namespace Display {

/* Child class of TEveTrack, that (partially) supports non-integer charges (in units of |e|).
 * It will display these fractional charges in the mouseover tooltip in the event display,
 * however the integer charge is still used to actually propagate the track (changing this sadly
 * requires some deeper changes to TEveTrackPropagator).
 */
class Track : public TEveTrack
{
public:
  Track();
  Track(TParticle* t, Int_t label, TEveTrackPropagator* prop=nullptr);

  virtual void SetCharge(Double_t charge);
  virtual Double_t GetCharge();
  virtual void SetStdTitle();

  // virtual void MakeTrack(Bool_t recurse=kTRUE); // NOTE: To actually support fractional charges, need to modify/inherit TEveTrackPropagator, and its Helix_t...

  ClassDef(Track, 1); // ROOT dictionary generation

protected:

  Double_t fChargeDouble; // could consider making it fCharge?

  // TODO: As a quick hack, consider setting effective charges and masses such that the track propagates as if its charge is non-integer?
  //       e.g. normalize charge to magnitude of 1, adjust mass accordingly.

};

} // namespace Display

#endif /* Track_h */