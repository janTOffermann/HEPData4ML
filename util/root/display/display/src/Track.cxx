#include "display/Track.h"
#include "TParticlePDG.h"

namespace Display
{

  //______________________________________________________________________________
  Track::Track() :
    TEveTrack(),
    fChargeDouble(0)
  {

  }

  Track::Track(TParticle* t, Int_t label, TEveTrackPropagator* prop) :
    TEveTrack(t, label, prop),
    fChargeDouble(0)
  {

    TParticlePDG* pdgp = t->GetPDG();
    if (pdgp){
        fPdg    = pdgp->PdgCode();
        fCharge = (Int_t) TMath::Nint(pdgp->Charge()/3);
        fChargeDouble = pdgp->Charge()/3.;
    }
  }

  //______________________________________________________________________________

  void Track::SetStdTitle(){
   TString idx(fIndex == kMinInt ? "<undef>" : Form("%d", fIndex));
   TString lbl(fLabel == kMinInt ? "<undef>" : Form("%d", fLabel));
   SetTitle(Form("Index=%s, Label=%s\nChg=%.3f, Pdg=%d\n"
                 "pT=%.3f, pZ=%.3f\nV=(%.3f, %.3f, %.3f)",
                 idx.Data(), lbl.Data(), fChargeDouble, fPdg,
                 fP.Perp(), fP.fZ, fV.fX, fV.fY, fV.fZ));
}

  void Track::SetCharge(Double_t charge){
    fChargeDouble = charge;
  }

  Double_t Track::GetCharge(){
    return fChargeDouble;
  }
}
