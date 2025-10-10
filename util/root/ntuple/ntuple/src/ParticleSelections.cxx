#include <ntuple/ParticleSelections.h>

//standard library includes
#include <algorithm> // std::transform. std::find

// ROOT includes
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TMath.h"
// #include "TLeafF.h"
#include "TLeafElement.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/FourVector.h"

using namespace std;

namespace NtupleProducer{

  FirstSelector::FirstSelector(Int_t status, Int_t pdgId, Bool_t hadronization){
    _status = status;
    _pdgId = pdgId;
    SetHadronization(hadronization);
  }

  void FirstSelector::SetHadronization(Bool_t hadronization){
    _hadronization = hadronization;
    if(!_hadronization && TMath::Abs(_pdgId) >= 1 && TMath::Abs(_pdgId) <= 5){
      _status = 1;
    }
  }

  void FirstSelector::Print(){
    cout << Form("FirstSelector: status = %i, pdgid = %i",_status,_pdgId) << endl;
  }

  vector<Int_t> FirstSelector::operator()(HepMC3::GenEvent* evt) const{
    // Find the first instance of a particle with _status and _pdgId.
    // Doing this in a sort of fancy way! - Jan
    vector<shared_ptr<HepMC3::GenParticle>> particles = evt->particles();

    // Build the predicate based on selection criteria -- fancy!
    auto matches = [this](const shared_ptr<HepMC3::GenParticle>& p) {
      if (p->pid() != _pdgId) return false;
      if (_status != -1 && p->status() != _status) return false;
      return true;
    };

    // Find the first matching particle
    auto it = find_if(particles.begin(), particles.end(), matches);

    if (it == particles.end()) {
      _selectionStatus = kFALSE;
      return {};
    }

    _selectionStatus = kTRUE;
    return {static_cast<Int_t>(distance(particles.begin(), it))};
  }

  vector<Int_t> MultiSelection::operator()(HepMC3::GenEvent* evt) const {
    vector<Int_t> result = {};

    for (const auto& selector : _selections) {
      auto indices = (*selector)(evt);
      result.insert(result.end(), indices.begin(), indices.end());
      if(!selector->GetSelectionStatus()) return {}; // failure
    }

    if (_enforceUnique) {
      std::sort(result.begin(), result.end());
      result.erase(std::unique(result.begin(), result.end()), result.end());
    }

    return result;
  }
}