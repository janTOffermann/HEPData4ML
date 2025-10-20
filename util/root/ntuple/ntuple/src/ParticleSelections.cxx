#include <ntuple/ParticleSelections.h>

//standard library includes
#include <algorithm> // std::transform. std::find

// ROOT includes
#include "TMath.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

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
      if (_status != 0 && p->status() != _status) return false;
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

  MultiSelection::MultiSelection(std::vector<BaseSelector*> selections, Bool_t enforceUnique)
  : _enforceUnique(enforceUnique) {
    for (auto selector : selections) {
        _selections.push_back(selector);
    }
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

  AlgoSelection::AlgoSelection(BaseSelectorAlgorithm* algorithm, Int_t n)
      : _N(n){
      _algorithm = unique_ptr<BaseSelectorAlgorithm>(algorithm);
  }

  vector<Int_t> AlgoSelection::operator()(HepMC3::GenEvent* evt) const {
    vector<Int_t> result = (*_algorithm)(evt);
    if(_N > 0 && result.size() > _N){ // truncate
      result.resize(_N);
    }
    return result;
  }

  // ---------------------------------------
  // Selection Algorithms, for AlgoSelection
  // ---------------------------------------

  vector<Int_t> GetDaughtersSingle(HepMC3::GenEvent* evt, Int_t idx){
    HepMC3::GenParticlePtr startingParticle = evt->particles().at(idx);
    //NOTE: For HepmC3::GenParticle::children(),
    // "Less efficient than via the vertex since return
    // must be by value (in case there is no vertex)".
    Int_t status = startingParticle->status();
    if(status == 1) return {};
    HepMC3::GenVertexPtr vertex = startingParticle->end_vertex();
    vector<HepMC3::GenParticlePtr> children = vertex->particles_out();
    vector<Int_t> indices = {};
    for(auto const p : children){
      // GenParticle.id() gives the HepMC3 index, which is *not* the index
      // of the particle w.r.t. particles(). It is, in fact, always 1 more,
      // because the particles are 1-indexed in the event listing.
      // So it should be safe to use id() - 1.
      indices.push_back(p->id() - 1);
    }
    return indices;
  }

  SelectDaughters::SelectDaughters(BaseSelector* selection){
    _selection = unique_ptr<BaseSelector>(selection);
  }

  vector<Int_t> SelectDaughters::operator()(HepMC3::GenEvent* evt) const{
    _status = kTRUE;
    vector<Int_t> startingParticleIndices = (*_selection)(evt);

    vector<Int_t> allDaughterIndices = {};

    for (Int_t startingIndex : startingParticleIndices){
      auto daughters = GetDaughtersSingle(evt, startingIndex);
      allDaughterIndices.insert(allDaughterIndices.end(), daughters.begin(), daughters.end());
    }

    if(allDaughterIndices.size() == 0){
      _status = kFALSE;
    }
    else{
      // enforce uniqueness
      std::sort(allDaughterIndices.begin(), allDaughterIndices.end());
      allDaughterIndices.erase(std::unique(allDaughterIndices.begin(), allDaughterIndices.end()), allDaughterIndices.end());
    }
    return allDaughterIndices;
  }

  SelectStableDaughters::SelectStableDaughters(BaseSelector* selection){
    _selection = unique_ptr<BaseSelector>(selection);
  }

  vector<Int_t> SelectStableDaughters::_GetStableDaughters(HepMC3::GenEvent* evt, Int_t idx) const{
    vector<Int_t> result = {};
    vector<Int_t> daughters = GetDaughtersSingle(evt, idx);
    for(Int_t j : daughters){
      Int_t status = evt->particles().at(j)->status();
      if(status == 1) result.push_back(j);
      else{
        vector<Int_t> grandDaughters = _GetStableDaughters(evt, j);
        result.insert(result.end(),grandDaughters.begin(),grandDaughters.end());
      }
    }
    return result;
  }

  vector<Int_t> SelectStableDaughters::operator()(HepMC3::GenEvent* evt) const{
    _status = kTRUE;
    vector<Int_t> startingParticleIndices = (*_selection)(evt);

    vector<Int_t> stableDaughterIndices = {};

    for (Int_t startingIndex : startingParticleIndices){
      auto daughters = _GetStableDaughters(evt,startingIndex);
      stableDaughterIndices.insert(stableDaughterIndices.end(), daughters.begin(), daughters.end());
    }

    if(stableDaughterIndices.size() == 0){
      _status = kFALSE;
    }
    else{
      // enforce uniqueness
      std::sort(stableDaughterIndices.begin(), stableDaughterIndices.end());
      stableDaughterIndices.erase(std::unique(stableDaughterIndices.begin(), stableDaughterIndices.end()), stableDaughterIndices.end());
    }
    return stableDaughterIndices;
  }
}