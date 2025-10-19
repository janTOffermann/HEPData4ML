#include <pythia_generator/EventFilter.h>

//standard library includes
#include <algorithm> // std::transform

// Pythia8 includes
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Pythia8/PythiaParallel.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenRunInfo.h"

using namespace std;

namespace PythiaGenerator{

  ParticleFilter::ParticleFilter(Int_t pdgId, Int_t statusHepMC, Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax){
    _pdgId = pdgId;
    _statusHepMC = statusHepMC;
    _ptMin = ptMin;
    _ptMax = ptMax;
    _etaMin = etaMin;
    _etaMax = etaMax;
  }

  Bool_t ParticleFilter::operator()(Pythia8::Pythia* pythiaPtr) const {
    Bool_t result = kFALSE;
    for (Int_t i = 1; i < pythiaPtr->event.size(); i++){ // skip entry 0, which isn't a particle
      if(!(pythiaPtr->event[i].id() == _pdgId)) continue;
      if(_statusHepMC >= 0 && !(pythiaPtr->event[i].statusHepMC() == _statusHepMC)) continue;
      if(_ptMin >= 0. && !(pythiaPtr->event[i].pT() > _ptMin)) continue;
      if(_ptMax > 0. && !(pythiaPtr->event[i].pT() < _ptMax)) continue;
      if(!(pythiaPtr->event[i].eta() > _etaMin)) continue;
      if(!(pythiaPtr->event[i].eta() < _etaMax)) continue;
      result = kTRUE;
      break;
    }
    return result;
  }


}