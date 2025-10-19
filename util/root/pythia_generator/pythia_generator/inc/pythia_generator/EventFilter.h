#ifndef EVENT_FILTER
#define EVENT_FILTER

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
// #include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif

// standard lib includes
#include <vector>
#include <utility> // std::pair
#include <tuple> // std::tuple

// // forward declarations for Pythia8
namespace Pythia8{
  class Pythia;
  class PythiaParallel;
  class Event;
  class Particle;
}


using namespace std;
namespace PythiaGenerator{

  class EventFilter{
    public:
      EventFilter(){};
      virtual ~EventFilter() = default;
      virtual Bool_t operator()(Pythia8::Pythia* pythiaPtr) const = 0;
  };

  class ParticleFilter : public EventFilter{
    public:
      ParticleFilter(Int_t pdgId, Int_t statusHepMC, Double_t ptMin=0., Double_t ptMax=-1., Double_t etaMin=-999., Double_t etaMax=999.);
      ~ParticleFilter(){};
      Bool_t operator()(Pythia8::Pythia* pythiaPtr) const;

    protected:
      Int_t _pdgId;
      Int_t _statusHepMC;
      Double_t _ptMin; // GeV
      Double_t _ptMax; // GeV
      Double_t _etaMin; // GeV
      Double_t _etaMax; // GeV
  };


}

#endif