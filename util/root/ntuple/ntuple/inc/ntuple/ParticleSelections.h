#ifndef PARTICLE_SELECTIONS
#define PARTICLE_SELECTIONS

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLeaf.h"

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#endif

// standard lib includes
#include <vector>
#include <utility> // std::pair
#include <tuple> // std::tuple
#include <memory> // shared_ptr, unique_ptr

// forward declarations for HepMC3
namespace HepMC3{
  class GenEvent;
}

using namespace std;
namespace NtupleProducer{

  class BaseSelector{
    public:
      BaseSelector(){};
      virtual ~BaseSelector() = default;
      virtual vector<Int_t> operator()(HepMC3::GenEvent* evt) const = 0;

      void SetN(Int_t n){_n = n;};
      Int_t GetN(){return _n;};

      Bool_t GetSelectionStatus(){return _selectionStatus;};
      Bool_t IsFixedLength(){return _fixedLength;};

    protected:
      mutable Bool_t _selectionStatus = kTRUE;
      Bool_t _fixedLength = kTRUE;
      Int_t _n = 1;
  };

  class FirstSelector : public BaseSelector{
    public:
      FirstSelector(Int_t status, Int_t pdgId, Bool_t hadronization=kTRUE);
      ~FirstSelector(){};

      void SetStatus(Int_t status){_status = status;};
      void SetPdgId(Int_t pdgId){_pdgId  = pdgId;};
      void SetHadronization(Bool_t hadronization);
      void Print();

      vector<Int_t> operator()(HepMC3::GenEvent* evt) const;

    protected:
      Int_t _status;
      Int_t _pdgId;
      Bool_t _hadronization;
  };

  // Note the use of "___Selector" for classes that return an Int_t,
  // versus "___Selection" for classes that return vector<Int_t>
  class MultiSelection : public BaseSelector{
    public:
      MultiSelection(std::vector<unique_ptr<BaseSelector>>&& selections, Bool_t enforceUnique=kFALSE)
      : _selections(std::move(selections)), _enforceUnique(enforceUnique) {};
      ~MultiSelection(){};

      vector<Int_t> operator()(HepMC3::GenEvent* evt) const;

    protected:
      vector<unique_ptr<BaseSelector>> _selections = {};
      Bool_t _enforceUnique = kFALSE;
  };

}

#endif