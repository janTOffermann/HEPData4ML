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
      MultiSelection(std::vector<BaseSelector*> selections, Bool_t enforceUnique=kFALSE);
      ~MultiSelection(){};

      vector<Int_t> operator()(HepMC3::GenEvent* evt) const;

    protected:
      vector<BaseSelector*> _selections = {};
      Bool_t _enforceUnique = kFALSE;
  };


  // ---------------------------------------
  // Selection Algorithms, for AlgoSelection
  // ---------------------------------------
  // TODO: Is this an OK way to organize things? Can't simply peel off into
  //       another header, because the algorithms use some selector stuff,
  //       and the selectors use some algorithms... (circular!)

  // Helper function
  vector<Int_t> GetDaughtersSingle(HepMC3::GenEvent* evt, Int_t idx);

  class BaseSelectorAlgorithm{
    public:
      BaseSelectorAlgorithm(){}
      virtual ~BaseSelectorAlgorithm() = default;

      Bool_t GetStatus(){return _status;};
      virtual vector<Int_t> operator()(HepMC3::GenEvent* evt) const = 0;

    protected:
      Bool_t _status = kFALSE;

  };

  class SelectDaughters : public BaseSelectorAlgorithm{
    /*
     * This algorithm selects the immediate daughters of
     * the particles selected by _selection.
     * No check is made on these particles's statuses, and
     * this algorithm is not recursive -- it only looks
     * for the immediate daughters.
     */
    public:
      SelectDaughters(BaseSelector* selection);
      ~SelectDaughters(){};

      vector<Int_t> operator()(HepMC3::GenEvent* evt);

    protected:
      unique_ptr<BaseSelector> _selection = 0;
  };

  class SelectStableDaughters : public BaseSelectorAlgorithm{
    /*
     * This algorithm selects the stable daughters of
     * the particles selected by _selection.
     * This will search the entire decay tree, so the
     * results are not necessarily the *immediate*
     * daughters of the results of _selection.
     */
    public:
      SelectStableDaughters(BaseSelector* selection);
      ~SelectStableDaughters(){};

      vector<Int_t> operator()(HepMC3::GenEvent* evt);

    protected:
      unique_ptr<BaseSelector> _selection = 0;
      vector<Int_t> _GetStableDaughters(HepMC3::GenEvent* evt, Int_t idx);
  };

  // Selector using the algorithms.
  // TODO: Can the headers be rearranged so that this is grouped with other selectors?
  class AlgoSelection : public BaseSelector{
    public:
      AlgoSelection(unique_ptr<BaseSelectorAlgorithm> algorithm, Int_t n=-1)
      : _algorithm(std::move(algorithm)), _N(n){};
      ~AlgoSelection(){};

      vector<Int_t> operator()(HepMC3::GenEvent* evt) const;


    protected:
      unique_ptr<BaseSelectorAlgorithm> _algorithm;
      Int_t _N = -1;
  };


}

#endif