#ifndef NTUPLE_PRODUCER
#define NTUPLE_PRODUCER

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

// #include <ntuple/ParticleSelections.h>

// forward declarations for HepMC3
namespace HepMC3{
  class GenEvent;
  class ReaderRootTree;
  class ReaderAscii;
}

// forward declarations for our own stuff,
// need this due to how PyROOT loads headers
namespace NtupleProducer {
    class BaseSelector;
}

using namespace std;
namespace NtupleProducer{

  struct FourMomentumData{ // a simple container for holding both (E,px,py,pz) and (pt,eta,phi,m) bases.
    Int_t N;

    /*
      * In principle, we know that the vector will always be length 4.
      * However we will do vector<vector<Double_t>> instead of vector<>
      * of some length-4 array struct since reading nested vectors with
      * ROOT is much easier (for the struct we'll need the class loaded,
      * otherwise it can't be done nicely with TTreeReader -- just see how
      * we have to handle the "Edges" of certain Delphes objects). -Jan
      */
    vector<vector<Double_t>> pmu;
    vector<vector<Double_t>> pmu_cyl;

    void Clear(){
      N = 0;
      pmu.clear();
      pmu_cyl.clear();
    }
  };

  struct TrackData{ // a simple container for holding some (typically) track-related data.
    Int_t N;

    vector<Double_t> d0;
    vector<Double_t> d0Error;
    vector<Double_t> z0;
    vector<Double_t> z0Error;
    vector<vector<Double_t>> Xd;

    void Clear(){
      N = 0;
      d0.clear();
      d0Error.clear();
      z0.clear();
      z0Error.clear();
      Xd.clear();
    }
  };

  struct CaloData{ // a simple container for holding some calorimeter-related data.
    Int_t N;

    vector<Double_t> Eem;
    vector<Double_t> Ehad;
    vector<Double_t> Etrack;

    // As oppose to doing a vector of "TwoVector" structs, we're keeping these
    // as vectors of vectors, so that in principle we can accomodate more complex
    // tower geometries than rectangles in (eta,phi). Of course that would require
    // some changes upstream in Delphes. - Jan
    vector<vector<Double_t>> edgesEta;
    vector<vector<Double_t>> edgesPhi;

    void Clear(){
      N = 0;
      Eem.clear();
      Ehad.clear();
      Etrack.clear();
      edgesEta.clear();
      edgesPhi.clear();
    }
  };

  struct ParticleData {
    Int_t N;
    FourMomentumData momentum;
    vector<Int_t> pdgId;
    vector<Int_t> indexHepMC;
    vector<vector<Double_t>> xmu_prod;

    vector<vector<Double_t>> xmu_decay;
    vector<Bool_t> isStable;

    void Clear(){
      N = 0;
      momentum.Clear();
      pdgId.clear();
      indexHepMC.clear();
      xmu_prod.clear();

      isStable.clear();
      xmu_decay.clear();
    }
  };

  struct DelphesReaderOutput{
    Int_t N; // object multiplicity
    FourMomentumData momentum;
    TrackData trackData;
    CaloData caloData;
    vector<vector<Double_t>> positionData;

    vector<Int_t> charge;
    vector<Int_t> pdgId;

    void Clear() {
      N = 0;
      momentum.Clear();
      trackData.Clear();
      caloData.Clear();
      charge.clear();
      pdgId.clear();
      positionData.clear();
    }
  };

  struct DelphesReaderData {

    // Output data
    DelphesReaderOutput output;

    // Input readers
    std::unique_ptr<TTreeReaderArray<Float_t>> pt;
    std::unique_ptr<TTreeReaderArray<Float_t>> eta;
    std::unique_ptr<TTreeReaderArray<Float_t>> phi;
    std::unique_ptr<TTreeReaderArray<Float_t>> mass;

    std::unique_ptr<TTreeReaderArray<Float_t>> d0;
    std::unique_ptr<TTreeReaderArray<Float_t>> d0Error;
    std::unique_ptr<TTreeReaderArray<Float_t>> z0;
    std::unique_ptr<TTreeReaderArray<Float_t>> z0Error;

    std::unique_ptr<TTreeReaderArray<Float_t>> xd;
    std::unique_ptr<TTreeReaderArray<Float_t>> yd;
    std::unique_ptr<TTreeReaderArray<Float_t>> zd;

    std::unique_ptr<TTreeReaderArray<Int_t>> charge; // NOTE: Not currently supporting fractional charges
    std::unique_ptr<TTreeReaderArray<Int_t>> pdgId;

    std::unique_ptr<TTreeReaderArray<Float_t>> Eem;
    std::unique_ptr<TTreeReaderArray<Float_t>> Ehad;
    std::unique_ptr<TTreeReaderArray<Float_t>> Etrack;

    // Edges of calorimeter cells -- this is a bit tricky because they
    // are fixed-length arrays within a collection.
    // We can't use TTreeReaderArray for this (yet), so we have to
    // handle this kind of branch the old-fashioned way. The extra twist
    // is that there are multiple kinds of Delphes objects that have an
    // Edges array, but they're not all the same length (typically 4, sometimes 2).
    // It'd also be nice to support arbitrary length, in case one has
    // a more complex cell geometry (though this would require changes upstream
    // in Delphes, and hopefully this would include moving to more modern
    // objects like std::vector...)
    Bool_t hasEdges = kFALSE;
    Int_t edgesSize = 4; // default -- can be adjusted if needed
    const static Int_t edgesMax = 400000; // unfortunately we need to set a max size -- so make it large (1.6MB)
    Float_t Edges[edgesMax];

    std::unique_ptr<TTreeReaderArray<Float_t>> T;
    std::unique_ptr<TTreeReaderArray<Float_t>> X;
    std::unique_ptr<TTreeReaderArray<Float_t>> Y;
    std::unique_ptr<TTreeReaderArray<Float_t>> Z;

    // will add other leaves as needed...

    void Clear() {
      output.Clear();
    }

  };

  class Converter{

    public:
      Converter();
      ~Converter();

      void SetInputFilesHepMC(vector<TString> files);
      void AddInputFileHepMC(TString file);

      void SetInputFilesDetector(vector<TString> files);
      void AddInputFileDetector(TString file);

      void SetDelphesObjects(vector<TString> objectNames){_delphesObjectNames = objectNames;};
      void AddDelphesObject(TString objectName){_delphesObjectNames.push_back(objectName);};
      void ResetDelphesObjects(){_delphesObjectNames.clear();};

      void AddTruthParticleSelector(TString selectionName, BaseSelector* selector);

      void Process(TString inputFileHepMC, TString inputFileDetector, TString outputFile);
      // ----
      void SetStableTruthParticleName(TString name){_truthParticleBranchPrefix = name;};
      void SetDelphesDefaultMass(TString branchName, Double_t mass){_delphesMassDefault[branchName] = mass;};
      void SetTreeName(TString name){_treeName = name;};
      TString GetTreeName(){return _treeName;};

    private:

      void _OpenHepMC3File(TString filename);
      void _OpenHepMC3FileRoot(TString filename);
      void _OpenHepMC3FileAscii(TString filename);
      Bool_t _ReadHepMCEvent();
      Bool_t _failedHepMC();

      void _OpenDelphesFile(TString filename);

      void _SetOutputFiles();
      void _SetOutputFiles(vector<TString> files){_outputFiles = files;};

      void _CreateHepMCBranches();
      void _CreateHepMC3BranchesSingle(TString particleCollectionName, ParticleData& data, Bool_t extra=kFALSE); // convenience func, used by _CreateHepMCBranches()

      void _CreateDelphesBranches();
      void _CreateDelphesBranch(TString inputBranchName); // for making a single branch
      void _DelphesMultiplicity(TString inputBranchName);
      void _DelphesMomentum(TString inputBranchName, vector<TString> attributes);
      void _DelphesD0Z0(TString inputBranchName, vector<TString> attributes);
      void _DelphesXd(TString inputBranchName, vector<TString> attributes);
      void _DelphesPdgIdCharge(TString inputBranchName, vector<TString> attributes);
      void _DelphesCalorimeter(TString inputBranchName, vector<TString> attributes);
      void _DelphesPosition(TString inputBranchName, vector<TString> attributes);

      void _FillStableParticles(); // fills the *stable* truth particles, which we always do
      void _FillTruthParticles(); // fills any user-specified truth particles
      void _FillDelphesObjects();
      void _IterateDelphesTree(Int_t entry);

      // Utility funcs
      Bool_t _CheckStringVector(vector<TString> v, TString target);

      // Input and output filenames -- can set these beforehand, then process them all.
      vector<TString> _inputFilesHepMC = {};
      vector<TString> _inputFilesDetector = {}; // e.g. Delphes files
      vector<TString> _outputFiles = {};
      Bool_t _hasDetectorFiles = kFALSE;

      // readers for HepMC3, and associated variables
      HepMC3::ReaderRootTree* _readerRoot = 0;
      HepMC3::ReaderAscii* _readerAscii = 0;
      Bool_t _rootMode = kTRUE; // kTRUE for ROOT, kFALSE for ASCII
      HepMC3::GenEvent* _evt = 0;
      // buffers for filling
      ParticleData _stableParticles;
      map<TString, ParticleData> _truthParticleStructs = {};

      // HepMC3 event record particle selectors
      map<TString, unique_ptr<BaseSelector>> _truthParticleSelectors = {};

      // reader for Delphes, and associated variables
      TTreeReader* _delphesReader = 0;
      TFile* _delphesFile = 0;
      TTree* _delphesTree = 0;
      vector<TString> _delphesObjectNames = {}; // which Delphes objects to copy over
      vector<TString> _delphesLeafNames = {};
      // buffers and variables related to filling output
      map<TString, Double_t> _delphesMassDefault = {};
      map<TString, Bool_t> _delphesFillXd = {};
      map<TString, Bool_t> _delphesAddedN = {};
      map<TString, Bool_t> _delphesIsTrack = {}; // keep track of what objects are "track-like"
      map<TString, std::unique_ptr<DelphesReaderData>> _delphesData = {};

      // variables associated with output ntuple
      TFile* _outputFile = 0;
      TTree* _outputTree = 0;
      TString _treeName = "hepdata4ml_tree";

      // event index -- for looping on the HepMC3/Delphes files, since that's done in a "while" loop (due to how HepMC3 works)
      ULong64_t _i = 0;

      // misc
      TString _truthParticleBranchPrefix = "StableTruthParticles";
  };
}

#endif