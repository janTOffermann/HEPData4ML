#ifndef NTUPLE_PRODUCER
#define NTUPLE_PRODUCER

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

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
  class ReaderRootTree;
  class ReaderAscii;
}

using namespace std;
namespace NtupleProducer{

  struct FourVector {
    Double_t data[4];

    // Convenience accessors
    Double_t& operator[](size_t i) { return data[i]; }
    const Double_t& operator[](size_t i) const { return data[i]; }
  };

  struct FourMomentumData{ // a simple container for holding both (E,px,py,pz) and (pt,eta,phi,m) bases.
    Int_t N;

    /*
      * Using vector<FourVector> instead of <vector<vector<Double_t>>
      * will be a bit more awkward for uproot access in Python,
      * but this should offer better compression since the inner
      * "vector" is fixed-length.
      */
    vector<FourVector> pmu;
    vector<FourVector> pmu_cyl;

    void Clear(){
      N = 0;
      pmu.clear();
      pmu_cyl.clear();
    }
  };

  struct TrackData{ // a simple container for holding both (E,px,py,pz) and (pt,eta,phi,m) bases.
    Int_t N;

    /*
      * Using vector<FourVector> instead of <vector<vector<Double_t>>
      * will be a bit more awkward for uproot access in Python,
      * but this should offer better compression since the inner
      * "vector" is fixed-length.
      */
    vector<Double_t> d0;
    vector<Double_t> d0Error;
    vector<Double_t> z0;
    vector<Double_t> z0Error;

    void Clear(){
      N = 0;
      d0.clear();
      d0Error.clear();
      z0.clear();
      z0Error.clear();
    }
  };

  struct TruthParticleData {
    Int_t N;
    FourMomentumData momentum;
    vector<Int_t> pdgId;
    vector<Int_t> indexHepMC;
    vector<FourVector> xmu_prod;

    void Clear(){
      N = 0;
      momentum.Clear();
      pdgId.clear();
      indexHepMC.clear();
      xmu_prod.clear();
    }
  };

  struct DelphesReaderData {

    // Output data
    Int_t N;
    FourMomentumData outputMomentum;
    TrackData trackData;

    // add more as needed...

    // Input readers
    std::unique_ptr<TTreeReaderArray<Float_t>> pt;
    std::unique_ptr<TTreeReaderArray<Float_t>> eta;
    std::unique_ptr<TTreeReaderArray<Float_t>> phi;
    std::unique_ptr<TTreeReaderArray<Float_t>> mass;

    std::unique_ptr<TTreeReaderArray<Float_t>> d0;
    std::unique_ptr<TTreeReaderArray<Float_t>> d0Error;
    std::unique_ptr<TTreeReaderArray<Float_t>> z0;
    std::unique_ptr<TTreeReaderArray<Float_t>> z0Error;

    // will add other leaves as needed...

    void Clear() {
      N = 0;
      outputMomentum.Clear();
      trackData.Clear();
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

      void Process(TString inputFileHepMC, TString inputFileDetector, TString outputFile);


      // ----
      void SetStableTruthParticleName(TString name){_truthParticleBranchPrefix = name;};
      void SetDelphesDefaultMass(TString branchName, Double_t mass){_delphesMassDefault[branchName] = mass;};


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

      void _CreateDelphesBranches();
      void _CreateDelphesBranch(TString inputBranchName); // for making a single branch

      void _FillStableParticles();
      void _FillDelphesObjects();

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
      TruthParticleData _stableParticles;
      vector<TruthParticleData> _truthParticleStructs = {};

      // reader for Delphes, and associated variables
      TTreeReader* _readerDelphes = 0;
      TFile* _fileDelphes = 0;
      TTree* _treeDelphes = 0;
      vector<TString> _delphesObjectNames = {}; // which Delphes objects to copy over
      vector<TString> _delphesLeafNames = {};
      // buffers for filling
      map<TString, Double_t> _delphesMassDefault = {};
      map<TString, std::unique_ptr<DelphesReaderData>> _delphesData = {};

      // variables associated with output ntuple
      TFile* _outputFile = 0;
      TTree* _outputTree = 0;

      // event index -- for looping on the HepMC3/Delphes files, since that's done in a "while" loop (due to how HepMC3 works)
      ULong64_t _i = 0;

      // misc
      TString _truthParticleBranchPrefix = "StableTruthParticles";


  };
}

#endif