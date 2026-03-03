#ifndef PILEUP
#define PILEUP

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
#include "TH1D.h"
#include "TRandomGen.h"
#include "Math/RotationZ.h"
#endif

// standard lib includes
#include <vector>
#include <utility> // std::pair
#include <tuple> // std::tuple
#include <random> // std::mt19937, std::random_device
#include <unordered_set>

// forward declarations for HepMC3
namespace HepMC3{
  class GenEvent;
  class FourVector;
  class WriterRootTree;
}

using namespace std;
namespace Pileup{

  struct FileRange {
      ULong_t start;
      ULong_t end;
      TString filename;
  };

  // A helper class, for using ROOT's RNG with things like std::shuffle
  class TRandomAdapter{
    public:
      typedef unsigned int result_type;

      TRandomAdapter(TRandom* rng) : fRandom(rng) {}
      static constexpr result_type min() { return 0; }
      static constexpr result_type max() { return UINT_MAX; }

      result_type operator()(){
        return fRandom->Integer(UINT_MAX);
      }

    private:
      TRandom* fRandom;
  };

  class PileupMixer{
    public:
      PileupMixer();
      ~PileupMixer();

      void AddPileupFile(TString filename){_pileupFilenames.push_back(filename);};
      void ClearPileupFiles(){_pileupFilenames.clear();};

      void operator()(TString inputFile, TString outputFile);
      void operator()(ULong_t nEvents, TString outputFile);

      // Setters
      void SetPileupFiles(vector<TString> filenames){_pileupFilenames = filenames;};
      void SetUseContiguousSampling(Bool_t flag){_useContiguousSampling = flag;}

      void SetStableOnly(Bool_t flag){_stableOnly = flag;};
      void SetHTCondorInfo(Bool_t flag, Int_t jobNumber, Int_t nJobs);
      void SetRNGSeed(Int_t seed);
      void SetAllowReuse(Bool_t flag){_allowReuse = flag;};
      void SetUsePhiRotations(Bool_t flag){_allowPhiRotations = flag;};

      void SetBeamSpotSigma(Double_t dt=0.16, Double_t dx=0.01, Double_t dy=0.01, Double_t dz=35.);
      void SetBatchSize(Int_t batchSize){_batchSize = batchSize;};

      Int_t GetRNGSeed(){return _rngSeed;};
      Int_t GetBatchSize(){return _batchSize;};

      // Various ways to initialize the mu distribution
      void InitMuDistribution(Double_t muAvg = 33.7, Double_t muSigma = 11.5);
      void InitMuDistribution(TH1D* muDistributionHistogram);
      TH1D* GetMuDistribution(){return _muDistribution;};

      void Initialize();

    protected:
      void _InitRNG();

      void _InitializeIndexMap();

      Int_t _SampleMuDistribution();
      void _PickEventIndices(Int_t nEvents);
      void _FetchEventSingleFile(const TString& filename, const vector<ULong_t>& localIndices, vector<HepMC3::GenEvent*> &events);
      vector<HepMC3::GenEvent*> _FetchEvents();
      void _AddPileupSingle(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation);
      void _AddPileupSingleStableOnly(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation);

      vector<HepMC3::GenEvent*> _CreatePileupEvents(Int_t nEvents);
      HepMC3::GenEvent* _CreatePileupEvent(const vector<HepMC3::GenEvent*>& pileup, const vector<vector<Double_t>>& pileupDisplacements, const vector<Double_t>& pileupPhiRotations);

      void _Flush(vector<HepMC3::GenEvent*> events);

      // Various helper functions
      vector<Double_t> _generateDisplacement();
      vector<Double_t> _rotateVectorPhi(HepMC3::FourVector* vector, const ROOT::Math::RotationZ& rotation);
      vector<Double_t> _rotateAndTranslateVector(HepMC3::FourVector* v, const vector<Double_t>& displacementCoordinates, const ROOT::Math::RotationZ& rotation);

      Bool_t _initialized = kFALSE; // when initialized, will actually count number of pileup events available
      Int_t _batchSize = 100;
      Bool_t _stableOnly = kTRUE;
      Bool_t _condorFlag = kFALSE;
      Int_t _jobNumber = -1;
      Int_t _nJobs = -1;
      Int_t _rngSeed = 0;
      Bool_t _allowReuse = kTRUE;
      Bool_t _allowPhiRotations = kTRUE;

      // Pileup event-related things
      vector<TString> _pileupFilenames = {};
      ULong_t _nPileupEvents = 0; // total number of pileup events in input
      vector<FileRange> _fileRanges; // sorted by start index
      Bool_t _indexingMapInitialized = kFALSE;
      std::unordered_set<ULong_t> _usedPileupIndices;
      Bool_t _useContiguousSampling = kFALSE;

      // Mu and chosen pileup event indices
      TH1D* _muDistribution = 0;
      Bool_t _muInitialized = kFALSE;
      Double_t _muAvg = 0.; // average mu when doing simple Gaussian initialization
      Double_t _muSigma = 0.;// sigma of mu when doing simple Gaussian initialization
      vector<ULong_t> _selectedGlobalPileupIndices = {};

      // Beam spot parameters
      vector<Double_t> _beamSpotSigma = {0.16, 0.01, 0.01, 35.};

      // Output writer -- making it a member variable so that we can keep the writer open.
      HepMC3::WriterRootTree* _writer = 0;

      //  RNG stuff
      TRandomMixMax17* _rng = 0;
      TRandomAdapter* _rngAdapter = 0;

      // Warning printout stuff
      Int_t _nWarning = 0;
      Int_t _nWarningMax = 5;
    };
}

#endif