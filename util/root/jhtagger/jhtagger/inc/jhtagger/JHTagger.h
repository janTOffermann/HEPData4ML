#ifndef ROOT_FJUTILS
#define ROOT_FJUTILS

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
// #include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif

// Fastjet includes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "fastjet/PseudoJet.hh"

// standard lib includes
#include <vector>

using namespace std;
namespace JHTagger{

  class JohnnyTagger{
    public:
      JohnnyTagger(){};
      JohnnyTagger(Double_t delta_p, Double_t delta_r, Double_t cos_theta_w_max);
      JohnnyTagger(Double_t delta_p, Double_t delta_r, Double_t cos_theta_w_max, Double_t top_mass_min, Double_t top_mass_max, Double_t W_mass_min, Double_t W_mass_max);
      ~JohnnyTagger();

      // Setters.
      void SetJetDeltaR(Double_t value){_R = value;};
      void SetDeltaP(Double_t value){_delta_p = value;};
      void SetDeltaR(Double_t value){_delta_r = value;};
      void SetCosThetaWMax(Double_t value){_cos_theta_w_max = value;};
      void SetTopMassMin(Double_t value){_top_mass_min = value;};
      void SetTopMassMax(Double_t value){_top_mass_max = value;};
      void SetWMassMin(Double_t value){_W_mass_min = value;};
      void SetWMassMax(Double_t value){_W_mass_max = value;};

      void SetTopMassRange(Double_t low, Double_t high){SetTopMassMin(low);SetTopMassMax(high);};
      void SetWMassRange(Double_t low, Double_t high){SetWMassMin(low);SetWMassMax(high);};

      // Getters.
      Double_t GetJetDeltaR(){return _R;};
      Double_t GetDeltaP(){return _delta_p;};
      Double_t GetDeltaR(){return _delta_r;};
      Double_t GetCosThetaWMax(){return _cos_theta_w_max;};
      Double_t GetTopMassMin(){return _top_mass_min;};
      Double_t GetTopMassMax(){return _top_mass_max;};
      Double_t GetWMassMin(){return _W_mass_min;};
      Double_t GetWMassMax(){return _W_mass_max;};
      Bool_t GetStatus(){return _status;};
      fastjet::PseudoJet* GetWCandidate(){return _vec;};

      // Some getters specific to fetching information on the W candidate's constituents.
      Int_t GetWCandidateNConstituents(){return _vec_constituents.size();};
      std::vector<Double_t> GetWCandidateConstituentsProperty(TString property="E");

      void InitializeCamAachAlgo();
      void CreateTagger(); // Construct the internal JH tagger. The various tagger settings need to have already been set!

      void TagJet(fastjet::PseudoJet jet); // Tag an input jet (which should have been clustered using the Cambridge-Aachen algorithm.)
      void TagJet(std::vector<Double_t> E, std::vector<Double_t> px, std::vector<Double_t> py, std::vector<Double_t> pz);

    private:

      // Get a vector of the W candidate's constituents.
      // We'll want other simpler methods for extracting the candidate info
      // via the PyROOT interface, since using this on the Python side seems to cause issues.
      void GetWCandidateConstituents(){_vec_constituents = _vec->constituents();};
      void ResetWCandidateConstituents(){_vec_constituents.clear();};

      Double_t _delta_p = 0.1;
      Double_t _delta_r = 0.19;
      Double_t _cos_theta_w_max = 0.7;
      Bool_t _status = kFALSE;
      fastjet::JHTopTagger* _tagger = 0;
      fastjet::PseudoJet* _vec = 0;
      std::vector<fastjet::PseudoJet> _vec_constituents = {};

      Double_t _R = 0.8;
      fastjet::JetDefinition* _jetdef = 0;

      // Settings for the mass selectors.
      Double_t _top_mass_min = 150.; // GeV
      Double_t _top_mass_max = 200.; // GeV
      Double_t _W_mass_min = 65.; // GeV
      Double_t _W_mass_max = 95.; // GeV



  };
}

#endif