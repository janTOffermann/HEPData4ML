#ifndef FASTJET_INTERFACE
#define FASTJET_INTERFACE

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
#include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif

// standard lib includes
#include <vector>
#include <utility> // std::pair
#include <tuple> // std::tuple
#include <map>
#include <memory>

// // forward declarations
namespace fastjet{
  class PseudoJet;
  class JetDefinition;
}

using namespace std;
namespace FastjetInterface{

  class JetFinder{
    public:
    JetFinder();
    ~JetFinder();

    ClassDef(JetFinder, 1);

    void SetJetAlgorithm(TString algorithmName);
    void SetRadius(Double_t radius=0.4);

    void SetInputs(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz);
    void SetInputs(const vector<Double_t> flattenedVectors);

    // Alternative input method -- ways of inputting/precomputing rapidity and phi, to speed up fastjet
    void SetInputsWithEtaPhiM(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m);
    void SetInputsWithEtaPhiM(const vector<Double_t> flattenedVectors, const vector<Double_t> flattenedEtaPhiM);

    void SetRapidityPhi(vector<Double_t> rapidity, vector<Double_t> phi){_rapidity=rapidity;_phi=phi;}; // directly set rapidity and phi for speeding up clustering

    void ClusterJets();
    void ClusterJets(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz);
    void ClusterJets(const vector<Double_t> flattenedVectors); // for flattened input

    // Various getters
    vector<Int_t> GetJetIndices(); // get the keys of _jets
    map<Int_t, vector<Double_t>> GetJetMomenta(); // format of each momentum will be (E, px, py, pz)
    map<Int_t, vector<Double_t>> GetJetMomentaCylindrical(); // format of each momentum will be (pt, eta, phi, m)
    map<Int_t, vector<vector<Double_t>>> GetJetConstituentMomenta(); // return vector of 4-momenta, each like (E, px, py, pz) -- possibly unused
    map<Int_t, vector<vector<Double_t>>> GetJetConstituentMomentaCylindrical(); // return vector of 4-momenta, each like (pt, eta, phi, m) -- possibly unused
    map<Int_t, vector<Int_t>> GetJetConstituentIndices();

    private:

    void _initializeAlgorithm();
    void _etaPhiMToRapidityPhi(vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m);
    void _etaPhiMToRapidityPhi(vector<Double_t> flattenedEtaPhiM);

    Double_t _radius;
    TString _algorithmName = "antikt";
    fastjet::JetDefinition* _jetdef = 0;

    /*
    * Containers for inputs; we use this paradigm instead of directly
    * passing arguments to ClusterJets(), because we might sometimes
    * supply additional info (rapidity/eta/phi) that can speed things up.
    */
    vector<unique_ptr<fastjet::PseudoJet>> _input_particles = {};
    vector<Double_t> _rapidity = {};
    vector<Double_t> _eta = {};
    vector<Double_t> _phi = {};
    vector<Double_t> _m = {};

    /*
     * We'll handle jets using a map, where the keys are their indices.
     * This should make it easy to uniquely identify jets, even if we
     * throw some out from the map.
     */
    map<Int_t, unique_ptr<fastjet::PseudoJet>> _jets = {};
    map<Int_t, vector<Int_t>> _jetConstituents = {}; // keep track of jet constituents' indices w.r.t. initial jet clustering input
  };
}

#endif