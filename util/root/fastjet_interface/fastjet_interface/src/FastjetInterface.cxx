#include <fastjet_interface/FastjetInterface.h>

//standard library includes
#include <algorithm> // std::transform

// fastjet includes
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;

namespace FastjetInterface{

  JetFinder::JetFinder(){};

  JetFinder::~JetFinder() = default;

  void JetFinder::SetRadius(Double_t radius){
    if(radius == _radius) return;
    _radius = radius;
    _initializeAlgorithm();
  }

  void JetFinder::SetJetAlgorithm(TString algorithmName){
    if(_algorithmName.EqualTo(algorithmName)) return;
    _algorithmName = algorithmName;
    _algorithmName.ToLower();
    _initializeAlgorithm();
  }

  void JetFinder::_initializeAlgorithm(){
    if(_jetdef != 0) delete _jetdef;

    if(_algorithmName.EqualTo("kt")) _jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm,_radius);
    else if(_algorithmName.Contains("cambridge") || _algorithmName.Contains("aachen") || _algorithmName.Contains("ca")){
      _jetdef = new fastjet::JetDefinition(fastjet::cambridge_aachen_algorithm,_radius);
    }
    else _jetdef = new fastjet::JetDefinition(fastjet::antikt_algorithm,_radius);
  }

  void JetFinder::SetInputs(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz){
    Size_t n = E.size();
    _input_particles.clear();
    _rapidity.clear();
    _eta.clear();
    _phi.clear();
    _m.clear();
    for(Size_t i = 0; i < n; i++){
      auto particle = make_unique<fastjet::PseudoJet>(px[i], py[i], pz[i], E[i]);
      particle->set_user_index(i);
      _input_particles.push_back(std::move(particle));
    }
  }

  void JetFinder::SetInputs(const vector<Double_t> flattenedVectors){
    Size_t n_elements = flattenedVectors.size();
    if((Int_t)n_elements % 4 != 0){
      std::cerr << "Error: flattenedVectors size must be a multiple of 4" << std::endl;
      return;
    }

    Size_t n_particles = n_elements / 4;

    vector<Double_t> E, px, py, pz;
    E.reserve(n_particles);
    px.reserve(n_particles);
    py.reserve(n_particles);
    pz.reserve(n_particles);

    for(Size_t i = 0; i < n_particles; i++){
      E.push_back(flattenedVectors[4*i]);
      px.push_back(flattenedVectors[4*i + 1]);
      py.push_back(flattenedVectors[4*i + 2]);
      pz.push_back(flattenedVectors[4*i + 3]);
    }

    SetInputs(E, px, py, pz);  // Reuse the other method
  }

  void JetFinder::_etaPhiMToRapidityPhi(vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m){
    //NOTE: Requires inputs to already be set
    vector<Double_t> rapidity;
    Size_t n = _input_particles.size();
    rapidity.reserve(n);
    ROOT::Math::PxPyPzEVector vec;
    for(Size_t i = 0; i < n; i++){
      if(m.at(i) == 0.) rapidity.push_back(eta.at(i));
      else{ // rapidity and pseudo-rapidity are not equal
        const fastjet::PseudoJet* particle = _input_particles[i].get();
        vec.SetCoordinates(particle->px(),particle->py(),particle->pz(),particle->E()); // TODO: Any advantage over calling particle->Rapidity()?
        rapidity.push_back(vec.Rapidity());
      }
    }
    SetRapidityPhi(rapidity,phi);
  }

  void JetFinder::_etaPhiMToRapidityPhi(vector<Double_t> flattenedEtaPhiM){
    Size_t n_elements = flattenedEtaPhiM.size();
    if((Int_t)n_elements % 3 != 0){
      std::cerr << "Error: flattenedEtaPhiM size must be a multiple of 3" << std::endl;
      return;
    }

    Size_t n_particles = n_elements / 3;

    vector<Double_t> eta,phi,m;
    eta.reserve(n_particles);
    phi.reserve(n_particles);
    m.reserve(  n_particles);

    for(Size_t i = 0; i < n_particles; i++){
      eta.push_back(flattenedEtaPhiM[3*i]);
      phi.push_back(flattenedEtaPhiM[3*i + 1]);
      m.push_back(  flattenedEtaPhiM[3*i + 2]);
    }

    _etaPhiMToRapidityPhi(eta,phi,m);  // Reuse the other method
  }

  void JetFinder::SetInputsWithEtaPhiM(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz, vector<Double_t> eta, vector<Double_t> phi, vector<Double_t> m){
    SetInputs(E,px,py,pz);
    _etaPhiMToRapidityPhi(eta,phi,m);
  }

  void JetFinder::SetInputsWithEtaPhiM(const vector<Double_t> flattenedVectors, const vector<Double_t> flattenedEtaPhiM){
    SetInputs(flattenedVectors);
    _etaPhiMToRapidityPhi(flattenedEtaPhiM);
  }

  void JetFinder::ClusterJets(vector<Double_t> E, vector<Double_t> px, vector<Double_t> py, vector<Double_t> pz){
    SetInputs(E,px,py,pz);
    ClusterJets();
  }

  void JetFinder::ClusterJets(){
    vector<fastjet::PseudoJet> particles;
    particles.reserve(_input_particles.size());
    Bool_t hasRapidityPhi = (_rapidity.size() > 0);

    for(Size_t i = 0; i < _input_particles.size(); i++){
      particles.push_back(*_input_particles.at(i));  // Dereference and copy
      if(hasRapidityPhi) particles.back().set_cached_rap_phi(_rapidity.at(i),_phi.at(i));
    }

    fastjet::ClusterSequence cs(particles, *_jetdef);
    vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

    // Now, get these jets' pointers -- a bit roundabout but necessary due to ROOT library loading forcing us
    // to use forward declarations of fastjet in the header (and thus just pointers for member variables!).
    _jets.clear(); // unique_ptr will handle deletion
    _jetConstituents.clear();

    // Store new jets and their constituent indices
    for(Size_t i = 0; i < jets.size(); i++){
      // Store the jet
      _jets[i] = make_unique<fastjet::PseudoJet>(jets[i]);

      // extract constituent indices
      vector<fastjet::PseudoJet> constituents = fastjet::sorted_by_pt(jets[i].constituents());
      vector<Int_t> constituent_indices;

      for(Size_t j = 0; j < constituents.size(); j++){
        constituent_indices.push_back(constituents[j].user_index());
      }

      _jetConstituents[i] = constituent_indices;
    }
  }

  void JetFinder::ClusterJets(const vector<Double_t> flattenedVectors){
    SetInputs(flattenedVectors);
    ClusterJets();
  }

  vector<Int_t> JetFinder::GetJetIndices(){
    vector<Int_t> indices;
    for(const auto& pair : _jets){
      indices.push_back(pair.first);
    }
    return indices;
  }

  map<Int_t, vector<Double_t>> JetFinder::GetJetMomenta(){
    map<Int_t, vector<Double_t>> momenta;
    for(const auto& pair : _jets){
      Int_t idx = pair.first;
      const fastjet::PseudoJet* jet = pair.second.get();

      vector<Double_t> p4 = {jet->E(), jet->px(), jet->py(), jet->pz()};
      momenta[idx] = p4;
    }
    return momenta;
  }

  map<Int_t, vector<Double_t>> JetFinder::GetJetMomentaCylindrical(){
    map<Int_t, vector<Double_t>> momenta;
    for(const auto& pair : _jets){
      Int_t idx = pair.first;
      const fastjet::PseudoJet* jet = pair.second.get();

      vector<Double_t> p4 = {jet->pt(), jet->eta(), jet->phi(), jet->m()};
      momenta[idx] = p4;
    }
    return momenta;
  }

  map<Int_t, vector<vector<Double_t>>> JetFinder::GetJetConstituentMomenta(){
    map<Int_t, vector<vector<Double_t>>> constituent_momenta;

    for(const auto& pair : _jetConstituents){
      Int_t jet_idx = pair.first;
      const vector<Int_t>& indices = pair.second;

      vector<vector<Double_t>> momenta;
      for(Size_t i = 0; i < indices.size(); i++){
        Int_t particle_idx = indices[i];
        const fastjet::PseudoJet* particle = _input_particles[particle_idx].get();

        vector<Double_t> p4 = {particle->E(), particle->px(),
                              particle->py(), particle->pz()};
        momenta.push_back(p4);
      }
      constituent_momenta[jet_idx] = momenta;
    }
    return constituent_momenta;
  }

  map<Int_t, vector<vector<Double_t>>> JetFinder::GetJetConstituentMomentaCylindrical(){
    map<Int_t, vector<vector<Double_t>>> constituent_momenta;

    for(const auto& pair : _jetConstituents){
      Int_t jet_idx = pair.first;
      const vector<Int_t>& indices = pair.second;

      vector<vector<Double_t>> momenta;
      for(Size_t i = 0; i < indices.size(); i++){
        Int_t particle_idx = indices[i];
        const fastjet::PseudoJet* particle = _input_particles[particle_idx].get();

        vector<Double_t> p4 = {particle->pt(), particle->eta(),
                              particle->phi(), particle->m()};
        momenta.push_back(p4);
      }
      constituent_momenta[jet_idx] = momenta;
    }
    return constituent_momenta;
  }

  map<Int_t, vector<Int_t>> JetFinder::GetJetConstituentIndices(){
    return _jetConstituents;
  }








}