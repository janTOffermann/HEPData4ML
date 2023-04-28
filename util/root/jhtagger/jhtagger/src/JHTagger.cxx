#include <jhtagger/JHTagger.h>

// Fastjet includes
#include "fastjet/Selector.hh"

using namespace std;

namespace JHTagger{

  JohnnyTagger::JohnnyTagger(Double_t delta_p, Double_t delta_r, Double_t cos_theta_w_max){
    SetDeltaP(delta_p);
    SetDeltaR(delta_r);
    SetCosThetaWMax(cos_theta_w_max);
    _vec = new fastjet::PseudoJet(0.,0.,0.,0.);
    return;
  }

  JohnnyTagger::JohnnyTagger(Double_t delta_p, Double_t delta_r, Double_t cos_theta_w_max,Double_t top_mass_min, Double_t top_mass_max, Double_t W_mass_min, Double_t W_mass_max){
    SetDeltaP(delta_p);
    SetDeltaR(delta_r);
    SetCosThetaWMax(cos_theta_w_max);
    SetTopMassRange(top_mass_min,top_mass_max);
    SetWMassRange(W_mass_min,W_mass_max);
    CreateTagger(); // with all the settings explicitly defined, we may as well create the tagger within this constructor
    return;
  }

  JohnnyTagger::~JohnnyTagger(){
    delete _vec;
    delete _tagger;
    delete _jetdef;
  }

  void JohnnyTagger::InitializeCamAachAlgo(){
    _jetdef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,_R);
    return;
  }

  void JohnnyTagger::CreateTagger(){
    _tagger = new fastjet::JHTopTagger(_delta_p,_delta_r,_cos_theta_w_max);
    _tagger->set_top_selector(fastjet::SelectorMassRange(_top_mass_min,_top_mass_max));
    _tagger->set_W_selector(fastjet::SelectorMassRange(_W_mass_min,_W_mass_max));

    // Also create the jet clusterer.
    InitializeCamAachAlgo();
  }

  void JohnnyTagger::TagJet(fastjet::PseudoJet jet){
    fastjet::PseudoJet top_candidate = _tagger->operator()(jet);
    delete _vec; // TODO: Is this necessary? My C++ is rusty but I wonder if constantly assigning "_vec = new ..." I have introduced a memory leak.
    _vec = 0;
    if (top_candidate != 0){
      _status = kTRUE;
      _vec = new fastjet::PseudoJet(top_candidate.structure_of<fastjet::JHTopTagger>().W());
      GetWCandidateConstituents();
    }
    else{
      _status = kFALSE;
      _vec = new fastjet::PseudoJet(0.,0.,0.,0.);
      ResetWCandidateConstituents();
    }
  }

  void JohnnyTagger::TagJet(std::vector<Double_t> E, std::vector<Double_t> px, std::vector<Double_t> py, std::vector<Double_t> pz){
    Size_t n = E.size();

    std::vector<fastjet::PseudoJet> particles;
    for(Size_t i = 0; i < n; i++){
      particles.push_back( fastjet::PseudoJet( px.at(i), py.at(i), pz.at(i), E.at(i) ) );
    }

    fastjet::ClusterSequence cs(particles, *_jetdef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    if(jets.size() == 0) return;

    fastjet::PseudoJet jet = jets.at(0);
    return TagJet(jet);
  }

  std::vector<Double_t> JohnnyTagger::GetWCandidateConstituentsProperty(TString property){
    std::vector<Double_t> values = {};
    for(fastjet::PseudoJet vec : _vec_constituents){
      if(property.EqualTo("E")) values.push_back(vec.E());
      else if(property.EqualTo("px")) values.push_back(vec.px());
      else if(property.EqualTo("py")) values.push_back(vec.py());
      else if(property.EqualTo("pz")) values.push_back(vec.pz());
      else if(property.EqualTo("pt")) values.push_back(vec.pt());
      else if(property.EqualTo("eta")) values.push_back(vec.eta());
      else if(property.EqualTo("phi")) values.push_back(vec.phi());
      else if(property.EqualTo("m")) values.push_back(vec.m());
      else if(property.EqualTo("y")) values.push_back(vec.rapidity());
      else{
        std::cout << Form("Error: Property \"%s\" passed to JohnnyTagger::GetWCandidateConstituentsProperty not understood",property.Data()) << std::endl;
        break;
      }

    }
    return values;
  }


}