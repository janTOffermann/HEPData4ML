#include <vectorcalcs/VectorCalcs.h>
using namespace std;

namespace VectorCalcs{
  Calculator::Calculator(){
    _vec1 = new ROOT::Math::PtEtaPhiMVector();
    _vec2 = new ROOT::Math::PtEtaPhiMVector();
    _vec_xyze = new ROOT::Math::PxPyPzEVector();
  }

  Calculator::~Calculator(){
    delete _vec1;
    delete _vec2;
    delete _vec_xyze;
  }

  vector<Double_t> Calculator::PtEtaPhiM_to_EPxPyPz(vector<Double_t> input_coords){
    return PtEtaPhiM_to_EPxPyPz(input_coords.at(0),input_coords.at(1),input_coords.at(2),input_coords.at(3));
  }

  vector<Double_t> Calculator::PtEtaPhiM_to_EPxPyPz(Double_t pt, Double_t eta, Double_t phi, Double_t m){
    _vec1->SetCoordinates(pt,eta,phi,m);
    vector<Double_t> result {_vec1->E(), _vec1->Px(), _vec1->Py(), _vec1->Pz()};
    return result;
  }

  vector<Double_t> Calculator::PxPyPzE_to_PtEtaPhiM(vector<Double_t> input_coords){
    return PxPyPzE_to_PtEtaPhiM(input_coords.at(0),input_coords.at(1),input_coords.at(2),input_coords.at(3));
  }

  vector<Double_t> Calculator::PxPyPzE_to_PtEtaPhiM(Double_t px, Double_t py, Double_t pz, Double_t e){
    _vec_xyze->SetCoordinates(px,py,pz,e);
    vector<Double_t> result {_vec_xyze->Pt(), _vec_xyze->Eta(), _vec_xyze->Phi(), _vec_xyze->M()};
    return result;
  }

  vector<Double_t> Calculator::PxPyPzE_to_PtEtaPhiM_Multi(vector<Double_t> input_vectors){
    // TODO: This is broken!
    Int_t nvecs = Int_t(input_vectors.size()) / 4;
    vector<Double_t> results(input_vectors.size(),0.);
    vector<Double_t> result_single(4,0.);

    for(Int_t i = 0; i < nvecs; i++){
      result_single = PxPyPzE_to_PtEtaPhiM(input_vectors[4 * i], input_vectors[4 * i + 1], input_vectors[4 * i + 2], input_vectors[4 * i + 3]);
      results[4 * i] = result_single[0];
      results[4 * i + 1] = result_single[1];
      results[4 * i + 2] = result_single[2];
      results[4 * i + 3] = result_single[3];
    }
    return results;
  }

  vector<Double_t> Calculator::PtEtaPhiM_to_EPxPyPz_Multi(vector<Double_t> input_vectors){
    Int_t nvecs = Int_t(input_vectors.size()) / 4;
    vector<Double_t> results(input_vectors.size(),0.);
    vector<Double_t> result_single(4,0.);

    for(Int_t i = 0; i < nvecs; i++){
      result_single = PtEtaPhiM_to_EPxPyPz(input_vectors[4 * i], input_vectors[4 * i + 1], input_vectors[4 * i + 2], input_vectors[4 * i + 3]);
      results[4 * i] = result_single[0];
      results[4 * i + 1] = result_single[1];
      results[4 * i + 2] = result_single[2];
      results[4 * i + 3] = result_single[3];
    }
    return results;
  }

  Double_t Calculator::DeltaR2(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
    _vec1->SetCoordinates(0.,eta1,phi1,0.);
    _vec2->SetCoordinates(0.,eta2,phi2,0.);
    return ROOT::Math::VectorUtil::DeltaR2(*_vec1,*_vec2);
  }

  vector<Double_t> Calculator::DeltaR2Vectorized(vector<Double_t> eta1, vector<Double_t> phi1, vector<Double_t> eta2, vector<Double_t> phi2){
    Int_t n = Int_t(eta1.size());
    Int_t m = Int_t(eta2.size());
    vector<Double_t> result(n * m, 0.);
    for(Int_t i = 0; i < n; i++){
      for(Int_t j = 0; j < m; j++){
        result[m * i + j] = DeltaR2(eta1.at(i),phi1.at(i),eta2.at(j),phi2.at(j));
      }
    }
    return result;
  }

}