#include <vectorcalcs/VectorCalcs.h>
using namespace std;

// -- Functions for converting a single vector --
vector<Double_t> VectorCalcs::PtEtaPhiM2PxPyPzEsingle(vector<Double_t> input_coords){
  ROOT::Math::PtEtaPhiMVector vec (0.,0.,0.,0.);
  vec.SetCoordinates(input_coords.at(0),input_coords.at(1),input_coords.at(2),input_coords.at(3));
  vector<Double_t> result{vec.Px(),vec.Py(),vec.Pz(),vec.E()};
  return result;
}

vector<Double_t> VectorCalcs::PtEtaPhiM2EPxPyPzsingle(vector<Double_t> input_coords){
  ROOT::Math::PtEtaPhiMVector vec (0.,0.,0.,0.);
  vec.SetCoordinates(input_coords.at(0),input_coords.at(1),input_coords.at(2),input_coords.at(3));
  vector<Double_t> result{vec.E(),vec.Px(),vec.Py(),vec.Pz()};
  return result;
}

vector<Double_t> VectorCalcs::PxPyPzE2PtEtaPhiMsingle(vector<Double_t> input_coords){
  ROOT::Math::PxPyPzEVector vec (0.,0.,0.,0.);
  vec.SetCoordinates(input_coords.at(0),input_coords.at(1),input_coords.at(2),input_coords.at(3));
  vector<Double_t> result{vec.Pt(),vec.Eta(),vec.Phi(),vec.M()};
  return result;
}

// -- Functions for converting multiple vectors, returned as vector of vectors (a bit cumbersome) --

vector<vector<Double_t>> VectorCalcs::PtEtaPhiM2PxPyPzE(vector<vector<Double_t>> input_vecs){
  vector<vector<Double_t>> result;
  for(auto vec : input_vecs) result.push_back(PtEtaPhiM2PxPyPzEsingle(vec));
  return result;
}

vector<vector<Double_t>> VectorCalcs::PtEtaPhiM2EPxPyPz(vector<vector<Double_t>> input_vecs){
  vector<vector<Double_t>> result;
  for(auto vec : input_vecs) result.push_back(PtEtaPhiM2EPxPyPzsingle(vec));
  return result;
}

vector<vector<Double_t>> VectorCalcs::PxPyPzE2PtEtaPhiM(vector<vector<Double_t>> input_vecs){
  vector<vector<Double_t>> result;
  for(auto vec : input_vecs) result.push_back(PxPyPzE2PtEtaPhiMsingle(vec));
  return result;
}

// -- Functions for converting multiple vectors, returns (flattened) 1D vector.
vector<Double_t> VectorCalcs::PtEtaPhiM2PxPyPzEflat(vector<Double_t> input_vecs){
  vector<Double_t> result;
  Int_t nvecs = Int_t(input_vecs.size()) / 4;
  for(Int_t i = 0; i < nvecs; i++){
    vector<Double_t> input_vec{input_vecs.at(4*i),input_vecs.at(4*i + 1),input_vecs.at(4*i + 2),input_vecs.at(4*i + 3)};
    for(Double_t val : PtEtaPhiM2PxPyPzEsingle(input_vec)) result.push_back(val);
  }
  return result;
}

vector<Double_t> VectorCalcs::PtEtaPhiM2EPxPyPzflat(vector<Double_t> input_vecs){
  vector<Double_t> result;
  Int_t nvecs = Int_t(input_vecs.size()) / 4;
  for(Int_t i = 0; i < nvecs; i++){
    vector<Double_t> input_vec{input_vecs.at(4*i),input_vecs.at(4*i + 1),input_vecs.at(4*i + 2),input_vecs.at(4*i + 3)};
    for(Double_t val : PtEtaPhiM2EPxPyPzsingle(input_vec)) result.push_back(val);
  }
  return result;
}

vector<Double_t> VectorCalcs::PxPyPzE2PtEtaPhiMflat(vector<Double_t> input_vecs){
  vector<Double_t> result;
  Int_t nvecs = Int_t(input_vecs.size()) / 4;
  for(Int_t i = 0; i < nvecs; i++){
    vector<Double_t> input_vec{input_vecs.at(4*i),input_vecs.at(4*i + 1),input_vecs.at(4*i + 2),input_vecs.at(4*i + 3)};
    for(Double_t val : PxPyPzE2PtEtaPhiMsingle(input_vec)) result.push_back(val);
  }
  return result;
}

// -- Functions for dR --
Double_t VectorCalcs::DeltaR2(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  ROOT::Math::PtEtaPhiMVector vec1;
  ROOT::Math::PtEtaPhiMVector vec2;
  vec1.SetCoordinates(0.,eta1,phi1,0.);
  vec2.SetCoordinates(0.,eta2,phi2,0.);
  return ROOT::Math::VectorUtil::DeltaR2(vec1,vec2);
}

// Find the dR^2 between all vecs in list 1 and all vecs in list 2. Returns a flattened list: [d_0_0, d_0_1, d_0_2... d_n_(m-1), d_n_m]
vector<Double_t> VectorCalcs::DeltaR2Vectorized(vector<Double_t> eta1, vector<Double_t> phi1, vector<Double_t> eta2, vector<Double_t> phi2){
  Int_t n = Int_t(eta1.size());
  Int_t m = Int_t(eta2.size());
  vector<Double_t> result;

  for(Int_t i = 0; i < n; i++){
    for(Int_t j = 0; j < m; j++){
      result.push_back(DeltaR2(eta1.at(i),phi1.at(i),eta2.at(j),phi2.at(j)));
    }
  }
  return result;
}