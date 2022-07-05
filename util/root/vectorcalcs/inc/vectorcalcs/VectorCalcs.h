#ifndef ROOT_VECTORCALCS
#define ROOT_VECTORCALCS

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#endif

// standard lib includes
#include <vector>

using namespace std;
namespace VectorCalcs{
  vector<Double_t> PtEtaPhiM2PxPyPzEsingle(vector<Double_t> input_coords, ROOT::Math::PtEtaPhiMVector* vec);
  vector<Double_t> PtEtaPhiM2EPxPyPzsingle(vector<Double_t> input_coords, ROOT::Math::PtEtaPhiMVector* vec);
  vector<Double_t> PxPyPzE2PtEtaPhiMsingle(vector<Double_t> input_coords, ROOT::Math::PxPyPzEVector* vec);

  vector<vector<Double_t>> PtEtaPhiM2PxPyPzE(vector<vector<Double_t>> input_vecs);
  vector<vector<Double_t>> PtEtaPhiM2EPxPyPz(vector<vector<Double_t>> input_vecs);
  vector<vector<Double_t>> PxPyPzE2PtEtaPhiM(vector<vector<Double_t>> input_vecs);

  vector<Double_t> PtEtaPhiM2PxPyPzEflat(vector<Double_t> input_vecs);
  vector<Double_t> PtEtaPhiM2EPxPyPzflat(vector<Double_t> input_vecs);
  vector<Double_t> PxPyPzE2PtEtaPhiMflat(vector<Double_t> input_vecs);

  Double_t DeltaR2(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
  vector<Double_t> DeltaR2Vectorized(vector<Double_t> eta1, vector<Double_t> phi1, vector<Double_t> eta2, vector<Double_t> phi2);
}

#endif