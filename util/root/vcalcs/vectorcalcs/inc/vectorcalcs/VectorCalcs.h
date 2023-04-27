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

  class Calculator{
    public:
      Calculator();
      ~Calculator();

      // Functions for simple coordinate conversions of a single 4-vector.
      // Two possible signatures -- either inputting a vector, or inputting all 4 coordinates.
      vector<Double_t> PxPyPzE_to_PtEtaPhiM(vector<Double_t> input_coords);
      vector<Double_t> PxPyPzE_to_PtEtaPhiM(Double_t px, Double_t py, Double_t pz, Double_t e);

      vector<Double_t> PtEtaPhiM_to_EPxPyPz(vector<Double_t> input_coords);
      vector<Double_t> PtEtaPhiM_to_EPxPyPz(Double_t pt, Double_t eta, Double_t phi, Double_t m);

      // Functions for converting a whole (flattened) list of 4-vectors,
      // like {pt_0, eta_0, phi_0, m_0, pt_1, eta_1, phi_1, m_1, pt_2 ...}.
      vector<Double_t> PxPyPzE_to_PtEtaPhiM_Multi(vector<Double_t> input_vectors);
      vector<Double_t> PtEtaPhiM_to_EPxPyPz_Multi(vector<Double_t> input_vectors);

      // Functions for calculating distances.
      Double_t DeltaR2(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
      // Find the dR^2 between all vecs in list 1 and all vecs in list 2. Returns a flattened list: [d_0_0, d_0_1, d_0_2... d_n_(m-1), d_n_m]
      vector<Double_t> DeltaR2Vectorized(vector<Double_t> eta1, vector<Double_t> phi1, vector<Double_t> eta2, vector<Double_t> phi2);

    private:
      ROOT::Math::PtEtaPhiMVector* _vec1 = 0;
      ROOT::Math::PtEtaPhiMVector* _vec2 = 0;
      ROOT::Math::PxPyPzEVector* _vec_xyze = 0;
  };
}

#endif