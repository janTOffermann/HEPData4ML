/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DisplayTestUtils_h
#define DisplayTestUtils_h

#include <vector>

#include "Rtypes.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

namespace Display{

  // Structure to hold jet four-vector components
  struct JetFourVector {
      Double_t pt;
      Double_t eta;
      Double_t phi;
      Double_t m;

      JetFourVector(Double_t pt_, Double_t eta_, Double_t phi_, Double_t m_)
          : pt(pt_), eta(eta_), phi(phi_), m(m_) {}
  };

  // Structure to hold track components
  struct Track {
      Double_t pt;
      Double_t eta;
      Double_t phi;
      Double_t m;
      Double_t x;
      Double_t y;
      Double_t z;
      Int_t pdgId;

      Track(Double_t pt_, Double_t eta_, Double_t phi_, Double_t m_, Double_t x_, Double_t y_, Double_t z_, Int_t pdgId_)
          : pt(pt_), eta(eta_), phi(phi_), m(m_), x(x_), y(y_), z(z_), pdgId(pdgId_) {}
  };

  std::vector<JetFourVector> generateRandomJets(Int_t nJets = 10, UInt_t seed = 0);
  std::vector<ROOT::Math::PtEtaPhiMVector> generateRandomJetVecs(Int_t nJets = 10, UInt_t seed = 0);

  std::vector<Track> generateRandomTracks(Int_t nTracks = 20, UInt_t seed = 0, std::vector<std::pair<Int_t, Double_t>> particleTypes={});
  std::vector<Track> generateRandomElectrons(Int_t nTracks = 20, UInt_t seed = 0);
  std::vector<Track> generateRandomMuons(Int_t nTracks = 20, UInt_t seed = 0);
  template<typename Container, typename Func>
  std::vector<Double_t> extractValues(const Container& objects, Func extractor){
      std::vector<Double_t> result;
      result.reserve(objects.size());
      for (const auto& obj : objects) {
          result.emplace_back(extractor(obj));
      }
      return result;
  };

  template<typename Container, typename Func>
  std::vector<Int_t> extractValuesInt(const Container& objects, Func extractor){
      std::vector<Int_t> result;
      result.reserve(objects.size());
      for (const auto& obj : objects) {
          result.emplace_back(extractor(obj));
      }
      return result;
  }

}

#endif /* DisplayTestUtils_h */
