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

#include "display/TestUtils.h"
#include <TRandom3.h>
#include <TMath.h>
#include<TDatabasePDG.h>

namespace Display{
//------------------------------------------------------------------------------

  // Function to generate random jets with realistic HEP distributions
  std::vector<JetFourVector> generateRandomJets(Int_t nJets, UInt_t seed) {
    TRandom3 rng(seed);  // Use TRandom3 for better quality random numbers
    std::vector<JetFourVector> jets;
    jets.reserve(nJets);

    for (Int_t i = 0; i < nJets; ++i) {
      // Generate pt with a realistic falling spectrum (power law-ish)
      // Most jets at low pt, fewer at high pt
      Double_t pt_min = 20.0;   // GeV, typical jet pt threshold
      Double_t pt_max = 500.0;  // GeV, reasonable upper limit
      Double_t alpha = 3.0;     // Power law exponent

      Double_t r = rng.Uniform(0, 1);
      Double_t pt = pt_min * TMath::Power((pt_max/pt_min), r) * TMath::Power(r, -1.0/alpha);

      // Generate eta uniformly in detector acceptance
      Double_t eta = rng.Uniform(-2.5, 2.5);  // Typical detector coverage

      // Generate phi uniformly over full range
      Double_t phi = rng.Uniform(-TMath::Pi(), TMath::Pi());

      // Generate mass - most jets are light, but allow for some heavy jets
      Double_t mass;
      if (rng.Uniform(0, 1) < 0.1) {
        // 10% chance of "heavy" jets (W/Z/top decay products)
        mass = rng.Uniform(60.0, 120.0);  // GeV
      } else {
        // Light jets (mostly gluon/light quark jets)
        mass = rng.Uniform(0.0, 20.0);    // GeV
      }

      jets.emplace_back(pt, eta, phi, mass);
    }
    return jets;
  }

  // Convenience function that returns TLorentzVector objects directly
  std::vector<ROOT::Math::PtEtaPhiMVector> generateRandomJetVecs(Int_t nJets, UInt_t seed) {
    auto jetData = generateRandomJets(nJets, seed);
    std::vector<ROOT::Math::PtEtaPhiMVector> jets;
    jets.reserve(nJets);

    for (const auto& jet : jetData) {
      ROOT::Math::PtEtaPhiMVector p4;
      p4.SetCoordinates(jet.pt, jet.eta, jet.phi, jet.m);
      jets.push_back(p4);
    }
    return jets;
  }

  std::vector<TrackVector> generateRandomTracks(Int_t nTracks, UInt_t seed, std::vector<std::pair<Int_t, Double_t>> particleTypes) {
    TRandom3 rng(seed);
    std::vector<TrackVector> tracks;
    tracks.reserve(nTracks);

    // Particle types and their relative probabilities
    if(particleTypes.size() == 0){
      particleTypes = {
        {211,  0.40},  // π± (40%)
        {321,  0.10},  // K± (10%)
        {2212, 0.08},  // proton (8%)
        {11,   0.10},  // electron (10%)
        {13,   0.08}   // muon (8%)
      };
    }

    // Get PDG database instance
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    for (Int_t i = 0; i < nTracks; ++i) {
      // Generate pt with realistic track spectrum (softer than jets)
      Double_t pt_min = 0.1;    // GeV, typical tracking threshold
      Double_t pt_max = 50.0;   // GeV, reasonable upper limit for isolated tracks
      Double_t alpha = 2.0;     // Softer spectrum than jets

      Double_t r = rng.Uniform(0, 1);
      Double_t pt = pt_min * TMath::Power((pt_max/pt_min), r) * TMath::Power(r, -1.0/alpha);

      // Generate eta uniformly in tracking acceptance
      Double_t eta = rng.Uniform(-2.5, 2.5);

      // Generate phi uniformly
      Double_t phi = rng.Uniform(-TMath::Pi(), TMath::Pi());

      // Select particle type based on probabilities
      Double_t particleRand = rng.Uniform(0, 1);
      Double_t cumProb = 0.0;
      Int_t basePdgId = 211; // fallback to pion

      for (const auto& particle : particleTypes) {
        cumProb += particle.second;
        if (particleRand <= cumProb) {
          basePdgId = particle.first;
          break;
        }
      }

      // Randomly assign charge (50/50 for charged particles)
      Int_t pdgId = basePdgId;
      if (abs(basePdgId) == 211 || abs(basePdgId) == 321 || abs(basePdgId) == 11 || abs(basePdgId) == 13) {
        if (rng.Uniform(0, 1) < 0.5) {
          pdgId = -abs(basePdgId);
        }
      }
      // For protons, make antiprotons less likely (more realistic)
      else if (abs(basePdgId) == 2212) {
        if (rng.Uniform(0, 1) < 0.33) {  // 1/3 chance for antiproton
          pdgId = -2212;
        }
      }

      // Get mass from PDG database
      TParticlePDG* particle = pdgDB->GetParticle(pdgId);
      Double_t mass = particle ? particle->Mass() : 0.13957; // fallback to pion mass

      // Generate innermost hit position (beam spot + detector geometry)
      Double_t x, y, z;

      // Realistic beam spot and primary vertex distribution
      // Beam spot: Gaussian in x,y with ~20 μm width, longer in z
      Double_t beamSpotX = rng.Gaus(0.0, 0.02);  // 20 μm RMS in x
      Double_t beamSpotY = rng.Gaus(0.0, 0.02);  // 20 μm RMS in y
      Double_t beamSpotZ = rng.Gaus(0.0, 50.0);  // 50 mm RMS in z

      // Inner detector radius (first layer of pixel/silicon detector)
      // Typical values: 30-50 mm for pixel detector inner radius
      Double_t innerRadius = 35.0;  // mm

      // Most tracks from primary vertex, some from secondary vertices
      if (rng.Uniform(0, 1) < 0.85) {
        // Primary vertex tracks - start near beam line
        x = beamSpotX + rng.Gaus(0.0, 0.1);  // Small additional smearing
        y = beamSpotY + rng.Gaus(0.0, 0.1);
        z = beamSpotZ + rng.Gaus(0.0, 10.0);
      }
      else{
        // Secondary vertex tracks (displaced vertices, decays)
        Double_t displaceR = rng.Uniform(1.0, 20.0);  // mm
        Double_t displaceAngle = rng.Uniform(0, 2*TMath::Pi());
        x = beamSpotX + displaceR * TMath::Cos(displaceAngle);
        y = beamSpotY + displaceR * TMath::Sin(displaceAngle);
        z = beamSpotZ + rng.Gaus(0.0, 30.0);
      }

      // Ensure track actually reaches the inner detector
      Double_t trackR = TMath::Sqrt(x*x + y*y);
      if (trackR < innerRadius){
        // Extrapolate to inner detector radius
        Double_t trackPhi = TMath::ATan2(y, x);
        x = innerRadius * TMath::Cos(trackPhi);
        y = innerRadius * TMath::Sin(trackPhi);
        // Adjust z based on track angle
        Double_t theta = 2.0 * TMath::ATan(TMath::Exp(-eta));
        Double_t deltaR = innerRadius - trackR;
        z += deltaR / TMath::Tan(theta);
      }

      tracks.emplace_back(pt, eta, phi, mass, x, y, z, pdgId);
    }
    return tracks;
  }

  std::vector<TrackVector> generateRandomElectrons(Int_t nTracks, UInt_t seed) {
    return generateRandomTracks(nTracks, seed, {{11,   1.0}});
  }

  std::vector<TrackVector> generateRandomMuons(Int_t nTracks, UInt_t seed){
    return generateRandomTracks(nTracks, seed, {{13,   1.0}});
  };
}
//------------------------------------------------------------------------------
