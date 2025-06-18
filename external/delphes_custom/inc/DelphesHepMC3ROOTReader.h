/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2021  Universite catholique de Louvain (UCL), Belgium
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

#ifndef DelphesHepMC3ROOTReader_h
#define DelphesHepMC3ROOTReader_h

/** \class DelphesHepMC3ROOTReader
 *
 *  Reads HepMC3 file, in ROOT format
 *
 *  \author Jan T. Offermann - Brown University
 * Based on the DelphesHepMC3Reader class, written by P. Demin (UCL)
 *
 */

#include <map>
#include <vector>

#include <stdio.h>

// ROOT includes
// #include "RTypes.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"
#include "TString.h"

// HepMC3 includes

// swapped in favor of ReaderRootTree -- still ROOT file, slightly different format.
// In principle, the file type used by ReaderRootTree supports random access, although
// it is not clear to me that this is directly supported by the HepMC ROOT IO classes. - Jan
// #include "HepMC3/ReaderRoot.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/GenEvent.h"

// DELPHES includes
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

class DelphesHepMC3ROOTReader
{
public:
  DelphesHepMC3ROOTReader();
  ~DelphesHepMC3ROOTReader();

  void SetInputFile(TString inputFile);

  void Clear();
  bool EventReady();

  Bool_t ReadEvent();

  bool ReadBlock(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  void AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
    TStopwatch *readStopWatch, TStopwatch *procStopWatch);

  void AnalyzeWeight(ExRootTreeBranch *branch);

  void Analyze(DelphesFactory *factory,
    TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);


private:

  void AnalyzeVertex(DelphesFactory *factory, int code, Candidate *candidate = 0);

  // void AnalyzeParticle(DelphesFactory *factory);
  void AnalyzeParticle(DelphesFactory *factory, std::shared_ptr<HepMC3::GenParticle> particle);


  void FinalizeParticles(TObjArray *allParticleOutputArray,
    TObjArray *stableParticleOutputArray,
    TObjArray *partonOutputArray);

  TString fInputFile;

  char *fBuffer;

  TDatabasePDG *fPDG;

  int fEventNumber, fMPI, fProcessID, fSignalCode, fVertexCounter, fParticleCounter;
  double fScale, fAlphaQCD, fAlphaQED;

  double fMomentumCoefficient, fPositionCoefficient;

  std::vector<double> fWeights;

  double fCrossSection, fCrossSectionError;

  int fID1, fID2;
  double fX1, fX2, fScalePDF, fPDF1, fPDF2;

  int fVertexCode, fVertexStatus;
  double fX, fY, fZ, fT;

  int fParticleCode, fPID, fParticleStatus, fOutVertexCode;
  double fPx, fPy, fPz, fE, fMass;

  std::vector<std::pair<TLorentzVector *, TObjArray *> > fVertices;
  std::vector<int> fParticles;

  std::map<int, int> fInVertexMap;
  std::map<int, int> fOutVertexMap;

  std::map<int, std::pair<int, int> > fMotherMap;
  std::map<int, std::pair<int, int> > fDaughterMap;

  HepMC3::ReaderRootTree* fReader;
  HepMC3::GenEvent fEvent;

};

#endif // DelphesHepMC3ROOTReader_h
