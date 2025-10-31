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

/** \class DelphesHepMC3ROOTReader
 *
 *  Reads HepMC3 file, in ROOT format
 *
 *  \author Jan T. Offermann - Brown University
 * Based on the DelphesHepMC3Reader class, written by P. Demin (UCL)
 *
 */

#include "DelphesHepMC3ROOTReader.h"

// Standard library includes
#include <iostream>
#include <sstream>
#include <stdexcept>

// ROOT includes
#include "TObjString.h"

// DELPHES includes
#include "classes/DelphesStream.h"

// HepMC3 includes
#include "HepMC3/Print.h"

using namespace std;

static const Int_t kBufferSize = 16384;

//---------------------------------------------------------------------------

DelphesHepMC3ROOTReader::DelphesHepMC3ROOTReader() :
  fInputFile(0), fBuffer(0), fPDG(0),
  fVertexCounter(-2), fParticleCounter(-1)
{
  fBuffer = new Char_t[kBufferSize];

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesHepMC3ROOTReader::~DelphesHepMC3ROOTReader(){
  if(fBuffer) delete[] fBuffer;
  if(fReader){
    fReader->close();
    delete fReader;
  }
  if(fReaderPileup){
    fReaderPileup->close();
    delete fReaderPileup;
  }
  // if(fEvent) delete fEvent;
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::SetInputFile(TString inputFile){
  fInputFile = inputFile;
  fReader = new HepMC3::ReaderRootTree(string(fInputFile.Data()));
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::SetInputPileupFile(TString inputPileupFile){
  if(inputPileupFile.EqualTo("")){
    fHasPileupFile = kFALSE;
    return;
  }
  fInputFilePileup = inputPileupFile;
  fReaderPileup = new HepMC3::ReaderRootTree(string(fInputFilePileup.Data()));
  fHasPileupFile = kTRUE;
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::Clear(){
  fWeights.clear();
  fMomentumCoefficient = 1.0;
  fPositionCoefficient = 1.0;
  fVertexCounter = -2;
  fParticleCounter = -1;
  fVertices.clear();
  fParticles.clear();
  fInVertexMap.clear();
  fOutVertexMap.clear();
  fMotherMap.clear();
  fDaughterMap.clear();
}

//---------------------------------------------------------------------------

Bool_t DelphesHepMC3ROOTReader::EventReady(){
  return kTRUE;
  //return (fVertexCounter == -1) && (fParticleCounter == 0);
}

//---------------------------------------------------------------------------

Bool_t DelphesHepMC3ROOTReader::ReadEvent(){
  // if(fEvent != 0) fEvent->clear(); // unnecessary; see GenEvent::read_data() (called by ReaderRootTree::read_event())
  fReader->read_event(fEvent);
  Bool_t status = !fReader->failed();

  if(fHasPileupFile){ // also need to read the pileup file
    fReaderPileup->read_event(fEventPileup);
    status = status && !fReaderPileup->failed();
  }

  return status;
}

void DelphesHepMC3ROOTReader::Analyze(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray){

  // Fetch the particles
  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles = fEvent.particles();

  // this creates "candidates" in the factory - and also puts the particle production vertices into fVertices
  for(auto particle: particles) AnalyzeParticle(factory,particle); // this "loads" the particles into the factory

  // If we have read in a pileup event, we also add that in.
  if(fHasPileupFile){
    std::vector<std::shared_ptr<HepMC3::GenParticle>> pileupParticles = fEventPileup.particles();
    for(auto particle: pileupParticles) AnalyzeParticle(factory,particle, kTRUE);
  }

  // capturing weights -- TODO: This might need some work, admittedly I don't fully understand this in the old code. -Jan
  fWeights = fEvent.weights();

  // Re-using some existing functions. This puts things into the TObjArrays.
  FinalizeParticles(allParticleOutputArray,stableParticleOutputArray,partonOutputArray);
}

void DelphesHepMC3ROOTReader::AnalyzeEvent(ExRootTreeBranch *branch, long long /*eventNumber*/,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  HepMCEvent *element;

  element = static_cast<HepMCEvent *>(branch->NewEntry());

  // Assign a bunch of things -- lets the structure of the code below stay the same. -Jan
  fEventNumber = fEvent.event_number();
  fProcessID = stoi(fEvent.attribute_as_string("signal_process_id")); // there *must* be a much better way to do this
  fMPI = stoi(fEvent.attribute_as_string("mpi")); // there *must* be a much better way to do this
  fCrossSection = fEvent.cross_section()->xsec();
  fCrossSectionError = fEvent.cross_section()->xsec_err();
  fScale = stof(fEvent.attribute_as_string("event_scale"));
  fAlphaQED = stof(fEvent.attribute_as_string("alphaQED"));
  fAlphaQCD = stof(fEvent.attribute_as_string("alphaQCD"));

  // some things that come from GenPdfInfo
  TString genPdfInfo = fEvent.attribute_as_string("GenPdfInfo");
  TObjArray* genPdfInfoArray = genPdfInfo.Tokenize(" ");
  fID1      = ((TObjString*)(genPdfInfoArray->At(0)))->String().Atoi();
  fID2      = ((TObjString*)(genPdfInfoArray->At(1)))->String().Atoi();
  fX1       = ((TObjString*)(genPdfInfoArray->At(2)))->String().Atof();
  fX2       = ((TObjString*)(genPdfInfoArray->At(3)))->String().Atof();
  fScalePDF = ((TObjString*)(genPdfInfoArray->At(4)))->String().Atof();
  fPDF1     = ((TObjString*)(genPdfInfoArray->At(5)))->String().Atof();
  fPDF2     = ((TObjString*)(genPdfInfoArray->At(6)))->String().Atof();
  delete genPdfInfoArray;

  element->Number = fEventNumber;

  element->ProcessID = fProcessID;
  element->MPI = fMPI;
  element->Weight = fWeights.size() > 0 ? fWeights[0] : 1.0;
  element->CrossSection = fCrossSection;
  element->CrossSectionError = fCrossSectionError;
  element->Scale = fScale;
  element->AlphaQED = fAlphaQED;
  element->AlphaQCD = fAlphaQCD;

  element->ID1 = fID1;
  element->ID2 = fID2;
  element->X1 = fX1;
  element->X2 = fX2;
  element->ScalePDF = fScalePDF;
  element->PDF1 = fPDF1;
  element->PDF2 = fPDF2;

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::AnalyzeWeight(ExRootTreeBranch *branch)
{
  Weight *element;
  vector<Double_t>::const_iterator itWeights;

  for(itWeights = fWeights.begin(); itWeights != fWeights.end(); ++itWeights){
    element = static_cast<Weight *>(branch->NewEntry());

    element->Weight = *itWeights;
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::AnalyzeVertex(DelphesFactory *factory, Int_t code, Candidate *candidate)
{
  Int_t index;
  TLorentzVector *position;
  TObjArray *array;
  vector<Int_t>::iterator itParticle;
  map<Int_t, Int_t>::iterator itVertexMap;

  itVertexMap = fOutVertexMap.find(code);
  if(itVertexMap == fOutVertexMap.end()){
    --fVertexCounter;

    index = fVertices.size();
    fOutVertexMap[code] = index;
    if(candidate && code > 0) fInVertexMap[code] = index;

    position = factory->New<TLorentzVector>();
    array = factory->NewArray();
    position->SetXYZT(0.0, 0.0, 0.0, 0.0);
    fVertices.push_back(make_pair(position, array));
  }
  else{
    index = itVertexMap->second;
    position = fVertices[index].first;
    array = fVertices[index].second;
  }

  if(candidate){
    array->Add(candidate);
  }
  else{
    position->SetXYZT(fX, fY, fZ, fT);
    for(itParticle = fParticles.begin(); itParticle != fParticles.end(); ++itParticle){
      fInVertexMap[*itParticle] = index;
    }
  }
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::AnalyzeParticle(DelphesFactory *factory, std::shared_ptr<HepMC3::GenParticle> particle, Bool_t isPileup){
  Candidate *candidate;

  candidate = factory->NewCandidate();

  candidate->PID = particle->pid();

  candidate->Status = particle->status();

  candidate->Mass = particle->generated_mass();

  candidate->Momentum.SetPxPyPzE(
    particle->momentum().px(),
    particle->momentum().py(),
    particle->momentum().pz(),
    particle->momentum().e()
  );

  candidate->D1 = particle->id();
  if(isPileup) candidate->IsPU = 1;

  auto prod_vtx = particle->production_vertex();
  Int_t prod_vertex_id = prod_vtx ? prod_vtx->id() : 0;

  // fills fVertices
  AnalyzeVertex(factory, prod_vertex_id, candidate);
}

//---------------------------------------------------------------------------

void DelphesHepMC3ROOTReader::FinalizeParticles(TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  TLorentzVector *position;
  TObjArray *array;
  Candidate *candidate;
  Candidate *candidateDaughter;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;
  map<Int_t, Int_t >::iterator itVertexMap;
  map<Int_t, pair<Int_t, Int_t> >::iterator itMotherMap;
  map<Int_t, pair<Int_t, Int_t> >::iterator itDaughterMap;
  size_t i;
  Int_t j, code, counter;

  counter = 0;
  for(i = 0; i < fVertices.size(); ++i){
    position = fVertices[i].first;
    array = fVertices[i].second;

    for(j = 0; j < array->GetEntriesFast(); ++j){
      candidate = static_cast<Candidate *>(array->At(j));

      candidate->Position = *position;
      if(fPositionCoefficient != 1.0){
        candidate->Position *= fPositionCoefficient;
      }

      if(fMomentumCoefficient != 1.0){
        candidate->Momentum *= fMomentumCoefficient;
      }

      candidate->M1 = i;

      itDaughterMap = fDaughterMap.find(i);
      if(itDaughterMap == fDaughterMap.end()){
        fDaughterMap[i] = make_pair(counter, counter);
      }
      else{
        itDaughterMap->second.second = counter;
      }

      code = candidate->D1;

      itVertexMap = fInVertexMap.find(code);
      if(itVertexMap == fInVertexMap.end()){
        candidate->D1 = -1;
      }
      else{
        code = itVertexMap->second;

        candidate->D1 = code;

        itMotherMap = fMotherMap.find(code);
        if(itMotherMap == fMotherMap.end()){
          fMotherMap[code] = make_pair(counter, -1);
        }
        else{
          itMotherMap->second.second = counter;
        }
      }

      allParticleOutputArray->Add(candidate);
      ++counter;

      pdgParticle = fPDG->GetParticle(candidate->PID);

      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;

      if(!pdgParticle) continue;

      pdgCode = TMath::Abs(candidate->PID);

      if(candidate->Status == 1){
        stableParticleOutputArray->Add(candidate);
      }
      else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15){
        partonOutputArray->Add(candidate);
      }
    }
  }

  for(j = 0; j < allParticleOutputArray->GetEntriesFast(); ++j){
    candidate = static_cast<Candidate *>(allParticleOutputArray->At(j));

    itMotherMap = fMotherMap.find(candidate->M1);
    if(itMotherMap == fMotherMap.end()){
      candidate->M1 = -1;
      candidate->M2 = -1;
    }
    else{
      candidate->M1 = itMotherMap->second.first;
      candidate->M2 = itMotherMap->second.second;
    }

    if(candidate->D1 < 0){
      candidate->D1 = -1;
      candidate->D2 = -1;
    }
    else{
      itDaughterMap = fDaughterMap.find(candidate->D1);
      if(itDaughterMap == fDaughterMap.end()){
        candidate->D1 = -1;
        candidate->D2 = -1;
        const TLorentzVector &decayPosition = candidate->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T());// decay position
      }
      else{
        candidate->D1 = itDaughterMap->second.first;
        candidate->D2 = itDaughterMap->second.second;
        candidateDaughter = static_cast<Candidate *>(allParticleOutputArray->At(candidate->D1));
        const TLorentzVector &decayPosition = candidateDaughter->Position;
        candidate->DecayPosition.SetXYZT(decayPosition.X(), decayPosition.Y(), decayPosition.Z(), decayPosition.T());// decay position
      }
    }
  }

    // cout << "FinalizeParticles: We have " << allParticleOutputArray->GetEntriesFast() << " entries in allParticleOutputArray" << std::endl;
    // cout << "FinalizeParticles: We have " << stableParticleOutputArray->GetEntriesFast() << " entries in stableParticleOutputArray" << std::endl;
    // cout << "FinalizeParticles: We have " << partonOutputArray->GetEntriesFast() << " entries in partonOutputArray" << std::endl;
    // cout << endl;

    // // print contents of stableParticleOutputArray
    // for(UInt_t i = 0; i < stableParticleOutputArray->GetEntriesFast(); i++){
    //   cout << "\t[" << i << "] " << stableParticleOutputArray->At(0)->
    // }



}

//---------------------------------------------------------------------------
