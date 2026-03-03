#include <pileup/Pileup.h>

//standard library includes
#include <algorithm> //  std::shuffle

// ROOT includes
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/WriterRootTree.h"

using namespace std;

namespace Pileup{

  PileupMixer::PileupMixer(){
    _InitRNG();
  }

  PileupMixer::~PileupMixer(){
    if(_muDistribution) delete _muDistribution;
    delete _rng;
    delete _rngAdapter;
  }

  void PileupMixer::Initialize(){
    _InitializeIndexMap();
  }

  void PileupMixer::_InitRNG(){
    if(_rng) delete _rng;
    _rng = new TRandomMixMax17(_rngSeed);
    _rngAdapter = new TRandomAdapter(_rng);
  }

  void PileupMixer::SetRNGSeed(Int_t rngSeed){
    _rngSeed = rngSeed;
    _InitRNG();
  }

  void PileupMixer::SetHTCondorInfo(Bool_t flag, Int_t jobNumber, Int_t nJobs){
    _condorFlag = flag;
    _jobNumber = jobNumber;
    _nJobs = nJobs;
  }

  void PileupMixer::SetBeamSpotSigma(Double_t dt, Double_t dx, Double_t dy, Double_t dz){
    _beamSpotSigma = {dt, dx, dy, dz};
  }

  void PileupMixer::_InitializeIndexMap(){
      _fileRanges.clear();
      _nPileupEvents = 0;

      for(TString filename : _pileupFilenames){
          TFile* f = new TFile(filename, "READ");
          TTree* t = dynamic_cast<TTree*>(f->Get("hepmc3_tree"));
          if(!t){ f->Close(); delete f; continue; }

          ULong_t nentries = t->GetEntries();
          _fileRanges.push_back({_nPileupEvents, _nPileupEvents + nentries - 1, filename});
          _nPileupEvents += nentries;
          f->Close();
          delete f;
      }
      // Already sorted by construction since we iterate files in order
      // and assign ranges sequentially, but sort explicitly for safety:
      std::sort(_fileRanges.begin(), _fileRanges.end(),
                [](const FileRange& a, const FileRange& b){ return a.start < b.start; });

      _indexingMapInitialized = kTRUE;
      _usedPileupIndices.clear();
  }

  void PileupMixer::InitMuDistribution(Double_t muAvg, Double_t muSigma){
    if(_muDistribution) delete _muDistribution;

    _muAvg = muAvg;
    _muSigma = muSigma;

    // Default initialization
    Int_t nBins = (Int_t)(_muAvg + 4. * _muSigma);
    TString name = "PileupOverlay_mu";
    _muDistribution = new TH1D(name,";#mu#Relative Count",nBins,0.,(Double_t)nBins);
    _muDistribution->SetDirectory(0);
    for(Int_t i = 0; i < nBins; i++){
      Double_t binCenter = _muDistribution->GetBinCenter(i+1);
      _muDistribution->SetBinContent(i+1, TMath::Gaus(binCenter,_muAvg,_muSigma));
    }
    _muInitialized = kTRUE;
    return;
  }

  void PileupMixer::InitMuDistribution(TH1D* muDistributionHistogram){
    _muDistribution = new TH1D(*muDistributionHistogram);
    _muDistribution->SetDirectory(0);
    _muDistribution->SetName("PileupOverlay_mu");
    _muDistribution->SetTitle("");
    _muDistribution->GetXaxis()->SetTitle(";#mu");
    _muDistribution->GetYaxis()->SetTitle(";#Relative Count");
    _muInitialized = kTRUE;
    return;
  }

  Int_t PileupMixer::_SampleMuDistribution(){
    if(!_muInitialized){
      return 0;
    }
    Int_t result = (Int_t)_muDistribution->GetRandom();
    if(result < 0) result = 0;
    return result;
  }

  void PileupMixer::_PickEventIndices(Int_t nEvents){
      _selectedGlobalPileupIndices.clear();
      _selectedGlobalPileupIndices.reserve(nEvents);

      if(_useContiguousSampling){
          // Pick a random start index and take nEvents contiguous indices, wrapping around.
          ULong_t start = _rng->Integer(_nPileupEvents);
          for(Int_t i = 0; i < nEvents; i++){
              _selectedGlobalPileupIndices.push_back((start + i) % _nPileupEvents);
          }
          // Sort so that per-file index lists are sequential, for I/O performance
          std::sort(_selectedGlobalPileupIndices.begin(), _selectedGlobalPileupIndices.end());
      }
      else if(_allowReuse){
          for(Int_t i = 0; i < nEvents; i++){
              _selectedGlobalPileupIndices.push_back(_rng->Integer(_nPileupEvents));
          }
          std::sort(_selectedGlobalPileupIndices.begin(), _selectedGlobalPileupIndices.end());
      }
      else{ // the "old" case -- this could potentially be much slower
          ULong_t nRemaining = _nPileupEvents - _usedPileupIndices.size();
          if((ULong_t)nEvents > nRemaining){
              if(_nWarning < _nWarningMax){
                  cout << "Warning: Running low on unused pileup events, resetting." << endl;
                  _nWarning++;
              }
              _usedPileupIndices.clear();
          }
          while((Int_t)_selectedGlobalPileupIndices.size() < nEvents){
              ULong_t candidate = _rng->Integer(_nPileupEvents);
              if(_usedPileupIndices.find(candidate) == _usedPileupIndices.end()){
                  _selectedGlobalPileupIndices.push_back(candidate);
                  _usedPileupIndices.insert(candidate);
              }
          }
          std::sort(_selectedGlobalPileupIndices.begin(), _selectedGlobalPileupIndices.end());
      }
  }

  void PileupMixer::_FetchEventSingleFile(const TString& filename, const vector<ULong_t>& localIndices, vector<HepMC3::GenEvent*> &events){
    HepMC3::ReaderRootTree reader(filename.Data()); // Note the use of HepMC3/ROOT format! Can optionally pass tree name but it is standard

    for(ULong_t idx : localIndices){
      HepMC3::GenEvent* evt = new HepMC3::GenEvent();
      Bool_t status = reader.read_event_at_index(*evt,idx);
      if(reader.failed() || !status){
        break;
      }
      events.push_back(evt);
    }
    reader.close();
    return;
  }

  vector<HepMC3::GenEvent*> PileupMixer::_FetchEvents(){
      vector<vector<ULong_t>> localIndicesByFile(_fileRanges.size());

      for(ULong_t globalIdx : _selectedGlobalPileupIndices){
          auto it = std::upper_bound(
              _fileRanges.begin(), _fileRanges.end(), globalIdx,
              [](ULong_t val, const FileRange& r){ return val < r.start; }
          );
          if(it == _fileRanges.begin()) continue;
          --it;
          localIndicesByFile[std::distance(_fileRanges.begin(), it)].push_back(globalIdx - it->start);
      }

      vector<HepMC3::GenEvent*> evts;
      for(ULong_t i = 0; i < _fileRanges.size(); i++){
          if(!localIndicesByFile[i].empty()){
              _FetchEventSingleFile(_fileRanges[i].filename, localIndicesByFile[i], evts);
          }
      }
      return evts;
  }

  vector<Double_t> PileupMixer::_generateDisplacement(){
    return {
      _rng->Gaus(0.,_beamSpotSigma[0]),
      _rng->Gaus(0.,_beamSpotSigma[1]),
      _rng->Gaus(0.,_beamSpotSigma[2]),
      _rng->Gaus(0.,_beamSpotSigma[3]),
    };
  }

  void PileupMixer::_AddPileupSingle(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation){
    // TODO: This has to be thoroughly tested!

    // Will overlay pileup events onto the main event "evt" using GenEvent::add_tree().
    // NOTE: I previously tried a different approach, where I fetch all the vertices
    //       from the pileup events, add incoming/outgoing particles, and then add
    //       these with GenEvent::add_vertex(). That approach seems *much* slower.

    // Fetch the particles.
    vector<HepMC3::GenParticlePtr> pileupParticles = pileup_evt->particles();

    // Now, we have to modify:
    // 1) The particles' momenta (rotation)
    // 2) The attached vertices' positions (translation + rotation)
    // The tricky part is handling the vertices, we don't want to
    // accidentally adjust the same vertex twice.
    ROOT::Math::RotationZ rotation(phiRotation);
    unordered_set<Int_t> accessedIDs;

    for(auto particle : pileupParticles){

      // Deal with the origin vertex.
      HepMC3::GenVertexPtr prodVertex = particle->production_vertex();
      if(prodVertex){
        Int_t id = prodVertex->id();
        // Only modify this vertex if we haven't modified it already.
        if(accessedIDs.find(id) == accessedIDs.end()){
          HepMC3::FourVector oldPosition = prodVertex->position();
          vector<Double_t> newPositionCoords = _rotateAndTranslateVector(&oldPosition, pileupDisplacement, rotation);
          HepMC3::FourVector newPosition(newPositionCoords[1],newPositionCoords[2],newPositionCoords[3],newPositionCoords[0]);
          prodVertex->set_position(newPosition);
          accessedIDs.insert(id);
        }
      }

      // Deal with the decay vertex.
      HepMC3::GenVertexPtr endVertex = particle->end_vertex();
      if(endVertex){
        Int_t id = endVertex->id();
        // Only modify this vertex if we haven't modified it already.
        if(accessedIDs.find(id) == accessedIDs.end()){
          HepMC3::FourVector oldPosition = endVertex->position();
          vector<Double_t> newPositionCoords = _rotateAndTranslateVector(&oldPosition, pileupDisplacement, rotation);
          HepMC3::FourVector newPosition(newPositionCoords[1],newPositionCoords[2],newPositionCoords[3],newPositionCoords[0]);
          endVertex->set_position(newPosition);
          accessedIDs.insert(id);
        }
      }

      // Deal with particle momentum -- rotation only.
      HepMC3::FourVector oldMomentum = particle->momentum();
      vector<Double_t> rotatedMomentumCoords = _rotateVectorPhi(&oldMomentum,rotation); // TODO: CHECK
      HepMC3::FourVector rotatedMomentum(rotatedMomentumCoords[1],rotatedMomentumCoords[2],rotatedMomentumCoords[3],rotatedMomentumCoords[0]);
      particle->set_momentum(rotatedMomentum);

    }
    evt->add_tree(pileupParticles);

    return;
  }

  void PileupMixer::_AddPileupSingleStableOnly(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation){
    // TODO: This has to be thoroughly tested!

    // Will overlay pileup events onto the main event "evt" using GenEvent::add_tree().
    // Will only take stable pileup particles, and will attach them to a single vertex.
    // As long as they're prompt, this should be OK! We typically don't care about the
    // whole history of the pileup events.
    // Fetch the pileup particles.
    vector<HepMC3::GenParticlePtr> pileupParticles = pileup_evt->particles();

    HepMC3::FourVector vtxPosition(pileupDisplacement[1],pileupDisplacement[2],pileupDisplacement[3],pileupDisplacement[0]);
    HepMC3::GenVertexPtr vtx = make_shared<HepMC3::GenVertex>(vtxPosition);
    ROOT::Math::RotationZ rotation(phiRotation);

    for(auto particle : pileupParticles){
      // Deal with particle momentum -- rotation only.
      if(particle->status() != 1) continue;
      HepMC3::GenParticlePtr p = make_shared<HepMC3::GenParticle>(particle->momentum(), particle->pid(),1);
      HepMC3::FourVector oldMomentum = p->momentum();
      vector<Double_t> rotatedMomentumCoords = _rotateVectorPhi(&oldMomentum,rotation); // TODO: CHECK
      HepMC3::FourVector rotatedMomentum(rotatedMomentumCoords[1],rotatedMomentumCoords[2],rotatedMomentumCoords[3],rotatedMomentumCoords[0]);
      p->set_momentum(rotatedMomentum);
      vtx->add_particle_out(p);
    }

    // cout << "\t_AddPileupSingleStableOnly: Calling add_vertex" << endl;
    evt->add_vertex(vtx);
    // cout << Form("\t\tevt->particles().size() = %i",(Int_t)evt->particles().size()) << endl;
    return;
  }

  vector<Double_t> PileupMixer::_rotateVectorPhi(HepMC3::FourVector* vector, const ROOT::Math::RotationZ& rotation){
    if(!_allowPhiRotations) return {vector->t(), vector->x(), vector->y(), vector->z()};

    ROOT::Math::XYZVector spatial(vector->x(), vector->y(), vector->z());
    auto rotatedVec = rotation * spatial;

    return {vector->t(), rotatedVec.X(), rotatedVec.Y(), rotatedVec.Z()};
  }

  vector<Double_t> PileupMixer::_rotateAndTranslateVector(HepMC3::FourVector* v, const vector<Double_t>& displacementCoordinates, const ROOT::Math::RotationZ& rotation){
    // rotation first, then translation
    vector<Double_t> rotatedCoords = _rotateVectorPhi(v,rotation);
    return {
      rotatedCoords[1] + displacementCoordinates[1],
      rotatedCoords[2] + displacementCoordinates[2],
      rotatedCoords[3] + displacementCoordinates[3],
      rotatedCoords[0] + displacementCoordinates[0],
    };
  }

vector<HepMC3::GenEvent*> PileupMixer::_CreatePileupEvents(Int_t nEvents){
    // This function takes in a number for the total number
    // of combined pileup events to produce, and makes them.
    // Batching things in this way -- not just making a single
    // combined event, but a bunch of them -- may be
    // more efficient than going event-by-event, since it
    // (may) reduce the total amount of I/O (incl. lookup
    // in the pileup files), at the cost of increased memory usage.

    // Sample mu for each input event, and keep track of how many
    // pileup events in total we're going to need to fetch.

    vector<HepMC3::GenEvent*> events = {};

    vector<Int_t> muValues = {};
    Int_t nPileupInBatch = 0;
    Int_t nEventsLocal = -1;
    for(Int_t i = 0; i < nEvents; i++){
      Int_t mu = _SampleMuDistribution();
      nPileupInBatch += mu;
      if(nPileupInBatch > _nPileupEvents){
        nPileupInBatch -= mu;

        nEventsLocal = i;
        break;
      }
      muValues.push_back(mu);
    }

    if(nEventsLocal == 0){
      cout << "Error: Not enough pileup events in input to produce a single event!" << endl;
      return {};
    }
    else if(nEventsLocal > 0){
      events = _CreatePileupEvents(nEvents - nEventsLocal);
      nEvents = nEventsLocal;
    }

    _PickEventIndices(nPileupInBatch);

    // Now fetch pileup events.
    vector<HepMC3::GenEvent*> pileupEvents = _FetchEvents();
    // Shuffle the pileup event vector.
    // (Without the shuffle, events are listed in blocks corresponding
    //  with the list of input files -- that maybe aren't random).
    shuffle(pileupEvents.begin(),pileupEvents.end(), *_rngAdapter);

    // Now loop thru the input events, and combine each with a set of pileup events,
    // using the mu values in muValues. The input events are modified in-place.
    Int_t muSum = 0;
    for(Int_t i = 0; i < nEvents; i++){
      Int_t mu = muValues.at(i);

      vector<HepMC3::GenEvent*> pileupEventsSingle;
      pileupEventsSingle.assign(pileupEvents.begin() + muSum, pileupEvents.begin() + muSum + mu);

      // Here, we'll generate displacements (based on beamspot) for the main event
      // as well as the pileup events, plus random phi rotations for the pileup events.
      vector<Double_t> mainDisplacement = _generateDisplacement();
      vector<vector<Double_t>> pileupDisplacements = {};
      for(Size_t i = 0; i < mu; i++) pileupDisplacements.push_back(_generateDisplacement());

      vector<Double_t> pileupPhiRotations(mu);
      generate(pileupPhiRotations.begin(),pileupPhiRotations.end(), [this](){return _rng->Uniform(2. * TMath::Pi());});

      HepMC3::GenEvent* evt = _CreatePileupEvent(
        pileupEventsSingle,
        pileupDisplacements,
        pileupPhiRotations
      );
      evt->set_event_number(i);
      events.push_back(evt);
      muSum += mu;
    }
    return events;
  }

  HepMC3::GenEvent* PileupMixer::_CreatePileupEvent(const vector<HepMC3::GenEvent*>& pileup, const vector<vector<Double_t>>& pileupDisplacements, const vector<Double_t>& pileupPhiRotations){
    // Create empty GenEvent
    HepMC3::GenEvent* evt = new HepMC3::GenEvent(); // defaults to GeV & mm units
    // cout << Form("_CreatePileupEvent: Created empty event, evt->particles().size() = %i",(Int_t)evt->particles().size()) << endl;

    // Loop over pileup events.
    for(Size_t i = 0; i < pileup.size(); i++){
      // cout << Form("\tAdding pileup %i/%i to the event", (Int_t)i + 1, (Int_t)pileup.size()) << endl;
      if(_stableOnly) _AddPileupSingleStableOnly(evt, pileup[i], pileupDisplacements[i], pileupPhiRotations[i]);
      else _AddPileupSingle(evt, pileup[i], pileupDisplacements[i], pileupPhiRotations[i]);
    }
    // cout << Form("_CreatePileupEvent: Reached the end of function, evt->particles().size() = %i\n",(Int_t)evt->particles().size()) << endl;

    return evt;
  }

  void PileupMixer::operator()(ULong_t nEvents, TString outputFile){

    // Make sure we have a mu distribution of some kind initialized.
    if(!_muInitialized){
      cout << "Warning: Mu distribution not initialized. Will fall back on default Gaussian." << endl;
      InitMuDistribution();
    }

    // Make sure we've indexed the pileup files.
    if(!_indexingMapInitialized) _InitializeIndexMap();

    // Prepare the writer.
    if(_writer) delete _writer;
    _writer = new HepMC3::WriterRootTree(outputFile.Data());

    ULong_t counter = 0;
    while(counter < nEvents){

      ULong_t batchSize = _batchSize;
      if(nEvents - counter < _batchSize) batchSize = nEvents - counter;

      // Create a batch of pileup events.
      vector<HepMC3::GenEvent*> eventBuffer = _CreatePileupEvents(batchSize);

      // flush the buffer.
      _Flush(eventBuffer);
      for(auto entry: eventBuffer) delete entry;
      counter += batchSize;
    }
    _writer->close();
    delete _writer; // safe to destroy it here

  }

  void PileupMixer::_Flush(vector<HepMC3::GenEvent*> events){
    for(HepMC3::GenEvent* evt : events){
      _writer->write_event(*evt);
    }
    return;
  }

  void PileupMixer::operator()(TString inputFile, TString outputFile){

    // Determine the number of events in the input file.
    TFile* f = new TFile(inputFile,"READ");
    TTree* t = (TTree*)f->Get("hepmc3_tree");
    ULong_t nEvents = t->GetEntries();
    f->Close();
    delete f;

    return this->operator()(nEvents,outputFile);
  }
}