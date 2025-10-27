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

  PileupOverlay::PileupOverlay(){
    _InitRNG();
  }

  PileupOverlay::~PileupOverlay(){
    if(_muDistribution) delete _muDistribution;
    delete _rng;
    delete _rngAdapter;
  }

  void PileupOverlay::Initialize(){
    _InitializeIndexMap();
  }

  void PileupOverlay::_InitRNG(){
    if(_rng) delete _rng;
    _rng = new TRandomMixMax17(_rngSeed);
    _rngAdapter = new TRandomAdapter(_rng);
  }

  void PileupOverlay::SetRNGSeed(Int_t rngSeed){
    _rngSeed = rngSeed;
    _InitRNG();
  }

  void PileupOverlay::SetHTCondorInfo(Bool_t flag, Int_t jobNumber, Int_t nJobs){
    _condorFlag = flag;
    _jobNumber = jobNumber;
    _nJobs = nJobs;
  }

  void PileupOverlay::SetBeamSpotSigma(Double_t dt, Double_t dx, Double_t dy, Double_t dz){
    _beamSpotSigma = {dt, dx, dy, dz};
  }

  void PileupOverlay::_InitializeIndexMap(){
    _indexingMap = {};
    _nPileupEvents = 0;

    ULong_t indexOffset = 0;
    for(TString filename : _pileupFilenames){
      // Assuming the input pileup files are ROOT format,
      // so we can figure out how many events they have easily
      // using ROOT. For ASCII this would be quite a mess, and
      // likely much slower.
      TFile* f = new TFile(filename,"READ");
      // cout << filename << endl;
      // f->ls();
      TTree* t = dynamic_cast<TTree*>(f->Get("hepmc3_tree"));
      if(!t){
        f->Close();
        continue;
      }
      ULong_t nentries = t->GetEntries();
      _indexingMap[filename] = make_pair(indexOffset, indexOffset + nentries - 1);
      indexOffset += nentries;
      _nPileupEvents += nentries;
      f->Close();
    }
    _indexingMapInitialized = kTRUE;

    // This is also the place where we'll initialize the index mask.
    _ResetMask();

    return;
  }

  void PileupOverlay::_ResetMask(){
    _pileupEventIndicesMask = vector<Bool_t>(_nPileupEvents, kTRUE);
    return;
  }

  void PileupOverlay::InitMuDistribution(Double_t muAvg, Double_t muSigma){
    if(_muDistribution) delete _muDistribution;

    _muAvg = muAvg;
    _muSigma = muSigma;

    // Default initialization
    Int_t nBins = (Int_t)(_muAvg + 4. * _muSigma);
    TString name = "PileupOverlay_mu";
    _muDistribution = new TH1D(name,";#mu#Relative Count",nBins,0.,(Double_t)nBins);
    for(Int_t i = 0; i < nBins; i++){
      Double_t binCenter = _muDistribution->GetBinCenter(i+1);
      _muDistribution->SetBinContent(i+1, TMath::Gaus(binCenter,_muAvg,_muSigma));
    }
    _muInitialized = kTRUE;
    return;
  }

  void PileupOverlay::InitMuDistribution(TH1D* muDistributionHistogram){
    _muDistribution = new TH1D(*muDistributionHistogram);
    _muDistribution->SetName("PileupOverlay_mu");
    _muDistribution->SetTitle("");
    _muDistribution->GetXaxis()->SetTitle(";#mu");
    _muDistribution->GetYaxis()->SetTitle(";#Relative Count");
    _muInitialized = kTRUE;
    return;
  }

  Int_t PileupOverlay::_SampleMuDistribution(){
    if(!_muInitialized){
      return 0;
    }
    Int_t result = (Int_t)_muDistribution->GetRandom();
    if(result < 0) result = 0;
    return result;
  }

  void PileupOverlay::_PickEventIndices(Int_t nEvents){

    // Do "reservoir sampling"
    _selectedGlobalPileupIndices.clear();
    _selectedGlobalPileupIndices.reserve(nEvents);

    Int_t seen = 0;
    for(ULong_t i = 0; i < _nPileupEvents; i++){
      if(!_pileupEventIndicesMask.at(i)) continue;
      if(seen < nEvents){
        _selectedGlobalPileupIndices.push_back(i);
      }
      else{
        // Int_t j = _mersenneTwister->
        Int_t j = _rng->Integer(seen + 1); // TODO: Check
        if(j < nEvents) _selectedGlobalPileupIndices[j] = i;
      }
      seen++;
    }

    if(_selectedGlobalPileupIndices.size() != nEvents){
      cout << "Warning: Running out of pileup, will recycle." << endl;
      _ResetMask();
      _PickEventIndices(nEvents);
    }
    std::sort(_selectedGlobalPileupIndices.begin(), _selectedGlobalPileupIndices.end());

    // If we don't allow reuse of indices, now we mask out the ones we picked.
    if(!_allowReuse) _UpdateMask();
  }

  void PileupOverlay::_UpdateMask(){
    for(ULong_t idx : _selectedGlobalPileupIndices) _pileupEventIndicesMask[idx] = kFALSE;
  }

  void PileupOverlay::_FetchEventSingleFile(const TString& filename, const vector<ULong_t>& localIndices, vector<HepMC3::GenEvent*> &events){
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

  map<TString, vector<ULong_t>> PileupOverlay::_GetLocalIndices(){
    map<TString, vector<ULong_t>> localIndexMapping = {};
    for(ULong_t globalIdx : _selectedGlobalPileupIndices){
      for (const auto& entry : _indexingMap){
        if(globalIdx >= entry.second.first && globalIdx <= entry.second.second){
          if(localIndexMapping.find(entry.first) == localIndexMapping.end()){
            localIndexMapping[entry.first] = {};
          }
          localIndexMapping[entry.first].push_back(globalIdx - entry.second.first);
          break;
        }
      }
    }
    return localIndexMapping;
  }

  vector<HepMC3::GenEvent*> PileupOverlay::_FetchEvents(){
    // Convert from the globalIndices to filenames and localIndices
    map<TString, vector<ULong_t>> indexMap = _GetLocalIndices(); // the global indices are sorted, so local index lists will be too
    vector<HepMC3::GenEvent*> evts = {};
    for (const auto& entry : indexMap){
      _FetchEventSingleFile(entry.first, entry.second, evts);
    }
    return evts;
  }

  void PileupOverlay::_CombineEventsWithPileup(vector<HepMC3::GenEvent*> &events){
    // This function takes in a vector of input events, and
    // adds pileup to all of them. Batching things this way
    // may be more efficient than going event-by-event, since
    // it (may) reduce the total amount of I/O (incl. lookup
    // in the pileup files), at the cost of increased memory usage.

    Size_t nEvents = events.size();

    // Sample mu for each input event, and keep track of how many
    // pileup events in total we're going to need to fetch.
    vector<Int_t> muValues = {};
    Int_t nPileupInBatch = 0;
    for(Size_t i = 0; i < nEvents; i++){
      Int_t mu = _SampleMuDistribution();
      muValues.push_back(mu);
      nPileupInBatch += mu;
    }

    // Edge case: nPileupInBatch > _nPileupEvents
    if(nPileupInBatch > _nPileupEvents){
      cout << "Error: Requesting more pileup events than are in the input pileup files." << endl;
      return;
    }

    // Fetch all the pileup indices we need at once.
    // The indices are sorted sequentially, which may help
    // with the file I/O (more likely to pick up multiple
    // events within a single basket/batch of the TTrees).
    _PickEventIndices(nPileupInBatch);

    // Now fetch pileup events.
    cout << "Fetching a batch of pilep events..." <<endl;
    vector<HepMC3::GenEvent*> pileupEvents = _FetchEvents();
    cout << "\tDone." << endl;
    // Shuffle the pileup event vector.
    // (Without the shuffle, events are listed in blocks corresponding
    //  with the list of input files -- that maybe aren't random).
    shuffle(pileupEvents.begin(),pileupEvents.end(), *_rngAdapter);

    // Now loop thru the input events, and combine each with a set of pileup events,
    // using the mu values in muValues. The input events are modified in-place.
    Int_t muSum = 0;
    for(Size_t i = 0; i < nEvents; i++){
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

      // displace the main event
      vector<HepMC3::GenVertexPtr> evtVertices = events.at(i)->vertices();
      for(HepMC3::GenVertexPtr vtx : evtVertices){
        const HepMC3::FourVector& oldPosition = vtx->position();
        HepMC3::FourVector new_position(
          oldPosition.x() + mainDisplacement[1],
          oldPosition.y() + mainDisplacement[2],
          oldPosition.z() + mainDisplacement[3],
          oldPosition.t() + mainDisplacement[0]
        );
        vtx->set_position(new_position);
      }

      cout << Form("Adding pileup to event %i/%i",(Int_t)i + 1, (Int_t)nEvents) << endl;
      _AddPileup(
        events.at(i),
        pileupEventsSingle,
        pileupDisplacements,
        pileupPhiRotations
      );
      muSum += mu;
    }
  }

  void PileupOverlay::_AddPileup(HepMC3::GenEvent* evt, const vector<HepMC3::GenEvent*>& pileup, const vector<vector<Double_t>>& pileupDisplacements, const vector<Double_t>& pileupPhiRotations){
    // Loop over pileup events.
    for(Size_t i = 0; i < pileup.size(); i++){
      // cout << Form("\tAdding pileup %i/%i to the event", (Int_t)i + 1, (Int_t)pileup.size()) << endl;
      _AddPileupSingleB(evt, pileup[i], pileupDisplacements[i], pileupPhiRotations[i]);
    }
    return;
  }

  vector<Double_t> PileupOverlay::_generateDisplacement(){
    return {
      _rng->Gaus(0.,_beamSpotSigma[0]),
      _rng->Gaus(0.,_beamSpotSigma[1]),
      _rng->Gaus(0.,_beamSpotSigma[2]),
      _rng->Gaus(0.,_beamSpotSigma[3]),
    };
  }

  void PileupOverlay::_AddPileupSingle(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation){
    // Create mapping from old vertex pointers to new vertex shared pointers
    unordered_map<HepMC3::ConstGenVertexPtr, HepMC3::GenVertexPtr> vertex_map;

    // Get all vertices from pileup event
    vector<HepMC3::GenVertexPtr> vertices_list = pileup_evt->vertices();
    vertex_map.reserve(vertices_list.size());

    ROOT::Math::RotationZ rotation(phiRotation);

    // Copy the pileup vertices, adding displacements to them.
    for (const auto& vertex : vertices_list) {
      HepMC3::FourVector old_position = vertex->position();

      // Apply rotation
      vector<Double_t> new_position_coords = _rotateVectorPhi(&old_position,rotation);

      // Apply displacement
      for(Size_t i = 0; i < 4; i++){
        new_position_coords[i] += pileupDisplacement[i];
      }

      // Note the (x, y, z, t) order for constructor
      HepMC3::FourVector new_position(new_position_coords[1], new_position_coords[2],
                              new_position_coords[3], new_position_coords[0]);

      auto new_vertex = make_shared<HepMC3::GenVertex>(new_position);
      new_vertex->set_status(vertex->status());

      // Store mapping
      vertex_map[vertex] = new_vertex;
    }

    // Copy particles from pileup event, and add them to the vertices (also from pileup event).
    for (const auto& particle : pileup_evt->particles()) {

      // Filter stable particles if requested.
      if(_stableOnly && particle->status() != 1) continue;

      // Get old momentum
      HepMC3::FourVector old_momentum = particle->momentum();

      // Apply rotation
      vector<Double_t> new_momentum_coords = _rotateVectorPhi(&old_momentum,rotation);

      // Create new particle. Note the (px, py, pz, e) order for constructor
      HepMC3::FourVector new_momentum(new_momentum_coords[1], new_momentum_coords[2],
                              new_momentum_coords[3], new_momentum_coords[0]);

      auto new_particle = make_shared<HepMC3::GenParticle>(new_momentum,particle->pid(),particle->status());
      new_particle->set_generated_mass(particle->generated_mass());

      // Set production vertex if it exists
      if (particle->production_vertex()) {
        auto it = vertex_map.find(particle->production_vertex());
        if(it != vertex_map.end()) it->second->add_particle_out(new_particle);
        // If not found, it's expected (e.g., beam particles) - no action needed
      }

      // Set end vertex if it exists
      if (particle->end_vertex()) {
        auto it = vertex_map.find(particle->end_vertex());
        if(it != vertex_map.end()) {
          it->second->add_particle_in(new_particle);
        }
        else{
          cout << "Warning: Could not find end vertex for particle " << particle->pid() << endl;
        }
      }
    }

    // Add all the (new) pileup vertices to the target event.
    // They have particles attached to them.
    for (const auto& pair : vertex_map) {
      evt->add_vertex(pair.second);
    }
    return;
  }



  void PileupOverlay::_AddPileupSingleB(HepMC3::GenEvent* evt, HepMC3::GenEvent* pileup_evt, const vector<Double_t>& pileupDisplacement, const Double_t& phiRotation){
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
    vector<Int_t> accessedIDs = {};

    for(auto particle : pileupParticles){

      // Deal with the origin vertex.
      HepMC3::GenVertexPtr prodVertex = particle->production_vertex();
      if(prodVertex){
        Int_t id = prodVertex->id();
        // Only modify this vertex if we haven't modified it already.
        if(find(accessedIDs.begin(),accessedIDs.end(),id) != accessedIDs.end()){
          HepMC3::FourVector oldPosition = prodVertex->position();
          vector<Double_t> newPositionCoords = _rotateAndTranslateVector(&oldPosition, pileupDisplacement, rotation);
          HepMC3::FourVector newPosition(newPositionCoords[1],newPositionCoords[2],newPositionCoords[3],newPositionCoords[0]);
          prodVertex->set_position(newPosition);
          accessedIDs.push_back(id);
        }
      }

      // Deal with the decay vertex.
      HepMC3::GenVertexPtr endVertex = particle->end_vertex();
      if(endVertex){
        Int_t id = endVertex->id();
        // Only modify this vertex if we haven't modified it already.
        if(find(accessedIDs.begin(),accessedIDs.end(),id) != accessedIDs.end()){
          HepMC3::FourVector oldPosition = prodVertex->position();
          vector<Double_t> newPositionCoords = _rotateAndTranslateVector(&oldPosition, pileupDisplacement, rotation);
          HepMC3::FourVector newPosition(newPositionCoords[1],newPositionCoords[2],newPositionCoords[3],newPositionCoords[0]);
          prodVertex->set_position(newPosition);
          accessedIDs.push_back(id);
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





















  vector<Double_t> PileupOverlay::_rotateVectorPhi(HepMC3::FourVector* vector, const ROOT::Math::RotationZ& rotation){
    if(!_allowPhiRotations) return {vector->t(), vector->x(), vector->y(), vector->z()};

    ROOT::Math::XYZVector spatial(vector->x(), vector->y(), vector->z());
    auto rotatedVec = rotation * spatial;

    return {vector->t(), rotatedVec.X(), rotatedVec.Y(), rotatedVec.Z()};
  }

  vector<Double_t> PileupOverlay::_rotateAndTranslateVector(HepMC3::FourVector* v, const vector<Double_t>& displacementCoordinates, const ROOT::Math::RotationZ& rotation){
    // rotation first, then translation
    vector<Double_t> rotatedCoords = _rotateVectorPhi(v,rotation);
    return {
      rotatedCoords[1] + displacementCoordinates[1],
      rotatedCoords[2] + displacementCoordinates[2],
      rotatedCoords[3] + displacementCoordinates[3],
      rotatedCoords[0] + displacementCoordinates[0],
    };
  }



  void PileupOverlay::operator()(TString inputFile, TString outputFile){

    // Make sure we have a mu distribution of some kind initialized.
    if(!_muInitialized){
      cout << "Warning: Mu distribution not initialized. Will fall back on default Gaussian." << endl;
      InitMuDistribution();
    }

    // Make sure we've indexed the pileup files.
    if(!_indexingMapInitialized) _InitializeIndexMap();

    // Open the input file.
    HepMC3::ReaderRootTree reader(inputFile.Data());
    Bool_t status = !reader.failed();

    // Prepare the writer.
    if(_writer) delete _writer;
    _writer = new HepMC3::WriterRootTree(outputFile.Data(),reader.run_info()); // TODO: Is it OK to just fetch the old run info?

    vector<HepMC3::GenEvent*> eventBuffer = {};

    ULong_t counter = 0;
    while(status){

      if((Int_t)eventBuffer.size() == _batchSize){
        // add in the pileup
        _CombineEventsWithPileup(eventBuffer);

        // flush the buffer
        _Flush(eventBuffer);

        for(auto entry: eventBuffer) delete entry;
        eventBuffer.clear();
      }

      HepMC3::GenEvent* evt = new HepMC3::GenEvent();
      status = reader.read_event(*evt);
      if(!status){
        break;
      }
      eventBuffer.push_back(evt);
      counter++;
    }

    // One more flush for any stragglers
    _CombineEventsWithPileup(eventBuffer);
    _Flush(eventBuffer);
    reader.close();
    _writer->close();
    delete _writer; // safe to destroy it here

  }

  void PileupOverlay::_Flush(vector<HepMC3::GenEvent*> events){
    for(HepMC3::GenEvent* evt : events){
      _writer->write_event(*evt);
    }
    return;
  }

}