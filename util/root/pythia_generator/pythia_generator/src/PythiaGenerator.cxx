#include <pythia_generator/PythiaGenerator.h>
#include <pythia_generator/PythiaToHepMC3.h>

//standard library includes
#include <algorithm> // std::transform

// Pythia8 includes
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "Pythia8/PythiaParallel.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenRunInfo.h"

using namespace std;

namespace PythiaGenerator{

  Generator::Generator(Bool_t parallel){
    createGenerator(parallel);
    _converter = new Pythia8ToHepMC3();
  }

  Generator::~Generator(){
    if (_parallel) delete _pythiaParallel;
    else delete _pythia;
    delete _converter;
    for (auto entry : _events) delete entry;
  }

  void Generator::createGenerator(Bool_t parallel){
    _parallel = parallel;
    if(_parallel){
      if(_instantiated) delete _pythiaParallel;
      _pythiaParallel = new Pythia8::PythiaParallel("",kFALSE);
    }
    else{
      if(_instantiated) delete _pythia;
      _pythia = new Pythia8::Pythia("",kFALSE); // avoid printing the banner
    }
    _instantiated = kTRUE;
  }

  void Generator::readString(TString string){
    if(!_parallel) _pythia->readString(string.Data());
    else _pythiaParallel->readString(string.Data());;
    return;
  }

  void Generator::setQuiet(){
    if(!_parallel) _pythia->readString("Print:quiet = on");
    else _pythiaParallel->readString("Print:quiet = on");
  }

  void Generator::init(){
    if(!_parallel) _pythia->init();
    else _pythiaParallel->init();
    _initialized = kTRUE;
    return;
  }

  void Generator::writeHepMC3File(TString filename){
    if(!_hepmcMode){
      cout << "PythiaGenerator::Generator::writeHepMc3File: Skipping write, HepMC3 mode is off." << endl;
      return; // nothing to write
    }
    if(_events.size() == 0){
      cout << "PythiaGenerator::Generator::writeHepMc3File: No events in buffer." << endl;
      return; // nothing to write
    }

    if(!(_hepmcRootMode || _hepmcAsciiMode)){
      cout << "PythiaGenerator::Generator::writeHepMc3File: Neither ROOT nor ASCII writing turned on." << endl;
      return; // nothing to write
    }

    if(_hepmcRootMode && _hepmcAsciiMode){
      cout << "PythiaGenerator::Generator::writeHepMc3File: Both ROOT nor ASCII writing turned on. This is currently unsupported." << endl;
      return; // nothing to write
    }

    if(_hepmcRootMode){
      HepMC3::WriterRootTree* writer = new HepMC3::WriterRootTree(
        filename.Data(),
        shared_ptr<HepMC3::GenRunInfo>(),
        kTRUE // append mode -- might not be available in HepMC3 via CVMFS, it's a new feature I added. -Jan
      );
      for(HepMC3::GenEvent* event : _events){
        writer->write_event(*event);
      }
      writer->close();
      delete writer;
    }
    else{ // Ascii mode
      HepMC3::WriterAscii* writer = new HepMC3::WriterAscii(filename.Data());
      for(HepMC3::GenEvent* event : _events){
        writer->write_event(*event);
      }
      writer->close();
      delete writer;
    }
    _ClearHepMC3Events();
  }

  void Generator::_ClearParticleContainers(){
    _pid.clear();
    _status.clear();
    _statusHepMC.clear();
    _p.clear();
    _vProd.clear();
    _hasVertex.clear();
    _mass.clear();
    _mother1.clear();
    _mother2.clear();
    _mothers.clear();
    _daughter1.clear();
    _daughter2.clear();
    _daughters.clear();
    return;
  }

  void Generator::_ClearContainers(){
    _ClearParticleContainers();

    _id1Pdf.clear();
    _id2Pdf.clear();
    _pdf1.clear();
    _pdf2.clear();
    _x1Pdf.clear();
    _x2Pdf.clear();
    _QFac.clear();
    _QRen.clear();
    _nMPI.clear();
    _code.clear();
    _alphaS.clear();
    _alphaEM.clear();

    _sigmaGen.clear();
    _sigmaErr.clear();
    _nWeights.clear();
    _weights.clear();

    _nParticlesPerEvent.clear();
    return;
  }

  void Generator::_ClearHepMC3Events(){
    for (auto entry : _events) delete entry;
    _events.clear();
  }

  void Generator::_FillArrays(Pythia8::Pythia* pythia){
    // Deal with particle-level information.
    vector<Int_t> pid = {};
    vector<Int_t> status = {};
    vector<Int_t> statusHepMC = {};
    vector<vector<Double_t>> p = {};
    vector<vector<Double_t>> vProd = {};
    vector<Bool_t> hasVertex = {};
    vector<Double_t> mass = {};
    vector<Int_t> mother1 = {};
    vector<Int_t> mother2 = {};
    vector<vector<Int_t>> mothers = {};

    vector<Int_t> daughter1 = {};
    vector<Int_t> daughter2 = {};
    vector<vector<Int_t>> daughters = {};

    vector<Int_t> col = {};
    vector<Int_t> acol = {};

    for (Int_t j = 1; j < pythia->event.size(); j++){ // Skip entry 0, which represents "the event as a whole"

      pid.push_back(pythia->event[j].id());
      status.push_back(pythia->event[j].status());
      statusHepMC.push_back(pythia->event[j].statusHepMC());

      vector<Double_t> momentum = {pythia->event[j].e(), pythia->event[j].px(), pythia->event[j].py(), pythia->event[j].pz()};
      p.push_back(momentum);

      vector<Double_t> vertex = {pythia->event[j].tProd(), pythia->event[j].xProd(), pythia->event[j].yProd(), pythia->event[j].zProd()};
      vProd.push_back(vertex);
      hasVertex.push_back(pythia->event[j].hasVertex());

      mass.push_back(pythia->event[j].m());

      // For mother and daughter info, convert from 1-indexing to 0-indexing
      mother1.push_back(pythia->event[j].mother1() - 1);
      mother2.push_back(pythia->event[j].mother2() - 1);
      mothers.push_back([&](){
        vector<Int_t> ml = pythia->event[j].motherList();
        std::transform(ml.begin(), ml.end(), ml.begin(), [](Int_t x) { return x - 1; });
        return ml;
      }());

      daughter1.push_back(pythia->event[j].daughter1() - 1);
      daughter2.push_back(pythia->event[j].daughter2() - 1);
      daughters.push_back([&](){
        vector<Int_t> dl = pythia->event[j].daughterList();
        std::transform(dl.begin(), dl.end(), dl.begin(), [](Int_t x) { return x - 1; });
        return dl;
      }());

      col.push_back(pythia->event[j].col());
      acol.push_back(pythia->event[j].acol());
    }
    _pid.push_back(pid);
    _status.push_back(status);
    _statusHepMC.push_back(statusHepMC);
    _p.push_back(p);
    _vProd.push_back(vProd);
    _hasVertex.push_back(hasVertex);
    _mass.push_back(mass);

    _mother1.push_back(mother1);
    _mother2.push_back(mother2);
    _mothers.push_back(mothers);

    _daughter1.push_back(daughter1);
    _daughter2.push_back(daughter2);
    _daughters.push_back(daughters);

    _col.push_back(col);
    _acol.push_back(acol);

    // Now, deal with some event information.
    _id1Pdf.push_back(pythia->info.id1pdf());
    _id2Pdf.push_back(pythia->info.id2pdf());

    _pdf1.push_back(pythia->info.pdf1());
    _pdf2.push_back(pythia->info.pdf2());

    _x1Pdf.push_back(pythia->info.x1pdf());
    _x2Pdf.push_back(pythia->info.x2pdf());

    _QFac.push_back(pythia->info.QFac());
    _QRen.push_back(pythia->info.QRen());

    _nMPI.push_back(pythia->info.nMPI());
    _code.push_back(pythia->info.code());

    _alphaS.push_back(pythia->info.alphaS());
    _alphaEM.push_back(pythia->info.alphaEM());

    _sigmaGen.push_back(pythia->info.sigmaGen());
    _sigmaErr.push_back(pythia->info.sigmaErr());

    Int_t nWeights = pythia->info.nWeights();
    _nWeights.push_back(nWeights);

    vector<Double_t> weights = {};
    for(Int_t j = 0; j < nWeights; j++){
      weights.push_back(pythia->info.weight(j));
    }
    _weights.push_back(weights);

  }

  void Generator::_FillHepMC3Event(Pythia8::Pythia* pythia){
    HepMC3::GenEvent* evt = new HepMC3::GenEvent();
    _converter->fill_next_event(*pythia,evt);
    _events.push_back(evt);
  }

  void Generator::generate(Int_t nevents, Bool_t refresh){
    if(_parallel){
      generateParallel(nevents,refresh);
      return;
    }

    if(!_initialized){
      cout << "Error: PythiaGenerator::Generator is not yet initialized." << endl;
      return;
    }

    // For array mode, we need to clear the arrays if requested.
    if(refresh && _arrayMode) _ClearContainers();

    for(Int_t i = 0; i < nevents; i++){
      Bool_t stat = _pythia->next();

      // In array mode, we fill the arrays with the various particle/event attributes.
      if(_arrayMode) _FillArrays(_pythia);

      // in HepMC3 mode, we will fill a vector of HepMC3 events, which we can later flush to a file.
      if(_hepmcMode) _FillHepMC3Event(_pythia);
    }

    // In array mode, we take the opportunity to compute the number of particles per event.
    if(_arrayMode) _nParticlesPerEvent = get2VectorCounts(_status);
    return;
  }

  void Generator::generateParallel(Int_t nevents, Bool_t refresh){
    if(!_initialized){
      cout << "Error: PythiaGenerator::Generator is not yet initialized." << endl;
      return;
    }
    // For array mode, we need to clear the arrays if requested.
    if(refresh && _arrayMode) _ClearContainers();

    // Generate the events using PythiaParallel
    _pythiaParallel->run(
      nevents, 
      [&](Pythia8::Pythia* pythiaPtr) {
        // In array mode, we fill the arrays with the various particle/event attributes.
      if(_arrayMode) _FillArrays(pythiaPtr);

      // in HepMC3 mode, we will fill a vector of HepMC3 events, which we can later flush to a file.
      if(_hepmcMode) _FillHepMC3Event(pythiaPtr);
    });

    // In array mode, we take the opportunity to compute the number of particles per event.
    if(_arrayMode) _nParticlesPerEvent = get2VectorCounts(_status);
    return;
  }

  vector<vector<Double_t>> Generator::getComponentArray(vector<vector<vector<Double_t>>> inputArray, Int_t index){
      vector<vector<Double_t>> outputArray;
      outputArray.reserve(inputArray.size());

      std::transform(inputArray.begin(), inputArray.end(), std::back_inserter(outputArray),
        [index](const auto& event) {  // Capture index in outer lambda
          vector<Double_t> components;
          components.reserve(event.size());
          std::transform(event.begin(), event.end(), std::back_inserter(components),
            [index](const auto& p) { return p.empty() ? 0.0 : p[index]; });  // Capture index in inner lambda
          return components;
        });
      return outputArray;
  }

  template<typename T>
  vector<T> Generator::get2VectorFlat(const vector<vector<T>>& nested_data) {
    vector<T> flat_data;
    for (const auto& event : nested_data) {
      for (const T& particle_data : event) {
        flat_data.push_back(particle_data);
      }
    }
    return flat_data;
  }

  template<typename T>
  vector<Int_t> Generator::get2VectorCounts(const vector<vector<T>>& nested_data) {
    vector<Int_t> counts;
    for (const auto& event : nested_data) {
      counts.push_back(event.size());
    }
    return counts;
  }

  template<typename T>
  vector<T> Generator::get3VectorFlat(const vector<vector<vector<T>>>& nested_data) {
    vector<T> flat_data;
    for (const auto& event : nested_data) {
      for (const auto& particle_data : event) {
        for (const T& value : particle_data) {
          flat_data.push_back(value);
        }
      }
    }
    return flat_data;
  }

  template<typename T>
  vector<Int_t> Generator::get3VectorCounts(const vector<vector<vector<T>>>& nested_data) {
    vector<Int_t> counts;
    for (const auto& event : nested_data) {
      for (const auto& particle_data : event) {
        counts.push_back(particle_data.size());
      }
    }
    return counts;
  }


  template vector<Int_t> Generator::get2VectorCounts<Int_t>(const vector<vector<Int_t>>& nested_data);
  template vector<Int_t> Generator::get2VectorCounts<Double_t>(const vector<vector<Double_t>>& nested_data);
  template vector<Int_t> Generator::get2VectorCounts<Bool_t>(const vector<vector<Bool_t>>& nested_data);
  template vector<Int_t> Generator::get2VectorFlat<Int_t>(const vector<vector<Int_t>>& nested_data);
  template vector<Double_t> Generator::get2VectorFlat<Double_t>(const vector<vector<Double_t>>& nested_data);
  template vector<Bool_t> Generator::get2VectorFlat<Bool_t>(const vector<vector<Bool_t>>& nested_data);

  template vector<Int_t> Generator::get3VectorCounts<Int_t>(const vector<vector<vector<Int_t>>>& nested_data);
  template vector<Int_t> Generator::get3VectorCounts<Double_t>(const vector<vector<vector<Double_t>>>& nested_data);
  template vector<Int_t> Generator::get3VectorFlat<Int_t>(const vector<vector<vector<Int_t>>>& nested_data);
  template vector<Double_t> Generator::get3VectorFlat<Double_t>(const vector<vector<vector<Double_t>>>& nested_data);
}