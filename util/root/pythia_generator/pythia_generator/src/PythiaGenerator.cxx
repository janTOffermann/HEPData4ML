#include <pythia_generator/PythiaGenerator.h>

//standard library includes
#include <algorithm> // std::transform

// Pythia8 includes
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

// Pythia8 plugin includes -- HepMC3 interface.
// NOTE: In general we don't require Pythia8
//       to be built against our local HepMC3 install,
//       which has some custom features -- we're just
//       going to use "standard" HepMC3/ROOT writing.
#include "Pythia8Plugins/HepMC3.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterRootTree.h"

using namespace std;

namespace PythiaGenerator{

  Generator::Generator(){
    _pythia = new Pythia8::Pythia("",kFALSE); // avoid printing the banner

    _converter = new HepMC3::Pythia8ToHepMC3();
  }

  Generator::~Generator(){
    delete _pythia;
    delete _converter;
    for (auto entry : _events) delete entry;
  }

  void Generator::readString(TString string){
    _pythia->readString(string.Data());
    return;
  }

  void Generator::setQuiet(){
    _pythia->readString("Print:quiet = on");
  }

  void Generator::init(){
    _pythia->init();
    _initialized = kTRUE;
    return;
  }

  void Generator::writeHepMC3File(TString filename){
    if(!_hepmcMode) return; // nothing to write

    HepMC3::WriterRootTree* writer = new HepMC3::WriterRootTree(filename.Data());

    for(HepMC3::GenEvent* event : _events){
      writer->write_event(*event);
    }
    writer->close();
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

  void Generator::_FillArrays(){
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

    for (Int_t j = 1; j < _pythia->event.size(); j++){ // Skip entry 0, which represents "the event as a whole"

      pid.push_back(_pythia->event[j].id());
      status.push_back(_pythia->event[j].status());
      statusHepMC.push_back(_pythia->event[j].statusHepMC());

      vector<Double_t> momentum = {_pythia->event[j].e(), _pythia->event[j].px(), _pythia->event[j].py(), _pythia->event[j].pz()};
      p.push_back(momentum);

      vector<Double_t> vertex = {_pythia->event[j].tProd(), _pythia->event[j].xProd(), _pythia->event[j].yProd(), _pythia->event[j].zProd()};
      vProd.push_back(vertex);
      hasVertex.push_back(_pythia->event[j].hasVertex());

      mass.push_back(_pythia->event[j].m());

      // For mother and daughter info, convert from 1-indexing to 0-indexing
      mother1.push_back(_pythia->event[j].mother1() - 1);
      mother2.push_back(_pythia->event[j].mother2() - 1);
      mothers.push_back([&](){
        vector<Int_t> ml = _pythia->event[j].motherList();
        std::transform(ml.begin(), ml.end(), ml.begin(), [](Int_t x) { return x - 1; });
        return ml;
      }());

      daughter1.push_back(_pythia->event[j].daughter1() - 1);
      daughter2.push_back(_pythia->event[j].daughter2() - 1);
      daughters.push_back([&](){
        vector<Int_t> dl = _pythia->event[j].daughterList();
        std::transform(dl.begin(), dl.end(), dl.begin(), [](Int_t x) { return x - 1; });
        return dl;
      }());

      col.push_back(_pythia->event[j].col());
      acol.push_back(_pythia->event[j].acol());
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
    _id1Pdf.push_back(_pythia->info.id1pdf());
    _id2Pdf.push_back(_pythia->info.id2pdf());

    _pdf1.push_back(_pythia->info.pdf1());
    _pdf2.push_back(_pythia->info.pdf2());

    _x1Pdf.push_back(_pythia->info.x1pdf());
    _x2Pdf.push_back(_pythia->info.x2pdf());

    _QFac.push_back(_pythia->info.QFac());
    _QRen.push_back(_pythia->info.QRen());

    _nMPI.push_back(_pythia->info.nMPI());
    _code.push_back(_pythia->info.code());

    _alphaS.push_back(_pythia->info.alphaS());
    _alphaEM.push_back(_pythia->info.alphaEM());

    _sigmaGen.push_back(_pythia->info.sigmaGen());
    _sigmaErr.push_back(_pythia->info.sigmaErr());

    Int_t nWeights = _pythia->info.nWeights();
    _nWeights.push_back(nWeights);

    vector<Double_t> weights = {};
    for(Int_t j = 0; j < nWeights; j++){
      weights.push_back(_pythia->info.weight(j));
    }
    _weights.push_back(weights);

  }

  void Generator::_FillHepMC3Events(){
    HepMC3::GenEvent* evt = new HepMC3::GenEvent();
    _converter->fill_next_event(*_pythia,evt);
    _events.push_back(evt);
  }

  void Generator::generate(Int_t nevents, Bool_t refresh){
    if(!_initialized){
      cout << "Error: PythiaGenerator::Generator is not yet initialized." << endl;
      return;
    }

    // For array mode, we need to clear the arrays if requested.
    if(refresh && _arrayMode) _ClearContainers();

    for(Int_t i = 0; i < nevents; i++){
      Bool_t stat = _pythia->next();

      // In array mode, we fill the arrays with the various particle/event attributes.
      if(_arrayMode) _FillArrays();
    }

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