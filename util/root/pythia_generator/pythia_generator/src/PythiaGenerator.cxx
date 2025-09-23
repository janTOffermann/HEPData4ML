#include <pythia_generator/PythiaGenerator.h>

//standard library includes
#include <algorithm> // std::transform

// Pythia8 includes
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

using namespace std;

namespace PythiaGenerator{

  Generator::Generator(){
    _pythia = new Pythia8::Pythia("",kFALSE); // avoid printing the banner
  }

  Generator::~Generator(){
    delete _pythia;
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
    return;
  }

  void Generator::Generate(Int_t nevents, Bool_t refresh){
    if(!_initialized){
      cout << "Error: PythiaGenerator::Generator is not yet initialized."
      return;
    }

    if(refresh) _ClearContainers();
    for(Int_t i = 0; i < nevents; i++){
      Bool_t stat = _pythia->next();

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
    return;
  }



}