#ifndef ROOT_FJUTILS
#define ROOT_FJUTILS

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
// #include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif

// standard lib includes
#include <vector>

// // forward declarations for Pythia8
namespace Pythia8{
  class Pythia;
  class Event;
  class Particle;
}

using namespace std;
namespace PythiaGenerator{

  class Generator{
    public:
      Generator();
      virtual ~Generator();

      void readString(TString string); // access to Pythia8's readString functionality
      void setQuiet();
      void init(); // access to Pythia8's initialization

      void Generate(Int_t nEvents = 1, Bool_t refresh=kTRUE); // generate events

      // Getters for entire particle property arrays, across events
      vector<vector<Int_t>> getParticleIDArray(){return _pid;};
      vector<vector<Int_t>> getStatusArray(){return _status;};
      vector<vector<Int_t>> getStatusHepMCArray(){return _statusHepMC;};

      vector<vector<vector<Double_t>>> getMomentumArray(){return _p;};
      vector<vector<Double_t>> getEnergyArray(){return getComponentArray(_p, 0);};
      vector<vector<Double_t>> getPxArray(){return getComponentArray(_p, 1);};
      vector<vector<Double_t>> getPyArray(){return getComponentArray(_p, 2);};
      vector<vector<Double_t>> getPzArray(){return getComponentArray(_p, 3);};

      vector<vector<vector<Double_t>>> getProdVertexArray(){return _vProd;};
      vector<vector<Double_t>> getProdTArray(){return getComponentArray(_vProd, 0);};
      vector<vector<Double_t>> getProdXArray(){return getComponentArray(_vProd, 1);};
      vector<vector<Double_t>> getProdYArray(){return getComponentArray(_vProd, 2);};
      vector<vector<Double_t>> getProdZArray(){return getComponentArray(_vProd, 3);};

      vector<vector<Bool_t>> getHasProdVertexArray(){return _hasVertex;};
      vector<vector<Double_t>> getMassArray(){return _mass;};
      vector<vector<Int_t>> getMother1Array(){return _mother1;};
      vector<vector<Int_t>> getMother2Array(){return _mother2;};
      vector<vector<vector<Int_t>>> getMotherArray(){return _mothers;};

      vector<vector<Int_t>> getDaughter1Array(){return _daughter1;};
      vector<vector<Int_t>> getDaughter2Array(){return _daughter2;};
      vector<vector<vector<Int_t>>> getDaughterArray(){return _daughters;};

      vector<vector<Int_t>> getColorArray(){return _col;};
      vector<vector<Int_t>> getAntiColorArray(){return _acol;};

      // Getters for particle property arrays corresponding with a single event
      vector<Int_t> getParticleID(Int_t i){return _pid.at(i);};
      vector<Int_t> getStatus(Int_t i){return _status.at(i);};
      vector<Int_t> getStatusHepMC(Int_t i){return _statusHepMC.at(i);};
      vector<vector<Double_t>> getMomentum(Int_t i){return _p.at(i);};
      vector<vector<Double_t>> getVertex(Int_t i){return _vProd.at(i);};
      vector<Bool_t> getHasVertex(Int_t i){return _hasVertex.at(i);};
      vector<Double_t> getMass(Int_t i){return _mass.at(i);};
      vector<Int_t> getMother1(Int_t i){return _mother1.at(i);};
      vector<Int_t> getMother2(Int_t i){return _mother2.at(i);};
      vector<vector<Int_t>> getMothers(Int_t i){return _mothers.at(i);};

      vector<Int_t> getDaughter1(Int_t i){return _daughter1.at(i);};
      vector<Int_t> getDaughter2(Int_t i){return _daughter2.at(i);};
      vector<vector<Int_t>> getDaughters(Int_t i){return _daughters.at(i);};

      vector<Int_t> getColor(Int_t i){return _col.at(i);};
      vector<Int_t> getAntiColor(Int_t i){return _acol.at(i);};

      // Getters for event-level information, for whole arrays across events
      vector<Int_t> getId1PdfArray(){return _id1Pdf;};
      vector<Int_t> getId2PdfArray(){return _id2Pdf;};
      vector<Double_t> getX1PdfArray(){return _x1Pdf;};
      vector<Double_t> getX2PdfArray(){return _x2Pdf;};
      vector<Double_t> getQFacArray(){return _QFac;};
      vector<Double_t> getQRenArray(){return _QRen;};
      vector<Int_t> getNMPIArray(){return _nMPI;};
      vector<Int_t> getCodeArray(){return _code;};
      vector<Double_t> getAlphaSArray(){return _alphaS;};
      vector<Double_t> getAlphaEMArray(){return _alphaEM;};
      vector<Double_t> getSigmaGenArray(){return _sigmaGen;};
      vector<Double_t> getSigmaErrArray(){return _sigmaErr;};
      vector<Int_t> getNWeightsArray(){return _nWeights;};
      vector<vector<Double_t>> getWeightsArray(){return _weights;};

      // Getters for event-level information, from a specific event.
      Int_t getId1Pdf(Int_t i){return _id1Pdf.at(i);};
      Int_t getId2Pdf(Int_t i){return _id2Pdf.at(i);};
      Double_t getXd1Pdf(Int_t i){return _x1Pdf.at(i);};
      Double_t getX2Pdf(Int_t i){return _x2Pdf.at(i);};
      Double_t getQFac(Int_t i){return _QFac.at(i);};
      Double_t getQRen(Int_t i){return _QRen.at(i);};
      Int_t getNMPI(Int_t i){return _nMPI.at(i);};
      Int_t getCode(Int_t i){return _code.at(i);};
      Double_t getAlphaS(Int_t i){return _alphaS.at(i);};
      Double_t getAlphaEM(Int_t i){return _alphaEM.at(i);};
      Double_t getSigmaGen(Int_t i){return _sigmaGen.at(i);};
      Double_t getSigmaErr(Int_t i){return _sigmaErr.at(i);};
      Int_t getNWeights(Int_t i){return _nWeights.at(i);};
      vector<Double_t> getWeights(Int_t i){return _weights.at(i);};


      ClassDef(Generator, 1);

    private:

      Bool_t _initialized = kFALSE;
      void _ClearParticleContainers();
      void _ClearContainers();

      // Underlying instance of Pythia8 generator
      Pythia8::Pythia* _pythia = 0;

      /*
       * Containers for storing particle-level information.
       * We use vectors of vectors, to represent jagged 2D arrays.
       * The 1st dim is the event index, the 2nd is the object index (particle/vertex).
       * We use zero-indexing -- and adjust the mother/daughter information accordingly.
       */
      vector<vector<Int_t>> _pid = {}; // Particle ID
      vector<vector<Int_t>> _status = {}; // Particle status (Pythia encoding)
      vector<vector<Int_t>> _statusHepMC = {}; // Particle status (HepMC3 encoding)
      vector<vector<vector<Double_t>>> _p = {}; // momentum (E, px, py, pz)
      vector<vector<vector<Double_t>>> _vProd = {}; // production vertex
      vector<vector<Bool_t>> _hasVertex = {}; // whether or not particle has a production vertex set
      vector<vector<Double_t>> _mass = {}; // Particle mass
      vector<vector<Int_t>> _mother1 = {}; // Particle mother 1
      vector<vector<Int_t>> _mother2 = {}; // Particle mother 2
      vector<vector<Int_t>> _daughter1 = {}; // Particle daughter 1
      vector<vector<Int_t>> _daughter2 = {}; // Particle daughter 2

      vector<vector<Int_t>> _col = {}; // color
      vector<vector<Int_t>> _acol = {}; // anti-color

      vector<vector<vector<Int_t>>> _mothers = {}; // Particle mothers (via motherList())
      vector<vector<vector<Int_t>>> _daughters = {}; // Particle daughters (via daughterList())

      // Containers for storing event-level information.
      vector<Int_t> _id1Pdf = {};
      vector<Int_t> _id2Pdf = {};
      vector<Double_t> _x1Pdf = {};
      vector<Double_t> _x2Pdf = {};
      vector<Double_t> _QFac = {};
      vector<Double_t> _QRen = {};
      vector<Int_t> _nMPI = {};
      vector<Int_t> _code = {};
      vector<Double_t> _alphaS = {};
      vector<Double_t> _alphaEM = {};

      vector<Double_t> _sigmaGen = {};
      vector<Double_t> _sigmaErr = {};
      vector<Int_t> _nWeights = {};
      vector<vector<Double_t>> _weights = {};

    // Methods
    vector<vector<Double_t>> getComponentArray(vector<vector<vector<Double_t>>> inputArray, Int_t index);
  };
}

#endif