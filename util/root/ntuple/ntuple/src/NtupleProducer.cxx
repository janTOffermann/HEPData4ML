#include <ntuple/NtupleProducer.h>
#include <ntuple/ParticleSelections.h>

//standard library includes
#include <algorithm> // std::transform. std::find

// ROOT includes
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TMath.h"
// #include "TLeafF.h"
#include "TLeafElement.h"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/FourVector.h"

using namespace std;

namespace NtupleProducer{

  Converter::Converter(){
    _evt = new HepMC3::GenEvent();
  }


  Converter::~Converter(){
    if(_readerAscii != 0) delete _readerAscii;
    if(_readerRoot != 0) delete _readerRoot;
    if(_evt != 0) delete _evt;
    if(_delphesReader != 0) delete _delphesReader;

    if(_delphesFile != 0){
      _delphesFile->Close(); // deletes _delphesTree?
      delete _delphesFile;
    }

    if(_outputFile != 0){
      _outputFile->Close(); // TODO: Is this needed?
      delete _outputFile;
    }

  }


  void Converter::SetInputFilesHepMC(vector<TString> files){
    _inputFilesHepMC = files;

    // Prepare the corresponding output filenames
    _SetOutputFiles();

  };

  void Converter::AddInputFileHepMC(TString file){
    _inputFilesHepMC.push_back(file);
    _SetOutputFiles();
  };

  void Converter::AddInputFileDetector(TString file){
    _inputFilesDetector.push_back(file);
    _hasDetectorFiles = kTRUE;
  }

  void Converter::SetInputFilesDetector(vector<TString> files){
    _inputFilesDetector = files;
    if(_inputFilesDetector.size() > 0) _hasDetectorFiles = kTRUE;
  };

  void Converter::_SetOutputFiles(){
    _outputFiles = {};

    for(TString inputFile : _inputFilesHepMC){
      TString outputFile = inputFile.ReplaceAll("hepmc.root","root"); // TODO: This is a bit fragile
      _outputFiles.push_back(outputFile);
    }
    return;
  }

  void Converter::_OpenHepMC3File(TString filename){
    // Figure out if this is an ASCII or ROOT file.
    // We'll do the simple thing, and assume ROOT only if
    // the extension is ".root".
    _rootMode = kFALSE;
    TObjArray* arr = filename.Tokenize(".");
    TString extension = ((TObjString*)arr->At(arr->GetEntries() - 1))->String();
    if(extension.EqualTo("root")) _rootMode = kTRUE;

    if(_rootMode) _OpenHepMC3FileRoot(filename);
    else _OpenHepMC3FileAscii(filename);
    delete arr;
  }

  void Converter::_OpenHepMC3FileRoot(TString filename){

    if(_readerRoot != 0){
      _readerRoot->close();
      delete _readerRoot;
    }
    _readerRoot = new HepMC3::ReaderRootTree(filename.Data());
  }

  void Converter::_OpenHepMC3FileAscii(TString filename){

    if(_readerAscii != 0){
      _readerAscii->close();
      delete _readerAscii;
    }
    _readerAscii = new HepMC3::ReaderAscii(filename.Data());
  }

  Bool_t Converter::_ReadHepMCEvent(){
    if(_rootMode) _readerRoot->read_event(*_evt);
    else _readerAscii->read_event(*_evt);
    return !_failedHepMC();
  }

  Bool_t Converter::_failedHepMC(){
    if(_rootMode) return _readerRoot->failed();
    return _readerAscii->failed();
  }

  void Converter::_CreateHepMCBranches(){

    TString branchName;
    // Branches for the stable truth particles
    branchName = Form("%s.N", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.N,Form("%s/I",branchName.Data()));

    branchName = Form("%s.Pmu", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.momentum.pmu);

    branchName = Form("%s.Pmu_cyl", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.momentum.pmu_cyl);

    branchName = Form("%s.PdgId", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.pdgId);

    branchName = Form("%s.HepMC3Index", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.indexHepMC);

    branchName = Form("%s.Production.Xmu", _truthParticleBranchPrefix.Data());
    _outputTree->Branch(branchName,&_stableParticles.xmu_prod);

    // TODO: Handling of truth particle selections

    return;
  }

  Bool_t Converter::_CheckStringVector(vector<TString> v, TString target){
    auto it = std::find(v.begin(), v.end(), target);
    return it!=v.end();
  }

  void Converter::_DelphesMultiplicity(TString inputBranchName){
    TString branchName = Form("%s.N", inputBranchName.Data());
    _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.N,Form("%s/I",branchName.Data()));
    return;
  }

  void Converter::_DelphesMomentum(TString inputBranchName, vector<TString> attributes){
    if(_CheckStringVector(attributes,"ET") || _CheckStringVector(attributes,"PT")){
      // Connect it to the input Delphes branches.
      if(_CheckStringVector(attributes,"ET")){
        _delphesData[inputBranchName]->pt = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.ET", inputBranchName.Data()));
      }
      else{
        _delphesData[inputBranchName]->pt = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.PT", inputBranchName.Data()));
      }
      // Assuming there's also Eta and Phi. (Should be safe, given how Delphes data is structured)
      _delphesData[inputBranchName]->eta = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Eta", inputBranchName.Data()));
      _delphesData[inputBranchName]->phi = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Phi", inputBranchName.Data()));

      if(_CheckStringVector(attributes,"Mass")){
        _delphesData[inputBranchName]->mass = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Mass", inputBranchName.Data()));
      }

      // Connect it to the output tree branches.
      TString branchName;
      branchName = Form("%s.Pmu", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.momentum.pmu);
      branchName = Form("%s.Pmu_cyl", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.momentum.pmu_cyl);
    }
    return;
  }

  void Converter::_DelphesD0Z0(TString inputBranchName, vector<TString> attributes){
    if(_CheckStringVector(attributes,"D0")){

      // Connect it to the input Delphes branches.
      _delphesData[inputBranchName]->d0 = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.D0", inputBranchName.Data()));
      _delphesData[inputBranchName]->d0Error = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.ErrorD0", inputBranchName.Data()));
      _delphesData[inputBranchName]->z0 = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.DZ", inputBranchName.Data())); // Note: Why does Delphes call it "DZ" and not "Z0"?
      _delphesData[inputBranchName]->z0Error = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.ErrorDZ", inputBranchName.Data()));

      // Connect it to the output tree branches.
      TString branchName;

      if(!_delphesAddedN[branchName]){
        branchName = Form("%s.N", inputBranchName.Data());
        _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.N,Form("%s/I",branchName.Data()));
        _delphesAddedN[branchName] = kTRUE;
      }


      branchName = Form("%s.D0", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.trackData.d0);
      branchName = Form("%s.D0.Error", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.trackData.d0Error);
      branchName = Form("%s.Z0", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.trackData.z0);
      branchName = Form("%s.Z0.Error", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.trackData.z0Error);
    }
    return;
  }

  // Point of nearest approach to z axis
  void Converter::_DelphesXd(TString inputBranchName, vector<TString> attributes){

    // There are two conditions on which we'll write these Xd, Yd and Zd values to output:
    // 1) The input Delphes object has these fields (as above).
    // 2) The input Delphes object has D0, Z0 and Phi fields.
    _delphesFillXd[inputBranchName] = kFALSE;
    if(_CheckStringVector(attributes,"Xd")){

      // Connect it to the input Delphes branches.
      _delphesData[inputBranchName]->xd = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Xd", inputBranchName.Data()));
      _delphesData[inputBranchName]->yd = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Yd", inputBranchName.Data()));
      _delphesData[inputBranchName]->zd = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Zd", inputBranchName.Data()));

      _delphesFillXd[inputBranchName] = kTRUE;
      _delphesIsTrack[inputBranchName] = kTRUE; // TODO: Should this just be identical to _delphesFillXd[inputBranchName]?
    }
    else if(_CheckStringVector(attributes,"D0") && _CheckStringVector(attributes,"DZ") && _CheckStringVector(attributes,"Phi"))
    {
      _delphesFillXd[inputBranchName] = kTRUE;
    }

    if(_delphesFillXd[inputBranchName]){

      // Connect it to the output tree branches.
      TString branchName = Form("%s.Xdi", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.trackData.Xd);
    }
    return;
  }

  void Converter::_DelphesPdgIdCharge(TString inputBranchName, vector<TString> attributes){
    if(_CheckStringVector(attributes,"PID")){
      _delphesData[inputBranchName]->pdgId = std::make_unique<TTreeReaderArray<Int_t>>(*_delphesReader,Form("%s.PID", inputBranchName.Data()));
      TString branchName = Form("%s.PdgId", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.pdgId);
    }
    if(_CheckStringVector(attributes,"Charge")){
      _delphesData[inputBranchName]->charge = std::make_unique<TTreeReaderArray<Int_t>>(*_delphesReader,Form("%s.Charge", inputBranchName.Data()));
      TString branchName = Form("%s.Charge", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.charge);
    }
  }

  void Converter::_DelphesCalorimeter(TString inputBranchName, vector<TString> attributes){
    if(_CheckStringVector(attributes,"Eem")){
      _delphesData[inputBranchName]->Eem = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Eem", inputBranchName.Data()));
      TString branchName = Form("%s.E.EM", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.Eem);
    }
    if(_CheckStringVector(attributes,"Ehad")){
      _delphesData[inputBranchName]->Ehad = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Ehad", inputBranchName.Data()));
      TString branchName = Form("%s.E.Hadronic", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.Ehad);
    }
    if(_CheckStringVector(attributes,"Etrk")){
      _delphesData[inputBranchName]->Etrack = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Etrk", inputBranchName.Data()));
      TString branchName = Form("%s.E.Track", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.Etrack);
    }
    // We handle Edges specially -- the name passed to TTreeReaderArray must include [N] suffix, where N is fixed array size

    if(_CheckStringVector(attributes,"Edges")){
      _delphesData[inputBranchName]->hasEdges = kTRUE;
      TLeaf* edgesLeaf = _delphesTree->GetLeaf(Form("%s.Edges", inputBranchName.Data()));
      _delphesData[inputBranchName]->edgesSize = edgesLeaf->GetLen();  // typically 4
      TString branchNameFull = Form("%s.Edges[%i]",inputBranchName.Data(),_delphesData[inputBranchName]->edgesSize);
      _delphesTree->SetBranchAddress(branchNameFull,&_delphesData[inputBranchName]->Edges);

      TString branchName;

      // Special case for edgesLeaf->GetLen() == 2
      if(_delphesData[inputBranchName]->edgesSize == 2){
        branchName = Form("%s.Edges.Eta", inputBranchName.Data());
        _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.edgesEta);
      }
      else{
        branchName = Form("%s.Edges.Eta", inputBranchName.Data());
        _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.edgesEta);
        branchName = Form("%s.Edges.Phi", inputBranchName.Data());
        _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.caloData.edgesPhi);
      }
    }
  }

  void Converter::_DelphesPosition(TString inputBranchName, vector<TString> attributes){
    if(_delphesIsTrack[inputBranchName]) return; // don't turn on reading if this is determined to be a track-type object; these variables might just always be zeros there (not filled by Delphes)
    if(_CheckStringVector(attributes,"X")){ // going to assume X,Y,Z,T all available
      _delphesData[inputBranchName]->T = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.T", inputBranchName.Data()));
      _delphesData[inputBranchName]->X = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.X", inputBranchName.Data()));
      _delphesData[inputBranchName]->Y = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Y", inputBranchName.Data()));
      _delphesData[inputBranchName]->Z = std::make_unique<TTreeReaderArray<Float_t>>(*_delphesReader,Form("%s.Z", inputBranchName.Data()));

      TString branchName = Form("%s.Xmu", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->output.positionData);
    }
  }


  void Converter::_CreateDelphesBranch(TString inputBranchName){

    TString branchName;

    // Determine what leaves this branch has.
    vector<TString> availableLeaves = {};
    vector<TString> attributes = {};
    for(TString leafName : _delphesLeafNames){
      if(leafName.Contains(inputBranchName + ".")){
        availableLeaves.push_back(leafName);
        TObjArray* tokens = leafName.Tokenize(".");
        attributes.push_back(((TObjString*)tokens->At(tokens->GetEntries()-1))->String());
        delete tokens;
      }
    }

    /*
     * Now we roughly mimic the logic of util/reconstruction/conversion.py .
     * We use our DelphesReaderData struct to read whichever branches are available.
     */


    // Create the DelphesReaderData object.
    _delphesData[inputBranchName] = std::make_unique<DelphesReaderData>();

    _delphesAddedN[inputBranchName] = kFALSE; // we'll have multiple opportunities to add a multiplicity branch; we only do it once
    _delphesIsTrack[inputBranchName] = kFALSE;

    // 0) Object multiplicity
    _DelphesMultiplicity(inputBranchName);

    // 1) Handling momentum
    _DelphesMomentum(inputBranchName,attributes);

    // 2) Handling D0/Z0
    _DelphesD0Z0(inputBranchName,attributes);

    // 3) Handling Xdi
    _DelphesXd(inputBranchName,attributes);

    // 4) Handling pdgId and charge
    _DelphesPdgIdCharge(inputBranchName,attributes);

    // 5) Handling calorimeter information
    _DelphesCalorimeter(inputBranchName,attributes);

    // 6) Handling position information (for non-track objects)
    _DelphesPosition(inputBranchName,attributes);


    return;
  }

  void Converter::_CreateDelphesBranches(){

    if(!_hasDetectorFiles) return;

    _delphesLeafNames.clear();

    // For handing Delphes objects, we read in all the available branches and their leaves.
    TObjArray* branchList = _delphesTree->GetListOfBranches();

    // Now, determine which branches we'll actually handle, based on what was requested
    // and what is actually available.
    vector<TString> branchNames = {};
    for(Int_t i = 0; i < branchList->GetEntries(); i++){
      TString name = branchList->At(i)->GetName();

      auto it = std::find(branchNames.begin(), branchNames.end(), name);
      if (it != branchNames.end()) continue;

      for(TString reqName : _delphesObjectNames){
        if(name.EqualTo(reqName)){
          branchNames.push_back(name);
          break;
        }
      }
    }

    // Determine what leaves are available.
    TObjArray* leaveList = _delphesTree->GetListOfLeaves();
    for(Int_t i = 0; i < leaveList->GetEntries(); i++){
      // Leaf names incl. the branch name, e.g.: Particle_, Particle.fUniqueID, Particle.fBits, Particle.PID ...
      TString leafName = leaveList->At(i)->GetName();
      // cout << "Checking leaf " << leafName << endl;

      for (TString reqName : branchNames){
        if(leafName.Contains(reqName + ".")){
          // cout << "\tContains " << reqName << endl;
          _delphesLeafNames.push_back(leafName);
          break;
        }
      }
    }

    // Now loop over the Delphes objects. For each, we determine
    // *how* to add them, i.e. what branches we'll be extracting and how.
    for(TString branchName : branchNames){
      _CreateDelphesBranch(branchName);
    }
  }

  void Converter::_OpenDelphesFile(TString filename){
    if(filename.EqualTo("")){
      _hasDetectorFiles = kFALSE;
      return;
    }
    _hasDetectorFiles = kTRUE;

    // If we have previously opened a Delphes file, close it.
    if(_delphesFile != 0){
      _delphesFile->Close();
      delete _delphesFile; // TODO: Is this necessary?
    }

    // Since we haven't built against Delphes, ROOT will generate
    // warnings about not having dictionaries for Delphes classes.
    // This is OK since we're reading the leaves.
    Int_t oldIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    _delphesFile = new TFile(filename,"READ");
    TString treeName = "Delphes";
    _delphesTree = (TTree*)_delphesFile->Get(treeName);
    _delphesReader = new TTreeReader(_delphesTree);
    gErrorIgnoreLevel = oldIgnoreLevel;

    //DEBUG
    // _delphesTree->Print();
    return;
  }

  void Converter::_FillStableParticles(){
    // Can loop in this simple way because the "selection algorithm" is hardcoded:
    // we simply take all particles where the HepMC3 status is equal to 1.
    _stableParticles.Clear();
    Int_t N = 0;
    for(Size_t i = 0; i < _evt->particles().size(); i++){
      shared_ptr<HepMC3::GenParticle> par = _evt->particles().at(i);
      if(par->status() != 1) continue;
      HepMC3::FourVector momentum = par->momentum();

      _stableParticles.pdgId.push_back(par->pid());
      _stableParticles.indexHepMC.push_back(par->id());
      _stableParticles.momentum.pmu.push_back({momentum.e(),momentum.px(),momentum.py(),momentum.pz()}); // note use of double-braces (moved from vector<vector<Double_t>> to vector<FourVector>)
      _stableParticles.momentum.pmu_cyl.push_back({momentum.pt(),momentum.eta(),momentum.phi(),momentum.m()});

      HepMC3::FourVector production_vertex = par->production_vertex()->position();
      _stableParticles.xmu_prod.push_back({production_vertex.t(),production_vertex.x(),production_vertex.y(),production_vertex.z()});
      N++;
    }
    _stableParticles.N = N;
    return;
  }

  void Converter::_FillDelphesObjects(){
    ROOT::Math::PtEtaPhiMVector v;
    for(auto& [branchName, data] : _delphesData){
        data->Clear();

        // TODO: Package these chunks up into their own functions
        if(data->pt){
          Int_t N = data->pt->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            Double_t pt = (*data->pt)[i];
            Double_t eta = (*data->eta)[i];
            Double_t phi = (*data->phi)[i];
            Double_t m   = data->mass ? (*data->mass)[i] :
              _delphesMassDefault.find(branchName) != _delphesMassDefault.end() ? _delphesMassDefault[branchName] :
              0.0;
            v.SetCoordinates(pt,eta,phi,m);
            data->output.momentum.pmu.push_back({v.E(), v.Px(), v.Py(), v.Pz()});
            data->output.momentum.pmu_cyl.push_back({pt, eta, phi, m});
          }
        }
        if(data->d0){
          Int_t N = data->d0->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.trackData.d0.push_back((*data->d0)[i]);
            data->output.trackData.d0Error.push_back((*data->d0Error)[i]);
            data->output.trackData.z0.push_back((*data->z0)[i]);
            data->output.trackData.z0Error.push_back((*data->z0Error)[i]);

          }
        }

        if(_delphesFillXd[branchName]){ // checking this way since there are two conditions under which we fill
          if(data->xd){
            Int_t N = data->xd->GetSize();
            data->output.N = N;
            for(Int_t i = 0; i < N; i++){
              data->output.trackData.Xd.push_back({(*data->xd)[i], (*data->yd)[i], (*data->zd)[i]});
            }
          }
          else{ // earlier setup should guarantee that the necessary inputs exist
            Int_t N = data->d0->GetSize();
            data->output.N = N;
            for(Int_t i = 0; i < N; i++){
              Double_t d0 = (*data->d0)[i];
              Double_t z0 = (*data->z0)[i];
              Double_t phi = (*data->phi)[i];
              Double_t xd = d0 * TMath::Cos(phi);
              Double_t yd = d0 * TMath::Sin(phi);
              data->output.trackData.Xd.push_back({xd, yd, z0});
            }
          }
        }

        if(data->pdgId){
          Int_t N = data->pdgId->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.pdgId.push_back((*data->pdgId)[i]);
          }
        }
        if(data->charge){
          Int_t N = data->charge->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.charge.push_back((*data->charge)[i]);
          }
        }

        if(data->Eem){
          Int_t N = data->Eem->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.caloData.Eem.push_back((*data->Eem)[i]);
          }
        }

        if(data->Ehad){
          Int_t N = data->Ehad->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.caloData.Ehad.push_back((*data->Ehad)[i]);
          }
        }

        if(data->Etrack){
          Int_t N = data->Etrack->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            data->output.caloData.Etrack.push_back((*data->Etrack)[i]);
          }
        }

        if(data->hasEdges){

          // TODO: Could consider something even safer, but we should have picked up N from something above.
          Int_t N = data->output.N; // assuming we got this already

          Int_t idx = 0;

          for(Int_t i = 0; i < N; i++) {

            if(idx + data->edgesSize >= data->edgesMax){
              cout << Form("NtupleProducer::Converter::_FillDelphesObjects(): Warning, more edges for branch %s.Edges than allowed by buffer, truncating.",branchName.Data()) << endl;
              cout << Form("\t(edgesMax = %i, but there are %i objects for event %llu, with %i edges each)",data->edgesMax,N,_i,data->edgesSize);
              break;
            }

            if(data->edgesSize == 2){
              vector<Double_t> rapidityEdges = {(Double_t)data->Edges[idx], (Double_t)data->Edges[idx + 1]};
              data->output.caloData.edgesEta.push_back(rapidityEdges);
            }
            else{
              vector<Double_t> etaEdges = {};
              vector<Double_t> phiEdges = {};
              for(Int_t j = 0; j < data->edgesSize; j++){
                Double_t val = (Double_t)data->Edges[idx + j];
                if(j < data->edgesSize / 2) etaEdges.push_back(val);
                else phiEdges.push_back(val);
              }
              data->output.caloData.edgesEta.push_back(etaEdges);
              data->output.caloData.edgesPhi.push_back(phiEdges);
            }
            idx += data->edgesSize;
          }
        }
        if(data->X){
          Int_t N = data->X->GetSize();
          data->output.N = N;
          for(Int_t i = 0; i < N; i++){
            Double_t t = (*data->T)[i];
            Double_t x = (*data->X)[i];
            Double_t y = (*data->Y)[i];
            Double_t z = (*data->Z)[i];

            data->output.positionData.push_back({t, x, y, z});
          }
        }
    }
  }

  // iterate both reader and tree (need latter for some old-fashioned branch access)
  void Converter::_IterateDelphesTree(Int_t entry){
    _delphesReader->SetEntry(entry);
    _delphesTree->GetEntry(entry);

  }

  void Converter::AddTruthParticleSelector(TString selectionName, BaseSelector* selector){
    _truthParticleSelectors[selectionName] = std::unique_ptr<BaseSelector>(selector);
  }

  void Converter::Process(TString inputFileHepMC, TString inputFileDetector, TString outputFile){

    // Reset the counter (unsure if we'll end up needing this at all)
    _i = 0;

    // Open the HepMC file -- can be either ASCII or ROOT format.
    cout << "Opening HepMC..." << endl;
    _OpenHepMC3File(inputFileHepMC); // NOTE: We don't know how many events are in the file; in principle can figure it out for ROOT but not ASCII
    cout << "\tDone." << endl;
    // Open the detector file, if present.
    // NOTE: For now, we assume it is a Delphes file.
    cout << "Opening Delphes..." << endl;
    _OpenDelphesFile(inputFileDetector);
    cout << "\tDone." << endl;

    // Create the output file

    if(_outputFile!=0){
      _outputFile->Close();
      delete _outputFile;
    }

    _outputFile = new TFile(outputFile,"RECREATE");
    _outputTree = new TTree("hepdata4ml_tree","");

    // Create branches for the HepMC3 information.
    _CreateHepMCBranches();

    // Create branches for the detector information.
    _CreateDelphesBranches();

    // Loop over events
    while(_ReadHepMCEvent()){ // fills _evt, will break the loop when we reach the end of the file
      if(_hasDetectorFiles) _IterateDelphesTree(_i);

      // Fetch the stable truth particle data, place it in output buffers.
      _FillStableParticles();

      // Fetch the Delphes object data, place it in output buffers.
      _FillDelphesObjects();


      //...

      // Write event to TTree
      _outputTree->Fill();
      _i++;
    }
    //...
    _outputFile->cd();
    _outputTree->Write();
    _outputFile->Close();
    delete _outputFile;
    _outputFile = 0;
  }

}