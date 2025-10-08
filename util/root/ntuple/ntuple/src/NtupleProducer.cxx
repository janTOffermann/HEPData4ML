#include <ntuple/NtupleProducer.h>

//standard library includes
#include <algorithm> // std::transform. std::find

// ROOT includes
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"

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
    if(_readerDelphes != 0) delete _readerDelphes;

    if(_fileDelphes != 0){
      _fileDelphes->Close(); // deletes _treeDelphes?
      delete _fileDelphes;
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
    // cout << "\nFor branch " << inputBranchName << ", we have attributes:" << endl;
    // for(TString attribute : attributes){
    //   cout << "\t" << attribute << endl;
    // }

    /*
     * Now we roughly mimic the logic of util/reconstruction/conversion.py .
     * We use our DelphesReaderData struct to read whichever branches are available.
     */


    // Create the DelphesReaderData object.
    _delphesData[inputBranchName] = std::make_unique<DelphesReaderData>();

    // 1) Handling momentum
    if(_CheckStringVector(attributes,"ET") || _CheckStringVector(attributes,"PT")){

      // Connect it to the input Delphes branches.
      if(_CheckStringVector(attributes,"ET")){
        _delphesData[inputBranchName]->pt = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.ET", inputBranchName.Data()));
      }
      else{
        _delphesData[inputBranchName]->pt = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.PT", inputBranchName.Data()));
      }
      // Assuming there's also Eta and Phi. (Should be safe, given how Delphes data is structured)
      _delphesData[inputBranchName]->eta = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.Eta", inputBranchName.Data()));
      _delphesData[inputBranchName]->phi = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.Phi", inputBranchName.Data()));

      if(_CheckStringVector(attributes,"Mass")){
        _delphesData[inputBranchName]->mass = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.Mass", inputBranchName.Data()));
      }

      // Connect it to the output tree branches.
      branchName = Form("%s.N", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->N,Form("%s/I",branchName.Data()));
      branchName = Form("%s.Pmu", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->outputMomentum.pmu);
      branchName = Form("%s.Pmu_cyl", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->outputMomentum.pmu_cyl);
    }

    // 2) Handling D0/Z0
    if(_CheckStringVector(attributes,"D0")){

        // Connect it to the input Delphes branches.
        _delphesData[inputBranchName]->d0 = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.D0", inputBranchName.Data()));
        _delphesData[inputBranchName]->d0Error = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.ErrorD0", inputBranchName.Data()));
        _delphesData[inputBranchName]->z0 = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.DZ", inputBranchName.Data())); // Note: Why does Delphes call it "DZ" and not "Z0"?
        _delphesData[inputBranchName]->z0Error = std::make_unique<TTreeReaderArray<Float_t>>(*_readerDelphes,Form("%s.ErrorDZ", inputBranchName.Data()));

      // Connect it to the output tree branches.
      branchName = Form("%s.D0", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->trackData.d0);
      branchName = Form("%s.D0.Error", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->trackData.d0Error);
      branchName = Form("%s.Z0", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->trackData.z0);
      branchName = Form("%s.Z0.Error", inputBranchName.Data());
      _outputTree->Branch(branchName,&_delphesData[inputBranchName]->trackData.z0Error);

    }


    return;
  }

  void Converter::_CreateDelphesBranches(){

    if(!_hasDetectorFiles) return;

    _delphesLeafNames.clear();

    // For handing Delphes objects, we read in all the available branches and their leaves.
    TObjArray* branchList = _treeDelphes->GetListOfBranches();

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
    TObjArray* leaveList = _treeDelphes->GetListOfLeaves();
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
    if(_fileDelphes != 0){
      _fileDelphes->Close();
      delete _fileDelphes; // TODO: Is this necessary?
    }

    // Since we haven't built against Delphes, ROOT will generate
    // warnings about not having dictionaries for Delphes classes.
    // This is OK since we're reading the leaves.
    Int_t oldIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    _fileDelphes = new TFile(filename,"READ");
    TString treeName = "Delphes";
    _treeDelphes = (TTree*)_fileDelphes->Get(treeName);
    _readerDelphes = new TTreeReader(_treeDelphes);
    gErrorIgnoreLevel = oldIgnoreLevel;

    //DEBUG
    _treeDelphes->Print();

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
      _stableParticles.momentum.pmu.push_back({{momentum.e(),momentum.px(),momentum.py(),momentum.pz()}}); // note use of double-braces (moved from vector<vector<Double_t>> to vector<FourVector>)
      _stableParticles.momentum.pmu_cyl.push_back({{momentum.pt(),momentum.eta(),momentum.phi(),momentum.m()}});

      HepMC3::FourVector production_vertex = par->production_vertex()->position();
      _stableParticles.xmu_prod.push_back({{production_vertex.t(),production_vertex.x(),production_vertex.y(),production_vertex.z()}});
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
          data->N = N;
          for(Int_t i = 0; i < N; i++){
            Double_t pt = (*data->pt)[i];
            Double_t eta = (*data->eta)[i];
            Double_t phi = (*data->phi)[i];
            Double_t m   = data->mass ? (*data->mass)[i] :
              _delphesMassDefault.find(branchName) != _delphesMassDefault.end() ? _delphesMassDefault[branchName] :
              0.0;
            v.SetCoordinates(pt,eta,phi,m);
            data->outputMomentum.pmu.push_back({{v.E(), v.Px(), v.Py(), v.Pz()}});
            data->outputMomentum.pmu_cyl.push_back({{pt, eta, phi, m}});
          }
        }
        if(data->d0){
          Int_t N = data->d0->GetSize();
          for(Int_t i = 0; i < N; i++){
            data->trackData.d0.push_back((*data->d0)[i]);
            data->trackData.d0Error.push_back((*data->d0Error)[i]);
            data->trackData.z0.push_back((*data->z0)[i]);
            data->trackData.z0Error.push_back((*data->z0Error)[i]);

          }
        }

        // etc.
    }
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
      if(_hasDetectorFiles) _readerDelphes->Next();

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