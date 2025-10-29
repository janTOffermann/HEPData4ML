/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
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


// Standard library includes
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <signal.h>

// ROOT includes
#include "TApplication.h"
#include "TROOT.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"
#include "TString.h"

// DELPHES includes
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
// #include "classes/DelphesHepMC3Reader.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

// custom includes
#include "DelphesHepMC3ROOTReader.h"

using namespace std;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

// Helper function to check if a string is a valid integer -- used for input parsing.
Bool_t isInteger(const char* str) {
    if (str == nullptr || *str == '\0') return kFALSE;

    // Handle optional minus sign
    if (*str == '-') str++;

    // Must have at least one digit after optional minus
    if (*str == '\0') return kFALSE;

    // Check that all remaining characters are digits
    while (*str) {
        if (*str < '0' || *str > '9') return kFALSE;
        str++;
    }
    return kTRUE;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesHepMC3";
  stringstream message;
  TString configFile = "";
  TString outputFile = "";
  TString inputFile = "";
  TString inputFilePileup = "";
  TFile *outputTFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0, *branchWeight = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  DelphesHepMC3ROOTReader *reader = 0;
  Int_t maxEvents, skipEvents;
  Long64_t eventCounter;

  Int_t rngDefaultSeed = -1;
  Bool_t foundDefaultSeed = kFALSE;

  if(argc < 4){
    cout << " Usage: " << appName
      << " config_file"
      << " output_file"
      << " input_file"
      << " [pileup_file]"
      << " [rng_seed]"
    << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file - input file in HepMC/ROOT format," << endl;
    cout << " pileup_file - input pileup file in HepMC/ROOT format [optional]," << endl;
    cout << " rng_seed - Integer to seed the random number generator [optional]." << endl;
    return 1;
  }

  if(argc == 4){
    if(isInteger(argv[argc-1])){ // config_file, output_file, input_file, rng_seed
      rngDefaultSeed = std::atoi(argv[argc-1]);
      foundDefaultSeed = kTRUE;
      cout << "** Using random seed default: " << rngDefaultSeed << " (can still be overwritten by settings in detector card!)" << endl;;
    }
    else{ // config_file, output_file, input_file, pileup_file
      inputFilePileup = argv[argc-1];
    }
  }

  else{ // argc >= 5
      inputFilePileup = argv[argc - 2];
      rngDefaultSeed = std::atoi(argv[argc-1]);
      foundDefaultSeed = kTRUE;
      cout << "** Using random seed default: " << rngDefaultSeed << " (can still be overwritten by settings in detector card!)" << endl;;
  }

  // Set the non-optional arguments
  configFile = argv[1];
  outputFile = argv[2];
  inputFile = argv[3];

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputTFile = TFile::Open(outputFile, "CREATE");

    if(outputTFile == NULL)
    {
      message << "can't create output file " << outputFile;
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputTFile, "Delphes");

    branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());
    branchWeight = treeWriter->NewBranch("Weight", Weight::Class());

    confReader = new ExRootConfReader;
    confReader->ReadFile(configFile.Data());

    maxEvents = confReader->GetInt("::MaxEvents", 0);
    skipEvents = confReader->GetInt("::SkipEvents", 0);

    if(maxEvents < 0)
    {
      throw runtime_error("MaxEvents must be zero or positive");
    }

    if(skipEvents < 0)
    {
      throw runtime_error("SkipEvents must be zero or positive");
    }

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);
    if(foundDefaultSeed) modularDelphes->SetDefaultRngSeed(rngDefaultSeed);

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    reader = new DelphesHepMC3ROOTReader();

    modularDelphes->InitTask();

    // Also check if there's a pileup file.

    cout << "** Reading " << inputFile << endl;

    if(inputFile == NULL)
    {
      message << "can't open " << argv[3];
      throw runtime_error(message.str());
    }

    reader->SetInputFile(inputFile);
    reader->SetInputPileupFile(inputFilePileup);

    // Loop over all objects
    eventCounter = 0;
    treeWriter->Clear();
    modularDelphes->Clear();
    reader->Clear();
    readStopWatch.Start();
    while((maxEvents <= 0 || eventCounter - skipEvents < maxEvents) && reader->ReadEvent() && !interrupted)
    {

      if(reader->EventReady()) // making this redundant
      {
        ++eventCounter;

        readStopWatch.Stop();

        if(eventCounter > skipEvents)
        {

          // Analyze() calls FinalizeParticles(), which must be called before modularDelphes::ProcessTask()
          reader->Analyze(factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray); // TODO

          procStopWatch.Start();
          modularDelphes->ProcessTask();
          procStopWatch.Stop();

          reader->AnalyzeEvent(branchEvent, eventCounter, &readStopWatch, &procStopWatch);
          reader->AnalyzeWeight(branchWeight);

          treeWriter->Fill();

          treeWriter->Clear();
        }

        modularDelphes->Clear();
        reader->Clear();

        readStopWatch.Start();
      }
      // progressBar.Update(ftello(inputFile), eventCounter);
    }

    // fseek(inputFile, 0L, SEEK_END);
    // progressBar.Update(ftello(inputFile), eventCounter, kTRUE);
    // progressBar.Finish();

    // if(inputFile != stdin) fclose(inputFile);

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete reader;
    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputTFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputTFile) delete outputTFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
