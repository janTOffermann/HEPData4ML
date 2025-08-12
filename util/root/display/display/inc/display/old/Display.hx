#ifndef ROOT_EVENT_DISPLAY
#define ROOT_EVENT_DISPLAY

// #include "BranchElement.h"

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
#include "TClonesArray.h"
#endif

// standard lib includes
#include <vector>
#include <map>

// // Forward declaration of some Delphes classes;
// // this circumvents issues when this header is
// // loaded by ROOT's gInterpreter for usage in PyROOT,
// // as gInterpreter doesn't know about headers from
// // external libraries (and I don't know how to fix that). - Jan
// class DelphesEventDisplay;
// class DelphesBranchBase;
// class Delphes3DGeometry;

// Delphes includes
#include "display/DelphesEventDisplay.h"
#include "display/Delphes3DGeometry.h"
#include "display/DelphesBranchElement.h"

using namespace std;
namespace EventDisplay{

  class EventDisplay : public DelphesEventDisplay{
    public:
      EventDisplay();
      ~EventDisplay();

      // // Functions for explicitly adding containers.
      // void AddCaloData(TString name, std::vector<BranchBase *> elements);
      // void AddTrackData(TString name, std::vector<BranchBase *> elements);

    protected:

      void BuildDetector(Delphes3DGeometry &det3D);

      std::map<TString,TClonesArray*> clone_arrays_;

      std::vector<DelphesBranchBase *> elements_;
  };
}
#endif