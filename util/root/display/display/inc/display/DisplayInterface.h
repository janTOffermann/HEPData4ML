#ifndef ROOT_EVENT_DISPLAY
#define ROOT_EVENT_DISPLAY

// ROOT includes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#include "TString.h"
// #include "Math/Vector4D.h"
// #include "Math/VectorUtil.h"
#endif


// standard lib includes
#include <vector>

using namespace std;
namespace Display{

  class EventDisplay;

  class DisplayInterface{
    public:
      DisplayInterface(TString delphesCardFilepath);
      ~DisplayInterface();

      void SetMode(Int_t mode);
      // void SetColorMode(Int_t mode);

      void DisplayEvent(TString delphesCardFilepath);
      EventDisplay* GetEventDisplay(){return eventDisplay_;}; // easy way to access a bunch of functions



      void Test(TString filename);

    private:
      EventDisplay* eventDisplay_;
      TString delphesCardFilepath_;

      void TestData();



  };
}

#endif