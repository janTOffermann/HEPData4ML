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
namespace EventDisplay{

  class Display{
    public:
      Display();
      ~Display();


      void DisplayEvent(TString delphesCardFilepath);


    private:

  };
}

#endif