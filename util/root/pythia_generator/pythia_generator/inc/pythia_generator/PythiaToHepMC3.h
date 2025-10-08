// This code is a modification of the "HepMC3.h" header file from the
// Pythia8 project. Some slight modifications have been made to turn
// things into pointers, to deal with how ROOT loads headers, and the
// definitions of functions have been moved into a corresponding source
// file. Otherwise,the underlying algorithms are the same. The original
// header comment is below. - J. T. Offermann

// HepMC3.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
//
// Author: HepMC 3 Collaboration, hepmc-dev@.cern.ch
// Based on the HepMC2 interface by Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch.
// Header file and function definitions for the Pythia8ToHepMC class,
// which converts a PYTHIA event record to the standard HepMC format.

#ifndef PythiaGenerator_HepMC3_H
#define PythiaGenerator_HepMC3_H

// Standard library includes
#include <string>

namespace Pythia8{
  class Pythia;
  class Settings;
  class Info;
  class Event;
}

namespace HepMC3{
  class GenEvent;
}

using namespace std;

namespace PythiaGenerator {
class Pythia8ToHepMC3 {

public:

  // Constructor and destructor
  Pythia8ToHepMC3(): m_internal_event_number(0), m_print_inconsistency(true),
    m_free_parton_warnings(true), m_crash_on_problem(false),
    m_convert_gluon_to_0(false), m_store_pdf(true), m_store_proc(true),
    m_store_xsec(true), m_store_weights(true) {}
  virtual ~Pythia8ToHepMC3() {}

  // The recommended method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Pythia& pythia, HepMC3::GenEvent* evt, int ievnum = -1 );
  bool fill_next_event( Pythia8::Pythia& pythia, HepMC3::GenEvent& evt){return fill_next_event( pythia, &evt);}

  // Alternative method to convert Pythia events into HepMC3 ones.
  bool fill_next_event( Pythia8::Event& pyev, HepMC3::GenEvent&evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0) {
    return fill_next_event(pyev, &evt, ievnum, pyinfo, pyset); }
  bool fill_next_event( Pythia8::Event& pyev, HepMC3::GenEvent* evt, int ievnum = -1,
    const Pythia8::Info* pyinfo = 0, Pythia8::Settings* pyset = 0);

  // Read out values for some switches.
  bool print_inconsistency()  const { return m_print_inconsistency; }
  bool free_parton_warnings() const { return m_free_parton_warnings; }
  bool crash_on_problem()     const { return m_crash_on_problem; }
  bool convert_gluon_to_0()   const { return m_convert_gluon_to_0; }
  bool store_pdf()            const { return m_store_pdf; }
  bool store_proc()           const { return m_store_proc; }
  bool store_xsec()           const { return m_store_xsec; }
  bool store_weights()        const { return m_store_weights; }

  // Set values for some switches.
  void set_print_inconsistency(bool b = true)  { m_print_inconsistency  = b; }
  void set_free_parton_warnings(bool b = true) { m_free_parton_warnings = b; }
  void set_crash_on_problem(bool b = false)    { m_crash_on_problem     = b; }
  void set_convert_gluon_to_0(bool b = false)  { m_convert_gluon_to_0   = b; }
  void set_store_pdf(bool b = true)            { m_store_pdf            = b; }
  void set_store_proc(bool b = true)           { m_store_proc           = b; }
  void set_store_xsec(bool b = true)           { m_store_xsec           = b; }
  void set_store_weights(bool b = true)        { m_store_weights        = b; }

private:

    // Try to send warning message to the logger if present, otherwise
  // send it to cout if print_inconsistency().
  bool warning(const Pythia8::Info * pyinfo, string loc,string message, string extraInfo = "");

  // Use of copy constructor is not allowed.
  Pythia8ToHepMC3( const Pythia8ToHepMC3& ) {}

  // Data members.
  int  m_internal_event_number;
  bool m_print_inconsistency, m_free_parton_warnings, m_crash_on_problem,
       m_convert_gluon_to_0, m_store_pdf, m_store_proc, m_store_xsec,
       m_store_weights;

};
} // end namespace PythiaGenerator

#endif // end PythiaGenerator_HepMC3_H
