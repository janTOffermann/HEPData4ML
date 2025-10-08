// This code is a modification of the "HepMC3.h" header file from the
// Pythia8 project -- I've copied function definitions from that header
// into this source file, since I want to have this as a "standalone"
// class (to avoid situations where the Pythia8 library I'm using might
// have been linked against a *different* HepMC3 library than the one I'm
// using throughout HepData4ML.). Below is the original block comment
// from the header file. - J. T. Offermann

// HepMC3.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
//
// Author: HepMC 3 Collaboration, hepmc-dev@.cern.ch
// Based on the HepMC2 interface by Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch.
// Header file and function definitions for the Pythia8ToHepMC class,
// which converts a PYTHIA event record to the standard HepMC format.

#include <pythia_generator/PythiaToHepMC3.h>

// Standard library includes
#include <vector>

// Pythia8 includes
#include "Pythia8/Pythia.h"
#include "Pythia8/HIInfo.h"

//HepMC3 includes
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenHeavyIon.h"
#include "HepMC3/GenPdfInfo.h"

using namespace std;

namespace PythiaGenerator{

  bool Pythia8ToHepMC3::fill_next_event( Pythia8::Pythia& pythia, HepMC3::GenEvent* evt, int ievnum){
      return fill_next_event( pythia.event, evt, ievnum, &pythia.info, &pythia.settings);
    }

  bool Pythia8ToHepMC3::fill_next_event( Pythia8::Event& pyev, HepMC3::GenEvent* evt, int ievnum,
    const Pythia8::Info* pyinfo, Pythia8::Settings* pyset){

    // 1. Error if no event passed.
    if (evt == nullptr) return warning(pyinfo,
      "Pythia8ToHepMC::fill_next_event", "passed null event");

    // Event number counter.
    if ( ievnum >= 0 ) {
      evt->set_event_number(ievnum);
      m_internal_event_number = ievnum;
    }
    else {
      evt->set_event_number(m_internal_event_number);
      ++m_internal_event_number;
    }

    // Set units to be GeV and mm, to agree with Pythia ones.
    evt->set_units(HepMC3::Units::GEV,HepMC3::Units::MM);

    // 1a. If there is a HIInfo object fill info from that.
    if ( pyinfo && pyinfo->hiInfo ) {
      auto ion = make_shared<HepMC3::GenHeavyIon>();
      ion->Ncoll_hard = pyinfo->hiInfo->nCollND();
      ion->Ncoll = pyinfo->hiInfo->nCollTot();
      ion->Npart_proj = pyinfo->hiInfo->nAbsProj() +
                        pyinfo->hiInfo->nDiffProj();
      ion->Npart_targ = pyinfo->hiInfo->nAbsTarg() +
                        pyinfo->hiInfo->nDiffTarg();
      ion->impact_parameter = pyinfo->hiInfo->b();
      evt->set_heavy_ion(ion);
    }

    // 2. Fill particle information.
    vector<HepMC3::GenParticlePtr> hepevt_particles;
    hepevt_particles.reserve( pyev.size() );
    for(int i = 0; i < pyev.size(); ++i) {
      hepevt_particles.push_back( std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector( pyev[i].px(), pyev[i].py(), pyev[i].pz(), pyev[i].e() ),
        pyev[i].id(), pyev[i].statusHepMC() ) );
      hepevt_particles[i]->set_generated_mass( pyev[i].m() );
    }

    // 3. Fill vertex information.
    vector<HepMC3::GenVertexPtr> vertex_cache;
    vector<HepMC3::GenParticlePtr> beam_particles;
    for (int i = 1; i < pyev.size(); ++i) {
      vector<int> mothers = pyev[i].motherList();
      sort(mothers.begin(),mothers.end());
      for (;;) {
        if (!mothers.empty() && mothers.front() == 0)
          mothers.erase(mothers.begin());
        else break;
      }
      if (mothers.size()) {
        HepMC3::GenVertexPtr prod_vtx = hepevt_particles[mothers[0]]->end_vertex();
        if (!prod_vtx) {
          prod_vtx = make_shared<HepMC3::GenVertex>();
          vertex_cache.push_back(prod_vtx);
          for (unsigned int j = 0; j < mothers.size(); ++j)
            prod_vtx->add_particle_in( hepevt_particles[mothers[j]] );
        }
        HepMC3::FourVector prod_pos( pyev[i].xProd(), pyev[i].yProd(),pyev[i].zProd(),
          pyev[i].tProd() );

        // Update vertex position if necessary.
        if (!prod_pos.is_zero() && prod_vtx->position().is_zero())
          prod_vtx->set_position( prod_pos );
        prod_vtx->add_particle_out( hepevt_particles[i] );
      } else beam_particles.push_back(hepevt_particles[i]);
    }

    // Reserve memory for the event.
    evt->reserve( hepevt_particles.size(), vertex_cache.size() );

    // Add particles and vertices in topological order.
    evt->add_tree( beam_particles );

    // Attributes should be set after adding the particles to event.
    for (int i = 0; i < pyev.size(); ++i) {
      /* TODO: Set polarization */
      // Colour flow uses index 1 and 2.
      int colType = pyev[i].colType();
      if (colType ==  -1 ||colType ==  1 || colType == 2) {
        int flow1 = 0, flow2 = 0;
        if (colType ==  1 || colType == 2) flow1 = pyev[i].col();
        if (colType == -1 || colType == 2) flow2 = pyev[i].acol();
        hepevt_particles[i]->add_attribute("flow1",
          make_shared<HepMC3::IntAttribute>(flow1));
        hepevt_particles[i]->add_attribute("flow2",
          make_shared<HepMC3::IntAttribute>(flow2));
      }
    }

    // If hadronization switched on then no final coloured particles.
    bool doHadr = (pyset == 0) ? m_free_parton_warnings
      : pyset->flag("HadronLevel:all") && pyset->flag("HadronLevel:Hadronize");

    // 4. Check for particles which come from nowhere, i.e. are without
    // mothers or daughters. These need to be attached to a vertex, or else
    // they will never become part of the event.
    for (int i = 1; i < pyev.size(); ++i) {

      // Check for particles not added to the event.
      // NOTE: We have to check if this step makes any sense in
      // the HepMC event standard.
      if ( hepevt_particles[i] == nullptr ||
        !hepevt_particles[i]->in_event()) {
        warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
          "found orphan particle", "i = " + Pythia8::toString(i));
        HepMC3::GenVertexPtr prod_vtx = make_shared<HepMC3::GenVertex>();
        prod_vtx->add_particle_out( hepevt_particles[i] );
        evt->add_vertex(prod_vtx);
      }

      // Also check for free partons (= gluons and quarks; not diquarks?).
      if ( doHadr && m_free_parton_warnings ) {
        if ( hepevt_particles[i]->pid() == 21
           && hepevt_particles[i]->end_vertex() == nullptr ) {
          warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
            "found gluon without end vertex", "i = " + Pythia8::toString(i));
          if ( m_crash_on_problem ) exit(1);
        }
        if ( abs(hepevt_particles[i]->pid()) <= 6
          && hepevt_particles[i]->end_vertex() == nullptr ) {
          warning(pyinfo, "Pythia8ToHepMC::fill_next_event",
            "found quark without end vertex", "i = " + Pythia8::toString(i));
          if ( m_crash_on_problem ) exit(1);
        }
      }
    }

    // 5. Store PDF, weight, cross section and other event information.
    // Flavours of incoming partons.
    if (m_store_pdf && pyinfo != 0) {
      int id1pdf = pyinfo->id1pdf();
      int id2pdf = pyinfo->id2pdf();
      if ( m_convert_gluon_to_0 ) {
        if (id1pdf == 21) id1pdf = 0;
        if (id2pdf == 21) id2pdf = 0;
      }

      // Store PDF information.
      HepMC3::GenPdfInfoPtr pdfinfo = make_shared<HepMC3::GenPdfInfo>();
      pdfinfo->set(id1pdf, id2pdf, pyinfo->x1pdf(), pyinfo->x2pdf(),
        pyinfo->QFac(), pyinfo->pdf1(), pyinfo->pdf2() );
      evt->set_pdf_info( pdfinfo );
    }

    // Store process code, scale, alpha_em, alpha_s.
    if (m_store_proc && pyinfo != 0) {
      evt->add_attribute("signal_process_id",
        std::make_shared<HepMC3::IntAttribute>( pyinfo->code()));
      evt->add_attribute("mpi",
        std::make_shared<HepMC3::IntAttribute>( pyinfo->nMPI()));
      evt->add_attribute("event_scale",
        std::make_shared<HepMC3::DoubleAttribute>(pyinfo->QRen()));
      evt->add_attribute("alphaQCD",
        std::make_shared<HepMC3::DoubleAttribute>(pyinfo->alphaS()));
      evt->add_attribute("alphaQED",
        std::make_shared<HepMC3::DoubleAttribute>(pyinfo->alphaEM()));
    }

    // Store event weights.
    if (m_store_weights && pyinfo != 0) {
      evt->weights().clear();
      for (int iWeight = 0; iWeight < pyinfo->numberOfWeights(); ++iWeight)
        evt->weights().push_back(pyinfo->weightValueByIndex(iWeight));
    }

    // Store cross-section information in pb.
    if (m_store_xsec && pyinfo != 0) {
      // First set atribute to event, such that
      // HepMC3::GenCrossSection::set_cross_section knows how many weights the
      // event has and sets the number of cross sections accordingly.
      HepMC3::GenCrossSectionPtr xsec = make_shared<HepMC3::GenCrossSection>();
      evt->set_cross_section(xsec);
      xsec->set_cross_section( pyinfo->sigmaGen() * 1e9,
        pyinfo->sigmaErr() * 1e9);
      // If multiweights with possibly different xsec, overwrite central value
      vector<double> xsecVec = pyinfo->weightContainerPtr->getTotalXsec();
      if (xsecVec.size() > 0) {
        for (unsigned int iXsec = 0; iXsec < xsecVec.size(); ++iXsec) {
          xsec->set_xsec(iXsec, xsecVec[iXsec]*1e9);
        }
      }
    }

    // Done.
    return true;
  }

  bool Pythia8ToHepMC3::warning(const Pythia8::Info * pyinfo, string loc, string message, string extraInfo) {
    if ( pyinfo )
      pyinfo->loggerPtr->warningMsg(loc, message, extraInfo);
    else if ( print_inconsistency() )
      cout << "Warning in " << loc << ": " << message << extraInfo << endl;
    return false;
  }


}