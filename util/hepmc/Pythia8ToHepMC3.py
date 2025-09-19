#=======================================
# This code is taken from the official
# HepMC3 repository, located at
# https://gitlab.cern.ch/hepmc/HepMC3.
# I have slightly modified it to improve
# the speed. - Jan T. Offermann
#=======================================
import sys
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from util.misc.timing import profile_method, profile_block

class Pythia8ToHepMC3:
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True

        # NOTE: Making the hanging particle/free parton check optional,
        #       it is unclear (from existing comments) if it is needed
        #       and it takes a small but measurable amount of time. - Jan
        self.m_hanging_particle_check = False

        self.setup = HepMCSetup(hepmc_dir,verbose=False)
        # self.setup.PrepHepMC()
        python_dir = self.setup.GetPythonDirectory()

        # uncache_hepmc3()
        prepend_to_pythonpath(python_dir)

    # The recommended method to convert Pythia events into HepMC ones
    @profile_method('fill_next_event1')
    def fill_next_event1(self, pythia, evt, ievnum):
        return self.fill_next_event(pythia.event, evt, ievnum, pythia.infoPython(), pythia.settings)

    # Alternative method to convert Pythia events into HepMC ones
    def fill_next_event(self, pyev, evt, ievnum, pyinfo, pyset):
        from pyHepMC3 import HepMC3 as hm
        # 1. Error if no event passed.
        with profile_block('Step 1'):
            if evt is None:
                print("Pythia8ToHepMC3::fill_next_event error - passed null event.")
                return False
            # Event number counter.
            if ievnum >= 0:
                evt.set_event_number(ievnum)
                self.m_internal_event_number = ievnum
            else:
                evt.set_event_number(self.m_internal_event_number)
                self.m_internal_event_number = self.m_internal_event_number + 1
            evt.set_units(hm.Units.GEV, hm.Units.MM)
        #        // 2. Fill particle information
        with profile_block('Step 2'):

            # Faster than the old loop method. - Jan
            hepevt_particles = [hm.GenParticle(
                hm.FourVector(prt.px(), prt.py(), prt.pz(), prt.e()),
                prt.id(), prt.statusHepMC()
            ) for prt in pyev]
            for i,prt in enumerate(hepevt_particles):
                prt.set_generated_mass(pyev[i].m())

        #        // 3. Fill vertex information and find beam particles.
        # For type compatibility
        with profile_block('Step 3a'):
            vertex_cache = hm.GenEvent().vertices()
            beam_particles = hm.GenEvent().particles()
            for i,prt in enumerate(pyev):
            # for i in range(0, pyev.size()):
                mothers = prt.motherList()
                if len(mothers) != 0:
                    prod_vtx = hepevt_particles[mothers[0]].end_vertex()
                    if prod_vtx is None:
                        prod_vtx = hm.GenVertex()
                        vertex_cache.append(prod_vtx)
                        for j in range(0, len(mothers)):
                            prod_vtx.add_particle_in(hepevt_particles[mothers[j]])

                    if prod_vtx.position().is_zero(): # NOTE: Restructured this w.r.t. original, only create prod_pos inside this if statment. Should be faster. - Jan
                        prod_pos = hm.FourVector(prt.xProd(), prt.yProd(), prt.zProd(), prt.tProd())
                        if(not prod_pos.is_zero()):
                            prod_vtx.set_position(prod_pos)

                    #                // Update vertex position if necessary
                    # if (not prod_pos.is_zero()) and prod_vtx.position().is_zero():
                    #     prod_vtx.set_position(prod_pos)
                    prod_vtx.add_particle_out(hepevt_particles[i])
                else:
                    beam_particles.append(hepevt_particles[i])

        with profile_block('Step 3b'):
            #        // Add particles and vertices in topological order
            if len(beam_particles) < 2:
                print("There are  ", len(beam_particles), "!=2 particles without mothers")
                if self.m_crash_on_problem:
                    sys.exit(1)
            evt.add_tree(beam_particles)

        with profile_block('Step 3c'):
            #        //Attributes should be set after adding the particles to event
            for i,prt in enumerate(pyev): # NOTE: A bit more unweildy than the old code, but marginally faster. - Jan
                colType = prt.colType()
                if(colType == -1):
                    hepevt_particles[i].add_attribute("flow1", hm.IntAttribute(0))
                    hepevt_particles[i].add_attribute("flow2", hm.IntAttribute(prt.acol()))
                elif(colType == 1):
                    hepevt_particles[i].add_attribute("flow1", hm.IntAttribute(prt.col()))
                    hepevt_particles[i].add_attribute("flow2", hm.IntAttribute(0))
                elif(colType == 2):
                    hepevt_particles[i].add_attribute("flow1", hm.IntAttribute(prt.col()))
                    hepevt_particles[i].add_attribute("flow2", hm.IntAttribute(prt.acol()))

            #        // If hadronization switched on then no final coloured particles.
            if pyset == None:
                doHadr = self.m_free_parton_warnings and pyset.flag("HadronLevel:Hadronize")
            else:
                doHadr = pyset.flag("HadronLevel:all") and pyset.flag("HadronLevel:Hadronize")

        #        // 4. Check for particles which come from nowhere, i.e. are without
        #        // mothers or daughters. These need to be attached to a vertex, or else
        #        // they will never become part of the event.
        with profile_block('Step 4'):
            for i,particle in enumerate(hepevt_particles):
                if(i == 0): continue

                if(self.m_hanging_particle_check):
                    #            // Check for particles not added to the event
                    #            // NOTE: We have to check if this step makes any sense in HepMC event standard
                    if not particle:
                        print("hanging particle ", i)
                        prod_vtx = hm.GenVertex()
                        prod_vtx.add_particle_out(particle)
                        evt.add_vertex(prod_vtx)

                #            // Also check for free partons (= gluons and quarks; not diquarks?).
                if doHadr and self.m_free_parton_warnings:
                    if particle.pid() == 21 and (particle.end_vertex() is None):
                        print("gluon without end vertex ", i)
                        if self.m_crash_on_problem:
                            sys.exit(1)
                    if abs(particle.pid()) <= 6 and (particle.end_vertex() is None):
                        print("quark without end vertex ", i)
                        if self.m_crash_on_problem:
                            sys.exit(1)

        #        // 5. Store PDF, weight, cross section and other event information.
        #        // Flavours of incoming partons.
        with profile_block('Step 5'):
            if self.m_store_pdf and pyinfo is not None:
                id1pdf = pyinfo.id1pdf()
                id2pdf = pyinfo.id2pdf()
                if self.m_convert_gluon_to_0:
                    if id1pdf == 21:
                        id1pdf = 0
                    if id2pdf == 21:
                        id2pdf = 0
                pdfinfo = hm.GenPdfInfo()
                pdfinfo.set(id1pdf, id2pdf, pyinfo.x1pdf(), pyinfo.x2pdf(), pyinfo.QFac(), pyinfo.pdf1(), pyinfo.pdf2())
                #            // Store PDF information.
                evt.set_pdf_info(pdfinfo)

            #        // Store process code, scale, alpha_em, alpha_s.
            if self.m_store_proc and pyinfo is not None:
                evt.add_attribute("mpi", hm.IntAttribute(pyinfo.nMPI()))
                evt.add_attribute("signal_process_id", hm.IntAttribute(pyinfo.code()))
                evt.add_attribute("event_scale", hm.DoubleAttribute(pyinfo.QRen()))
                evt.add_attribute("alphaQCD", hm.DoubleAttribute(pyinfo.alphaS()))
                evt.add_attribute("alphaQED", hm.DoubleAttribute(pyinfo.alphaEM()))

            #        // Store cross-section information in pb.
            if self.m_store_xsec and pyinfo is not None:
                xsec = hm.GenCrossSection()
                xsec.set_cross_section(pyinfo.sigmaGen() * 1e9, pyinfo.sigmaErr() * 1e9)
                evt.set_cross_section(xsec)

            #        // Store event weights.
            if self.m_store_weights and pyinfo is not None:
                evt.weights().clear()
                for iweight in range(0, pyinfo.nWeights()):
                    evt.weights().append(pyinfo.weight(iweight))

        #        // Done.
        return True
