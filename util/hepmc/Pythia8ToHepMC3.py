#=======================================
# This code is taken from the official
# HepMC3 repository, located at
# https://gitlab.cern.ch/hepmc/HepMC3.
# I have slightly modified it to improve
# the speed. - Jan T. Offermann
#=======================================
import sys
import numpy as np
import awkward as ak
import pythia8 as pyth
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from util.misc.timing import profile_method, profile_block

class PythiaToHepMC:
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True
        self.debug = False

        # NOTE: Making the hanging particle/free parton check optional,
        #       it is unclear (from existing comments) if it is needed
        #       and it takes a small but measurable amount of time. - Jan
        self.m_hanging_particle_check = False

        self.setup = HepMCSetup(hepmc_dir,verbose=False)
        # self.setup.PrepHepMC()
        python_dir = self.setup.GetPythonDirectory()

        # uncache_hepmc3()
        prepend_to_pythonpath(python_dir)

    def SetDebug(self,val:bool):
        self.debug = val

    # The recommended method to convert Pythia events into HepMC ones
    @profile_method('fill_next_event1')
    def fill_next_event1(self, pythia, evt, ievnum):
        return self.fill_next_event(pythia.event, evt, ievnum, pythia.infoPython(), pythia.settings)

    # Alternative method to convert Pythia events into HepMC ones
    def fill_next_event(self, pyev, evt, ievnum, pyinfo, pyset):
        from pyHepMC3 import HepMC3 as hm
        # 1. Error if no event passed.
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

        # Faster than the old loop method. - Jan
        hepevt_particles = [hm.GenParticle(
            hm.FourVector(prt.px(), prt.py(), prt.pz(), prt.e()),
            prt.id(), prt.statusHepMC()
        ) for prt in pyev]
        for i,particle in enumerate(hepevt_particles):
            particle.set_generated_mass(pyev[i].m())

        if(self.debug):
            for i,prt in enumerate(pyev):
                if(i == 0): continue # status 11 appears
                print('[{}]'.format(i-1),prt.status(), prt.id(), [x-1 for x in prt.motherList()], '->',prt.statusHepMC())
                print('\tdaughters = ',[x-1 for x in prt.daughterList()])


            debug_dict = {}
            for prt in pyev:
                if(prt.status() not in debug_dict.keys()):
                    debug_dict[prt.status()] = []
                debug_dict[prt.status()].append(prt.statusHepMC())
                debug_dict[prt.status()] = list(set(debug_dict[prt.status()]))

            for key,val in debug_dict.items():
                print('[{}] -> {}'.format(key,val))

        #        // 3. Fill vertex information and find beam particles.
        # For type compatibility
        vertex_cache = hm.GenEvent().vertices()
        beam_particles = hm.GenEvent().particles()

        for i,prt in enumerate(pyev):
            k = i - 1
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

        #        // Add particles and vertices in topological order
        if len(beam_particles) < 2:
            print("There are  ", len(beam_particles), "!=2 particles without mothers")
            if self.m_crash_on_problem:
                sys.exit(1)
        evt.add_tree(beam_particles)

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

class PythiaToHepMCBatch:
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True

        self.setup = HepMCSetup(hepmc_dir, verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

        self.code_converter = StatusCodeConverter()

    @profile_method('fill_batch_events_no_info')
    def fill_batch_events_no_info(self, awkward_batch, start_event_num=None):
        """
        Vectorized Pythia8->HepMC3 - process all events simultaneously.

        Args:
            awkward_batch: Awkward array containing batch of Pythia8 events
            start_event_num: Starting event number (optional)

        Returns:
            List of HepMC3 GenEvent objects
        """
        from pyHepMC3 import HepMC3 as hm

        # NOTE: need to trim off the first "particle" from each event, this is actually a
        # pseudo-particle with status=11, which represents "the event as a whole"

        # TODO: Code functions, but not quite as expected: batch_size = 1 always.
        batch_size = len(awkward_batch)

        # Set up event numbering
        if start_event_num is not None:
            event_numbers = list(range(start_event_num, start_event_num + batch_size))
        else:
            event_numbers = list(range(self.m_internal_event_number,
                                     self.m_internal_event_number + batch_size))
            self.m_internal_event_number += batch_size

        # Extract data for all events
        all_prt_data = awkward_batch['prt'][:,1:] # trimming off entry 0

        # Get flattened particle data across all events
        all_p = ak.flatten(all_prt_data['p'])
        all_px = np.asarray(all_p['px'])
        all_py = np.asarray(all_p['py'])
        all_pz = np.asarray(all_p['pz'])
        all_e = np.asarray(all_p['e'])

        # Flatten other particle properties
        all_mass = np.asarray(ak.flatten(all_prt_data['m']))
        all_pid = np.asarray(ak.flatten(all_prt_data['id']))
        all_status = np.asarray(ak.flatten(all_prt_data['status']))
        all_mother1 = np.asarray(ak.flatten(all_prt_data['mother1']))
        all_mother2 = np.asarray(ak.flatten(all_prt_data['mother2']))
        all_col = np.asarray(ak.flatten(all_prt_data['col']))
        all_acol = np.asarray(ak.flatten(all_prt_data['acol']))

        # Handle vProd with None values.
        all_vProd_flat = ak.flatten(all_prt_data['vProd'])
        all_x_prod_masked = ak.to_numpy(all_vProd_flat['px'], allow_missing=True)
        all_y_prod_masked = ak.to_numpy(all_vProd_flat['py'], allow_missing=True)
        all_z_prod_masked = ak.to_numpy(all_vProd_flat['pz'], allow_missing=True)
        all_t_prod_masked = ak.to_numpy(all_vProd_flat['e'],  allow_missing=True)

        all_x_prod = np.ma.filled(all_x_prod_masked, 0.0)
        all_y_prod = np.ma.filled(all_y_prod_masked, 0.0)
        all_z_prod = np.ma.filled(all_z_prod_masked, 0.0)
        all_t_prod = np.ma.filled(all_t_prod_masked, 0.0)

        # # TODO: TEMPORARY HACK!
        # # See: https://gitlab.com/Pythia8/releases/-/issues/634
        # # This hack probably doesn't actually fix things, the last vertex will be assigned the wrong info.
        # # But it may help the current unit tests pass, as a way of isolating any other issues.
        # all_x_prod = np.roll(all_x_prod,-1,axis=-1)
        # all_y_prod = np.roll(all_y_prod,-1,axis=-1)
        # all_z_prod = np.roll(all_z_prod,-1,axis=-1)
        # all_t_prod = np.roll(all_t_prod,-1,axis=-1)

        all_vProd_is_none = all_x_prod_masked.mask if hasattr(all_x_prod_masked, 'mask') else np.zeros(len(all_x_prod), dtype=bool)

        # Get event boundaries (number of particles per event)
        particles_per_event = ak.num(all_prt_data)
        event_starts = np.concatenate([[0], np.cumsum(particles_per_event)[:-1]])
        event_ends = np.cumsum(particles_per_event)

        # Produce the mother lists. These are used in two separate places.
        all_mother_lists = [
            self._build_mother_lists(
                all_pid[event_starts[i]:event_ends[i]],
                all_status[event_starts[i]:event_ends[i]],
                all_mother1[event_starts[i]:event_ends[i]],
                all_mother2[event_starts[i]:event_ends[i]],
                zero_index=True
            )
            for i in range(batch_size)
        ]

        # Pythia8->HepMC3 status conversion
        # TODO: Is there a more robust way to do it than this function?
        all_hepmc_status = self._pythia_to_hepmc_status(all_status, all_pid, all_mother_lists, particles_per_event)

        # Create hepmc3 particles
        total_particles = len(all_px)
        all_particles = []

        for i in range(total_particles):
            particle = hm.GenParticle(
                hm.FourVector(all_px[i], all_py[i], all_pz[i], all_e[i]),
                int(all_pid[i]),
                int(all_hepmc_status[i])
            )
            particle.set_generated_mass(all_mass[i])
            all_particles.append(particle)

        # Split particles back into individual events and build event structures
        hepmc_events = [hm.GenEvent() for i in range(batch_size)] # pre-allocate event list; might speed things up slightly
        for i in range(batch_size):
            event = hepmc_events[i]
            event.set_event_number(event_numbers[i])
            event.set_units(hm.Units.GEV, hm.Units.MM)

            # Get particles for this event
            start_idx = event_starts[i]
            end_idx = event_ends[i]
            event_particles = all_particles[start_idx:end_idx]

            mother_lists = all_mother_lists[i]
            # mother_lists = self._build_mother_lists(
            #     all_pid[start_idx:end_idx],
            #     all_status[start_idx:end_idx],
            #     all_mother1[start_idx:end_idx],
            #     all_mother2[start_idx:end_idx],
            #     zero_index=True
            # )

            # Build vertices for this event
            self._build_event_vertices(
                event, event_particles,
                mother_lists,
                all_pid[start_idx:end_idx],
                all_status[start_idx:end_idx],
                all_x_prod[start_idx:end_idx], all_y_prod[start_idx:end_idx],
                all_z_prod[start_idx:end_idx], all_t_prod[start_idx:end_idx],
                all_vProd_is_none[start_idx:end_idx]
            )

            # Set color attributes for this event
            self._set_color_attributes(
                event_particles,
                all_col[start_idx:end_idx], all_acol[start_idx:end_idx]
            )

        return hepmc_events

    @profile_method('fill_batch_events')
    def fill_batch_events(self, awkward_batch, python_info, start_event_num=None):
        """
        Enhanced version that also handles PDF info, cross-sections, etc.
        from the awkward batch structure.

        Args:
            awkward_batch: Batch of events with 'prt' and 'info' fields
            start_event_num: Starting event number (optional)
        """
        from pyHepMC3 import HepMC3 as hm

        hepmc_events = self.fill_batch_events_no_info(awkward_batch, start_event_num)

        # Add event info from the 'info' field
        for i, evt in enumerate(hepmc_events):
            info_data = awkward_batch[i]['info']
            self._add_event_info(evt, info_data, python_info)

        return hepmc_events

    def _pythia_to_hepmc_status(self, all_status, all_pid, all_mother_lists, particles_per_event):
        """
        Convert the Pythia8 status codes to HepMC3 format.
        """
        all_hepmc_status = np.zeros_like(all_status)

        current_idx = 0
        for i,event_size in enumerate(particles_per_event):
            end_idx = current_idx + event_size

            # Extract event-specific arrays
            event_status = all_status[current_idx:end_idx]
            event_pid = all_pid[current_idx:end_idx]
            mother_lists = all_mother_lists[i]

            # Convert status for this event
            event_hepmc_status = self.code_converter.Convert(
                event_status, event_pid, mother_lists
            )

            all_hepmc_status[current_idx:end_idx] = event_hepmc_status
            current_idx = end_idx

        return all_hepmc_status

    def _build_event_vertices(self, event, event_particles, mother_lists,
                                        pid,status,
                                       x_prod, y_prod, z_prod, t_prod, vProd_is_none):
        """
        Build vertices for a single event using vectorized data.
        """
        from pyHepMC3 import HepMC3 as hm

        n_particles = len(event_particles)
        vertex_cache = hm.GenEvent().vertices()
        beam_particles = hm.GenEvent().particles()

        for i in range(n_particles):
            # Reconstruct mother list from mother1/mother2 (adjust for global indexing)
            mothers_list = mother_lists[i]
            # print(i, 'mothers_list = ',mothers_list)

            if len(mothers_list) > 0:
                # Find or create production vertex
                prod_vtx = None

                with profile_block('PythiaToHepMCBatch._build_event_vertices - ML1'):
                    for mother_idx in mothers_list:
                        if event_particles[mother_idx].end_vertex():
                            prod_vtx = event_particles[mother_idx].end_vertex()
                            break

                if prod_vtx is None:
                    prod_vtx = hm.GenVertex()

                    vertex_cache.append(prod_vtx)

                    with profile_block('PythiaToHepMCBatch._build_event_vertices - ML2'):
                        for mother_idx in mothers_list:
                            prod_vtx.add_particle_in(event_particles[mother_idx])

                # Set vertex position if available
                if prod_vtx.position().is_zero():
                    if not vProd_is_none[i]:
                        # if not (x_prod[i] == 0 and y_prod[i] == 0 and z_prod[i] == 0 and t_prod[i] == 0):
                        prod_pos = hm.FourVector(x_prod[i], y_prod[i], z_prod[i], t_prod[i])
                        prod_vtx.set_position(prod_pos)

                # if prod_vtx.position().is_zero(): # NOTE: Restructured this w.r.t. original, only create prod_pos inside this if statment. Should be faster. - Jan
                #     prod_pos = hm.FourVector(prt.xProd(), prt.yProd(), prt.zProd(), prt.tProd())
                #     if(not prod_pos.is_zero()):
                #         prod_vtx.set_position(prod_pos)
                #         print(' > Set vertex position ', prod_vtx.position().t(), prod_vtx.position().x(), prod_vtx.position().y(), prod_vtx.position().z())


                prod_vtx.add_particle_out(event_particles[i])
            else:
                beam_particles.append(event_particles[i])

        # Add particles to event
        if len(beam_particles) < 2:
            print(f"Warning: {len(beam_particles)} != 2 beam particles")
            if self.m_crash_on_problem:
                return False

        event.add_tree(beam_particles)
        return True

    def _set_color_attributes(self, event_particles, col, acol):
        """
        Set color flow attributes for all particles in an event.
        """
        from pyHepMC3 import HepMC3 as hm

        for i, particle in enumerate(event_particles):
            # Determine color type from col/acol values
            if col[i] == 0 and acol[i] > 0:
                # Anti-triplet
                particle.add_attribute("flow1", hm.IntAttribute(0))
                particle.add_attribute("flow2", hm.IntAttribute(int(acol[i])))
            elif col[i] > 0 and acol[i] == 0:
                # Triplet
                particle.add_attribute("flow1", hm.IntAttribute(int(col[i])))
                particle.add_attribute("flow2", hm.IntAttribute(0))
            elif col[i] > 0 and acol[i] > 0:
                # Octet (gluon)
                particle.add_attribute("flow1", hm.IntAttribute(int(col[i])))
                particle.add_attribute("flow2", hm.IntAttribute(int(acol[i])))

    def _add_event_info(self, evt, info_data, pyinfo=None):
        """Add PDF, cross-section, and other event information from awkward array."""
        from pyHepMC3 import HepMC3 as hm

        # PDF information
        if self.m_store_pdf:

            id1pdf = int(info_data['id1'])
            id2pdf = int(info_data['id2'])

            if self.m_convert_gluon_to_0:
                if id1pdf == 21:
                    id1pdf = 0
                if id2pdf == 21:
                    id2pdf = 0

            pdfinfo = hm.GenPdfInfo()
            pdfinfo.set(id1pdf, id2pdf,
                        float(info_data['x1']), float(info_data['x2']),
                        float(np.sqrt(info_data['Q2Fac'])),
                        float(info_data['pdf1']), float(info_data['pdf2']))
            evt.set_pdf_info(pdfinfo)

        # Note: Some of these things are fetched from the "info_data" that
        #       is part of Pythia8's nextBatch() output, whereas some things
        #       come from pyinfo that is fetched after the full batch is
        #       generated. I default to the former when possible, in case
        #       these are variables that change event-by-event (some probably do!). - Jan
        if(self.m_store_proc and pyinfo is not None):
            evt.add_attribute("mpi", hm.IntAttribute(pyinfo.nMPI())) # <- I think this is constant across events
            evt.add_attribute("signal_process_id", hm.IntAttribute(pyinfo.code())) # <- I think this is constant across events
            evt.add_attribute("event_scale", hm.DoubleAttribute(float(np.sqrt(info_data['Q2Ren'])))) # NOTE: QRen -> Q2Ren
            evt.add_attribute("alphaQCD", hm.DoubleAttribute(float(info_data['alphaS'])))
            evt.add_attribute("alphaQED", hm.DoubleAttribute(float(info_data['alphaEM'])))

        # Cross-section - also fetched from pyinfo.
        # TODO: Is this OK with batch processing? - Jan
        if self.m_store_xsec:
            xsec = hm.GenCrossSection()
            xsec.set_cross_section(pyinfo.sigmaGen() * 1e9, pyinfo.sigmaErr() * 1e9)
            evt.set_cross_section(xsec)

        # Weights - this might need adjustment
        if self.m_store_weights and 'weights' in info_data.fields:
            evt.weights().clear()
            weights = info_data['weights']
            if hasattr(weights, '__iter__'):
                for weight in weights:
                    evt.weights().append(float(weight))
            else:
                evt.weights().append(float(weights))

    def _build_mother_lists(self, pid_array, status_array, mother1_array, mother2_array, zero_index=True):
        """
        For an event, produce a "motherList-style" array of mother
        indices of particles. Defaults to 0-indexing.
        Based on: https://gitlab.com/Pythia8/releases/-/blob/master/src/Event.cc#L189
        """
        mother_lists = [np.empty(0,dtype=int) for _ in pid_array]

        for i in range(len(mother_lists)):
            mother1 = mother1_array[i]
            mother2 = mother2_array[i]
            status = status_array[i]
            abs_status = np.abs(status)

            if(not zero_index):
                mom1 = mother1
                mom2 = mother2

            else:
                mom1 = mother1 - 1
                mom2 = mother2 - 1


            if(abs_status == 11 or abs_status == 12):
                continue
            elif(mother1 == 0 and mother2 == 0):
                if(not zero_index):
                    mother_lists[i] = np.zeros(1) # NOTE: When does this happen?
                else:
                    mother_lists[i] = np.full(1,-1)
            elif(mother2 == 0 or mother2 == mother1):
                mother_lists[i] = np.array([mom1])

            elif((abs_status > 80 and abs_status < 90) or (abs_status > 100 and abs_status < 107)):
                mother_lists[i] = np.arange(mom1, mom2 + 1)
            else:
                if(mom2 > mom1):
                    mother_lists[i] = np.array([mom1, mom2])
                else:
                    mother_lists[i] = np.array([mom2, mom1])
        return mother_lists

class StatusCodeConverter:
    """
    Converts Pythia8 status codes to HepMC3 ones.
    Useful if using Pythia8's awkward array output,
    where one doesn't have a pythia8.Particle object
    to call statusHepMC() on.
    NOTE: This is quite ad-hoc, based on testing with
          Pythia8 plus its source code.
    """
    def __init__(self):
        self.hadron_lookup = {}
        self.debug = False
        self.status_lookup = {}

    def SetDebug(self,val:bool):
        self.debug = val

    def Convert(self, status_array, pid_array, mother_lists):
        """
        Based on https://gitlab.com/Pythia8/releases/-/blob/master/src/Event.cc#L383.
        Some changes since we have mother info instead of daughter info
        directly available from the Pythia.nextBatch() output.
        """

        n_particles = len(status_array)
        hepmc_status = np.zeros(n_particles, dtype=int)


        # Build daughter lists from mother information
        daughters = [[] for _ in range(n_particles)]

        for i in range(n_particles):
            mother_list = mother_lists[i]

            for mother_idx in mother_list:
                daughters[mother_idx].append(i)
        daughters = [list(set(x)) for x in daughters]

        # Now apply the Pythia8 statusHepMC logic
        for i in range(n_particles):
            status = status_array[i]
            particle_id = pid_array[i]

            # Positive codes are final particles
            if status > 0:
                hepmc_status[i] = 1
                continue

            # Status -12 are beam particles
            if status == -12:
                hepmc_status[i] = 4
                continue

            # Hadrons, muons, taus that decay normally are status 2
            if self.is_hadron(particle_id) or abs(particle_id) == 13 or abs(particle_id) == 15:
                if len(daughters[i]) > 0:
                    # Check first daughter
                    first_daughter_idx = daughters[i][0]
                    # Particle should not decay into itself (e.g. Bose-Einstein)
                    if pid_array[first_daughter_idx] != particle_id:
                        daughter_status = abs(status_array[first_daughter_idx])
                        if 90 < daughter_status < 95:
                            hepmc_status[i] = 2

                            continue

            # Other acceptable negative codes as their positive counterpart
            if -200 <= status <= -11:
                hepmc_status[i] = -status
                continue

            # Unacceptable codes as 0 (though this might cause issues)
            hepmc_status[i] = 0
        return hepmc_status

    def is_hadron(self,pid):
        """
        Check if particle ID corresponds to a hadron.
        """

        try:
            return self.hadron_lookup[pid]
        except:
            particle_data = pyth.ParticleDataEntry(pid)
            self.hadron_lookup[pid] = particle_data.isHadron()
        return self.hadron_lookup[pid]

    def Print(self):
        for key,val in self.status_lookup.items():
            print('[{}] -> {}'.format(key,val))

class PythiaWrapperToHepMCBatch:
    """
    Similar to PythiaToHepMCBatch, but this class
    interfaces with our custom Pythia8 interface.
    This also provides batches of events as awkward arrays,
    but they have things like particles' motherList and
    daughterList pre-computed, which makes the conversion
    here simpler.

    NOTE: This is currently far slower than using nextBatch()
          plus the PythiaToHepMCBatch converter, but maybe it
          offers a possibility for a speedup if restructured?
          It has the potential advantage of not needing to
          compute motherList/daughterList on the Python side,
          the latter being necessary for correctly deducing
          the statusHepMC for the particles.
    """
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True

        self.setup = HepMCSetup(hepmc_dir, verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

        self.code_converter = StatusCodeConverter()

    @profile_method('fill_batch_events_no_info')
    def fill_batch_events_no_info(self, awkward_batch, start_event_num=None):
        """
        Vectorized Pythia8->HepMC3 - process all events simultaneously.

        Args:
            awkward_batch: Awkward array containing batch of Pythia8 events
            start_event_num: Starting event number (optional)

        Returns:
            List of HepMC3 GenEvent objects
        """
        from pyHepMC3 import HepMC3 as hm

        # NOTE: need to trim off the first "particle" from each event, this is actually a
        # pseudo-particle with status=11, which represents "the event as a whole"

        # TODO: Code functions, but not quite as expected: batch_size = 1 always.
        batch_size = len(awkward_batch)

        # Set up event numbering
        if start_event_num is not None:
            event_numbers = list(range(start_event_num, start_event_num + batch_size))
        else:
            event_numbers = list(range(self.m_internal_event_number,
                                     self.m_internal_event_number + batch_size))
            self.m_internal_event_number += batch_size

        # Extract data for all events
        all_prt_data = awkward_batch['prt']

        with profile_block('fill_batch_events_no_info: A'):
            # Get flattened particle data across all events
            # all_p = ak.flatten(all_prt_data['p'])
            all_px = np.asarray(ak.flatten(all_prt_data['p']['px']))
            all_py = np.asarray(ak.flatten(all_prt_data['p']['py']))
            all_pz = np.asarray(ak.flatten(all_prt_data['p']['pz']))
            all_e = np.asarray(ak.flatten(all_prt_data['p']['e']))

        with profile_block('fill_batch_events_no_info: B'):
        # Flatten other particle properties
            all_mass = np.asarray(ak.flatten(all_prt_data['m']))
            all_pid = np.asarray(ak.flatten(all_prt_data['id']))
            all_hepmc_status = np.asarray(ak.flatten(all_prt_data['status'])) # HepMC status by default
            # all_mother1 = np.asarray(ak.flatten(all_prt_data['mother1']))
            # all_mother2 = np.asarray(ak.flatten(all_prt_data['mother2']))
            all_col = np.asarray(ak.flatten(all_prt_data['col']))
            all_acol = np.asarray(ak.flatten(all_prt_data['acol']))

        with profile_block('fill_batch_events_no_info: C'):
        # all_vProd_flat = ak.flatten(all_prt_data['vProd'])
            all_x_prod = ak.to_numpy(ak.flatten(all_prt_data['vProd']['x']))
            all_y_prod = ak.to_numpy(ak.flatten(all_prt_data['vProd']['y']))
            all_z_prod = ak.to_numpy(ak.flatten(all_prt_data['vProd']['z']))
            all_t_prod = ak.to_numpy(ak.flatten(all_prt_data['vProd']['t']))
            all_vProd_is_none = ak.flatten(all_prt_data['vProdStatus'])

        with profile_block('fill_batch_events_no_info: D'):
            all_mother_lists = all_prt_data['motherList'] # keeping as awkward array since its jagged

        with profile_block('fill_batch_events_no_info: E'):
            # Get event boundaries (number of particles per event)
            particles_per_event = [len(x['id']) for x in all_prt_data]
            # particles_per_event = ak.num(all_prt_data)
            event_starts = np.concatenate([[0], np.cumsum(particles_per_event)[:-1]])
            event_ends = np.cumsum(particles_per_event)

        with profile_block('fill_batch_events_no_info: F'):
            # Create hepmc3 particles
            total_particles = len(all_px)
            all_particles = []

            for i in range(total_particles):
                particle = hm.GenParticle(
                    hm.FourVector(all_px[i], all_py[i], all_pz[i], all_e[i]),
                    int(all_pid[i]),
                    int(all_hepmc_status[i])
                )
                particle.set_generated_mass(all_mass[i])
                all_particles.append(particle)

        with profile_block('fill_batch_events_no_info: G'):
            # Split particles back into individual events and build event structures
            hepmc_events = [hm.GenEvent() for i in range(batch_size)] # pre-allocate event list; might speed things up slightly
            for i in range(batch_size):

                with profile_block('fill_batch_events_no_info: G1'):


                    event = hepmc_events[i]
                    event.set_event_number(event_numbers[i])
                    event.set_units(hm.Units.GEV, hm.Units.MM)

                    # Get particles for this event
                    start_idx = event_starts[i]
                    end_idx = event_ends[i]
                    event_particles = all_particles[start_idx:end_idx]

                    mother_lists = all_mother_lists[i]
                    # mother_lists = self._build_mother_lists(
                    #     all_pid[start_idx:end_idx],
                    #     all_status[start_idx:end_idx],
                    #     all_mother1[start_idx:end_idx],
                    #     all_mother2[start_idx:end_idx],
                    #     zero_index=True
                    # )
                with profile_block('fill_batch_events_no_info: G2'):

                    # Build vertices for this event
                    self._build_event_vertices(
                        event, event_particles,
                        mother_lists,
                        all_x_prod[start_idx:end_idx], all_y_prod[start_idx:end_idx],
                        all_z_prod[start_idx:end_idx], all_t_prod[start_idx:end_idx],
                        all_vProd_is_none[start_idx:end_idx]
                    )
                with profile_block('fill_batch_events_no_info: G3'):

                    # Set color attributes for this event
                    self._set_color_attributes(
                        event_particles,
                        all_col[start_idx:end_idx], all_acol[start_idx:end_idx]
                    )

        return hepmc_events

    @profile_method('fill_batch_events')
    def fill_batch_events(self, awkward_batch, start_event_num=None):
        """
        Enhanced version that also handles PDF info, cross-sections, etc.
        from the awkward batch structure.

        Args:
            awkward_batch: Batch of events with 'prt' and 'info' fields
            start_event_num: Starting event number (optional)
        """
        from pyHepMC3 import HepMC3 as hm

        hepmc_events = self.fill_batch_events_no_info(awkward_batch, start_event_num)

        # Add event info from the 'info' field
        for i, evt in enumerate(hepmc_events):
            info_data = awkward_batch[i]['info']
            self._add_event_info(evt, info_data)

        return hepmc_events

    def _build_event_vertices(self, event, event_particles, mother_lists,
                                       x_prod, y_prod, z_prod, t_prod, vProd_is_none):
        """
        Build vertices for a single event using vectorized data.
        """
        from pyHepMC3 import HepMC3 as hm

        n_particles = len(event_particles)
        vertex_cache = hm.GenEvent().vertices()
        beam_particles = hm.GenEvent().particles()

        for i in range(n_particles):
            # Reconstruct mother list from mother1/mother2 (adjust for global indexing)
            with profile_block('PythiaWrapperToHepMCBatch._build_event_vertices - Block1'):
                mothers_list = mother_lists[i]
            # print(i, 'mothers_list = ',mothers_list)
            if len(mothers_list) > 0:
                with profile_block('PythiaWrapperToHepMCBatch._build_event_vertices - Block2'):

                    # Find or create production vertex
                    prod_vtx = None

                    for mother_idx in mothers_list:
                        if event_particles[mother_idx].end_vertex():
                            prod_vtx = event_particles[mother_idx].end_vertex()
                            break

                with profile_block('PythiaWrapperToHepMCBatch._build_event_vertices - Block3'):

                    if prod_vtx is None:
                        prod_vtx = hm.GenVertex()

                        vertex_cache.append(prod_vtx)
                        with profile_block('PythiaWrapperToHepMCBatch._build_event_vertices - ML2'):
                            for mother_idx in mothers_list:
                                prod_vtx.add_particle_in(event_particles[mother_idx])

                with profile_block('PythiaWrapperToHepMCBatch._build_event_vertices - Block4'):

                    # Set vertex position if available
                    if prod_vtx.position().is_zero():
                        if not vProd_is_none[i]:
                            prod_pos = hm.FourVector(x_prod[i], y_prod[i], z_prod[i], t_prod[i])
                            prod_vtx.set_position(prod_pos)

                    prod_vtx.add_particle_out(event_particles[i])
            else:
                beam_particles.append(event_particles[i])

        # Add particles to event
        if len(beam_particles) < 2:
            print(f"Warning: {len(beam_particles)} != 2 beam particles")
            if self.m_crash_on_problem:
                return False

        event.add_tree(beam_particles)
        return True

    def _set_color_attributes(self, event_particles, col, acol):
        """
        Set color flow attributes for all particles in an event.
        """
        from pyHepMC3 import HepMC3 as hm

        for i, particle in enumerate(event_particles):
            # Determine color type from col/acol values
            if col[i] == 0 and acol[i] > 0:
                # Anti-triplet
                particle.add_attribute("flow1", hm.IntAttribute(0))
                particle.add_attribute("flow2", hm.IntAttribute(int(acol[i])))
            elif col[i] > 0 and acol[i] == 0:
                # Triplet
                particle.add_attribute("flow1", hm.IntAttribute(int(col[i])))
                particle.add_attribute("flow2", hm.IntAttribute(0))
            elif col[i] > 0 and acol[i] > 0:
                # Octet (gluon)
                particle.add_attribute("flow1", hm.IntAttribute(int(col[i])))
                particle.add_attribute("flow2", hm.IntAttribute(int(acol[i])))

    def _add_event_info(self, evt, info_data):
        """Add PDF, cross-section, and other event information from awkward array."""
        from pyHepMC3 import HepMC3 as hm

        # PDF information
        if self.m_store_pdf:

            id1pdf = int(info_data['id1'])
            id2pdf = int(info_data['id2'])

            if self.m_convert_gluon_to_0:
                if id1pdf == 21:
                    id1pdf = 0
                if id2pdf == 21:
                    id2pdf = 0

            pdfinfo = hm.GenPdfInfo()
            pdfinfo.set(id1pdf, id2pdf,
                        float(info_data['x1']), float(info_data['x2']),
                        float(info_data['QFac']),
                        float(info_data['pdf1']), float(info_data['pdf2']))
            evt.set_pdf_info(pdfinfo)

        if(self.m_store_proc):
            evt.add_attribute("mpi", hm.IntAttribute(info_data['nMPI'])) # <- I think this is constant across events
            evt.add_attribute("signal_process_id", hm.IntAttribute(info_data['code'])) # <- I think this is constant across events
            evt.add_attribute("event_scale", hm.DoubleAttribute(float(info_data['QRen']))) # NOTE: QRen -> Q2Ren
            evt.add_attribute("alphaQCD", hm.DoubleAttribute(float(info_data['alphaS'])))
            evt.add_attribute("alphaQED", hm.DoubleAttribute(float(info_data['alphaEM'])))

        # Cross-section
        if self.m_store_xsec:
            xsec = hm.GenCrossSection()
            xsec.set_cross_section(info_data['sigmaGen'] * 1e9, info_data['sigmaErr'] * 1e9)
            evt.set_cross_section(xsec)

        # Weights - this might need adjustment
        if self.m_store_weights and 'weights' in info_data.fields:
            evt.weights().clear()
            weights = info_data['weights']
            if hasattr(weights, '__iter__'):
                for weight in weights:
                    evt.weights().append(float(weight))
            else:
                evt.weights().append(float(weights))