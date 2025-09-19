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
        # with profile_block('Step 1'):
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
        # with profile_block('Step 2'):

        # Faster than the old loop method. - Jan
        hepevt_particles = [hm.GenParticle(
            hm.FourVector(prt.px(), prt.py(), prt.pz(), prt.e()),
            prt.id(), prt.statusHepMC()
        ) for prt in pyev]
        for i,prt in enumerate(hepevt_particles):
            prt.set_generated_mass(pyev[i].m())

        #        // 3. Fill vertex information and find beam particles.
        # For type compatibility
        # with profile_block('Step 3a'):
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

        # with profile_block('Step 3b'):
        #        // Add particles and vertices in topological order
        if len(beam_particles) < 2:
            print("There are  ", len(beam_particles), "!=2 particles without mothers")
            if self.m_crash_on_problem:
                sys.exit(1)
        evt.add_tree(beam_particles)

        # with profile_block('Step 3c'):
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
        # with profile_block('Step 4'):
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
        # with profile_block('Step 5'):
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























##################################################

class PythiaToHepMCBatchV2:
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True
        self.m_hanging_particle_check = False

        self.setup = HepMCSetup(hepmc_dir, verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

    @profile_method('_fill_batch_events_no_info')
    def _fill_batch_events_no_info(self, awkward_batch, start_event_num=None):
        """
        Convert a batch of Pythia8 events (as awkward array) to HepMC3 events.
        Optimized version with pre-extraction of batch data.
        
        Args:
            awkward_batch: Awkward array containing batch of Pythia8 events
            start_event_num: Starting event number (optional)
            
        Returns:
            List of HepMC3 GenEvent objects
        """
        from pyHepMC3 import HepMC3 as hm
        
        batch_size = len(awkward_batch)
        hepmc_events = []
        
        # Set up event numbering
        if start_event_num is not None:
            event_numbers = list(range(start_event_num, start_event_num + batch_size))
        else:
            event_numbers = list(range(self.m_internal_event_number, 
                                     self.m_internal_event_number + batch_size))
            self.m_internal_event_number += batch_size
        
        # Pre-extract ALL batch data at once to reduce awkward array overhead
        with profile_block('Pre-extract batch data'):
            batch_prt = awkward_batch['prt']
            
            # Extract momentum data for all events
            batch_p = batch_prt['p']
            batch_px = [np.asarray(batch_p['px'][i]) for i in range(batch_size)]
            batch_py = [np.asarray(batch_p['py'][i]) for i in range(batch_size)]
            batch_pz = [np.asarray(batch_p['pz'][i]) for i in range(batch_size)]
            batch_e = [np.asarray(batch_p['e'][i]) for i in range(batch_size)]
            
            # Extract particle properties for all events
            batch_mass = [np.asarray(batch_prt['m'][i]) for i in range(batch_size)]
            batch_pid = [np.asarray(batch_prt['id'][i]) for i in range(batch_size)]
            batch_status = [np.asarray(batch_prt['status'][i]) for i in range(batch_size)]
            
            # Extract mother relationships for all events
            batch_mother1 = [np.asarray(batch_prt['mother1'][i]) for i in range(batch_size)]
            batch_mother2 = [np.asarray(batch_prt['mother2'][i]) for i in range(batch_size)]
            
            # Extract color info for all events
            batch_col = [np.asarray(batch_prt['col'][i]) for i in range(batch_size)]
            batch_acol = [np.asarray(batch_prt['acol'][i]) for i in range(batch_size)]
            
            # Extract vertex data for all events (handle None values)
            batch_vProd = batch_prt['vProd']
            batch_vProd_extracted = []
            for i in range(batch_size):
                vProd_i = batch_vProd[i]
                x_prod_masked = ak.to_numpy(vProd_i['px'], allow_missing=True)
                y_prod_masked = ak.to_numpy(vProd_i['py'], allow_missing=True) 
                z_prod_masked = ak.to_numpy(vProd_i['pz'], allow_missing=True)
                t_prod_masked = ak.to_numpy(vProd_i['e'],  allow_missing=True)
                
                x_prod = np.ma.filled(x_prod_masked, 0.0)
                y_prod = np.ma.filled(y_prod_masked, 0.0)
                z_prod = np.ma.filled(z_prod_masked, 0.0)
                t_prod = np.ma.filled(t_prod_masked, 0.0)
                
                vProd_is_none = x_prod_masked.mask if hasattr(x_prod_masked, 'mask') else np.zeros(len(x_prod), dtype=bool)
                
                batch_vProd_extracted.append((x_prod, y_prod, z_prod, t_prod, vProd_is_none))
        
        # Process events with pre-extracted data
        with profile_block('Process pre-extracted events'):
            for i, event_num in enumerate(event_numbers):
                hepmc_event = hm.GenEvent()
                if self._fill_from_extracted_arrays(
                    hepmc_event, event_num,
                    batch_px[i], batch_py[i], batch_pz[i], batch_e[i],
                    batch_mass[i], batch_pid[i], batch_status[i],
                    batch_mother1[i], batch_mother2[i],
                    batch_col[i], batch_acol[i],
                    *batch_vProd_extracted[i]
                ):
                    hepmc_events.append(hepmc_event)
                else:
                    print(f"Failed to convert event {event_num}")
                    
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
        
        hepmc_events = self._fill_batch_events_no_info(awkward_batch, start_event_num)
        
        # Add event info from the 'info' field
        with profile_block('Add event info'):
            for i, evt in enumerate(hepmc_events):
                info_data = awkward_batch[i]['info']
                self._add_event_info(evt, info_data)
        
        return hepmc_events

    def _fill_from_extracted_arrays(self, evt, event_num, px, py, pz, e, mass, pid, pythia_status,
                                   mother1, mother2, col, acol, x_prod, y_prod, z_prod, t_prod, vProd_is_none):
        """
        Fill HepMC3 event from pre-extracted numpy arrays.
        This should be faster than the awkward array version.
        """
        from pyHepMC3 import HepMC3 as hm
        
        with profile_block('Event setup'):
            evt.set_event_number(event_num)
            evt.set_units(hm.Units.GEV, hm.Units.MM)
        
        n_particles = len(px)
        
        # Convert status codes vectorially 
        with profile_block('Convert status codes'):
            hepmc_status = np.array([self._pythia_to_hepmc_status(stat) for stat in pythia_status])
        
        # Create particles
        with profile_block('Create particles'):
            hepevt_particles = []
            for i in range(n_particles):
                particle = hm.GenParticle(
                    hm.FourVector(px[i], py[i], pz[i], e[i]),
                    int(pid[i]), 
                    int(hepmc_status[i])
                )
                particle.set_generated_mass(mass[i])
                hepevt_particles.append(particle)
        
        # Build vertex structure (same logic as before)
        with profile_block('Build vertices'):
            vertex_cache = hm.GenEvent().vertices()
            beam_particles = hm.GenEvent().particles()
            
            for i in range(n_particles):
                # Reconstruct mother list from mother1/mother2 (Pythia8 uses 1-based indexing)
                mothers_list = []
                if mother1[i] > 0:  # mother1 > 0 indicates valid mother
                    mothers_list.append(mother1[i] - 1)  # Convert to 0-based indexing
                if mother2[i] > 0 and mother2[i] != mother1[i]:  # mother2 different from mother1
                    mothers_list.append(mother2[i] - 1)
                # Handle range case: if mother2 > mother1, all particles in range are mothers
                elif mother2[i] > mother1[i]:
                    for m in range(mother1[i], mother2[i] + 1):
                        if m - 1 not in mothers_list:  # Avoid duplicates
                            mothers_list.append(m - 1)
                
                if len(mothers_list) > 0:
                    # Find or create production vertex
                    prod_vtx = None
                    for mother_idx in mothers_list:
                        if 0 <= mother_idx < len(hepevt_particles) and hepevt_particles[mother_idx].end_vertex():
                            prod_vtx = hepevt_particles[mother_idx].end_vertex()
                            break
                    
                    if prod_vtx is None:
                        prod_vtx = hm.GenVertex()
                        vertex_cache.append(prod_vtx)
                        for mother_idx in mothers_list:
                            if 0 <= mother_idx < len(hepevt_particles):
                                prod_vtx.add_particle_in(hepevt_particles[mother_idx])
                    
                    # Set vertex position if available and vertex position is currently zero
                    # This follows the logic from the original converter
                    if prod_vtx.position().is_zero():
                        # Only try to set position if vProd is not None for this particle
                        if not vProd_is_none[i]:
                            # Check if the production vertex is non-zero
                            if not (x_prod[i] == 0 and y_prod[i] == 0 and z_prod[i] == 0 and t_prod[i] == 0):
                                prod_pos = hm.FourVector(x_prod[i], y_prod[i], z_prod[i], t_prod[i])
                                prod_vtx.set_position(prod_pos)
                    
                    prod_vtx.add_particle_out(hepevt_particles[i])
                else:
                    beam_particles.append(hepevt_particles[i])
        
        # Add particles to event
        with profile_block('Add to event'):
            if len(beam_particles) < 2:
                print(f"Warning: {len(beam_particles)} != 2 beam particles")
                if self.m_crash_on_problem:
                    return False
            
            evt.add_tree(beam_particles)
        
        # Set particle attributes (color flow)
        with profile_block('Set attributes'):
            for i in range(n_particles):
                particle = hepevt_particles[i]
                
                # Determine color type from col/acol values
                if col[i] == 0 and acol[i] == 0:
                    # No color
                    pass
                elif col[i] == 0 and acol[i] > 0:
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
        
        return True
                
    def _pythia_to_hepmc_status(self, pythia_status):
        """
        Convert Pythia8 status code to HepMC3 status code.
        Based on the standard conversion used in Pythia8ToHepMC3.
        
        This is a simplified version - you might want to make it more comprehensive
        based on your specific needs.
        """
        # Common Pythia8 to HepMC3 status conversions:
        # Pythia8 -> HepMC3
        # -11, -12: incoming beam particles -> 4 (beam particle)
        # 11, 12: outgoing beam particles after hard interaction -> 4
        # 21-29: particles from hard process -> 3 (outgoing)
        # 31-39, 41-49, etc: particles from showering/hadronization -> 1 (final state)
        # 51-59: particles from hadron/tau decays -> 1 (final state) or 2 (intermediate)
        # 61-69: particles from beam remnants -> 1 (final state)
        # 71-79, 81-89: particles from multiple interactions -> varies
        
        if pythia_status in [-11, -12, 11, 12]:
            return 4  # Beam particle
        elif 21 <= pythia_status <= 29:
            return 3  # Outgoing from hard process
        elif pythia_status > 0 and (
            (31 <= pythia_status <= 39) or 
            (41 <= pythia_status <= 49) or
            (51 <= pythia_status <= 59) or
            (61 <= pythia_status <= 69) or
            (71 <= pythia_status <= 79) or
            (81 <= pythia_status <= 89)
        ):
            return 1  # Final state particle
        elif pythia_status < 0:
            return 2  # Intermediate state
        else:
            # Default mapping for other cases
            if pythia_status > 0:
                return 1  # Final state
            else:
                return 2  # Intermediate

    def _add_event_info(self, evt, info_data):
        """Add PDF, cross-section, and other event information from awkward array."""
        from pyHepMC3 import HepMC3 as hm
        
        # You'll need to check what fields are available in info_data
        # This is a template based on typical Pythia8 info content
        
        try:
            # PDF information - adjust field names as needed
            if self.m_store_pdf and 'pdf' in info_data.fields:
                pdf_info = info_data['pdf']
                id1pdf = int(pdf_info['id1'])
                id2pdf = int(pdf_info['id2'])
                
                if self.m_convert_gluon_to_0:
                    if id1pdf == 21:
                        id1pdf = 0
                    if id2pdf == 21:
                        id2pdf = 0
                
                pdfinfo = hm.GenPdfInfo()
                pdfinfo.set(id1pdf, id2pdf, 
                           float(pdf_info['x1']), float(pdf_info['x2']), 
                           float(pdf_info['QFac']), 
                           float(pdf_info['pdf1']), float(pdf_info['pdf2']))
                evt.set_pdf_info(pdfinfo)
            
            # Process information - adjust field names as needed  
            if self.m_store_proc:
                if 'nMPI' in info_data.fields:
                    evt.add_attribute("mpi", hm.IntAttribute(int(info_data['nMPI'])))
                if 'code' in info_data.fields:
                    evt.add_attribute("signal_process_id", hm.IntAttribute(int(info_data['code'])))
                if 'QRen' in info_data.fields:
                    evt.add_attribute("event_scale", hm.DoubleAttribute(float(info_data['QRen'])))
                if 'alphaS' in info_data.fields:
                    evt.add_attribute("alphaQCD", hm.DoubleAttribute(float(info_data['alphaS'])))
                if 'alphaEM' in info_data.fields:
                    evt.add_attribute("alphaQED", hm.DoubleAttribute(float(info_data['alphaEM'])))
            
            # Cross-section - adjust field names as needed
            if self.m_store_xsec and 'sigmaGen' in info_data.fields:
                xsec = hm.GenCrossSection()
                sigma_gen = float(info_data['sigmaGen'])
                sigma_err = float(info_data['sigmaErr']) if 'sigmaErr' in info_data.fields else 0.0
                xsec.set_cross_section(sigma_gen * 1e9, sigma_err * 1e9)  # Convert to pb
                evt.set_cross_section(xsec)
            
            # Weights - this might need adjustment based on how weights are stored
            if self.m_store_weights and 'weights' in info_data.fields:
                evt.weights().clear()
                weights = info_data['weights']
                if hasattr(weights, '__iter__'):
                    for weight in weights:
                        evt.weights().append(float(weight))
                else:
                    evt.weights().append(float(weights))
                    
        except Exception as e:
            print(f"Warning: Could not add event info: {e}")
            # Continue without crashing


###############################


# v1
class PythiaToHepMCBatchV1:
    """
    A class for converting batches of Pythia8 events, that
    are produced using the Pythia8 Python interface's
    nextBatch() function -- which produces awkward arrays.
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
        self.m_hanging_particle_check = False

        self.setup = HepMCSetup(hepmc_dir, verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

    @profile_method('_fill_batch_events_no_info')
    def _fill_batch_events_no_info(self, awkward_batch, start_event_num=None):
        """
        Convert a batch of Pythia8 events (as awkward array) to HepMC3 events.
        
        Args:
            awkward_batch: Awkward array containing batch of Pythia8 events
            start_event_num: Starting event number (optional)
            
        Returns:
            List of HepMC3 GenEvent objects
        """
        from pyHepMC3 import HepMC3 as hm
        
        batch_size = len(awkward_batch)
        hepmc_events = []
        
        # Set up event numbering
        if start_event_num is not None:
            event_numbers = list(range(start_event_num, start_event_num + batch_size))
        else:
            event_numbers = list(range(self.m_internal_event_number, 
                                     self.m_internal_event_number + batch_size))
            self.m_internal_event_number += batch_size
        
        with profile_block('Batch conversion'):
            for i, (pyev_awkward, event_num) in enumerate(zip(awkward_batch, event_numbers)):
                hepmc_event = hm.GenEvent()
                if self._fill_single_event(pyev_awkward, hepmc_event, event_num):
                    hepmc_events.append(hepmc_event)
                else:
                    print(f"Failed to convert event {event_num}")
                    
        return hepmc_events

    @profile_method('PythiaToHepMCBatch._fill_single_event')
    def _fill_single_event(self, pyev_awkward, evt, event_num):
        """
        Fill a single HepMC3 event from awkward array data.
        Works with the Pythia8 nextBatch() structure.
        """
        from pyHepMC3 import HepMC3 as hm
        
        with profile_block('Event setup'):
            evt.set_event_number(event_num)
            evt.set_units(hm.Units.GEV, hm.Units.MM)
        
        # Extract particle data from the actual awkward structure
        with profile_block('Extract particle data'):
            prt_data = pyev_awkward['prt']  # Particle data
            n_particles = len(prt_data)
            
            # Convert to numpy arrays for faster access
            # Momentum 4-vectors
            p = prt_data['p']
            px = np.array(p['px'])
            py = np.array(p['py'])
            pz = np.array(p['pz'])
            e = np.array(p['e'])
            
            # Particle properties
            mass = np.array(prt_data['m'])
            pid = np.array(prt_data['id'])
            pythia_status = np.array(prt_data['status'])
            
            # Production vertex (x,y,z,t format) - make sure to handle None values!
            vProd_data = prt_data['vProd']
            x_prod_masked = ak.to_numpy(vProd_data['px'], allow_missing=True)
            y_prod_masked = ak.to_numpy(vProd_data['py'], allow_missing=True) 
            z_prod_masked = ak.to_numpy(vProd_data['pz'], allow_missing=True)
            t_prod_masked = ak.to_numpy(vProd_data['e'],  allow_missing=True)

            x_prod = np.ma.filled(x_prod_masked, 0.0)
            y_prod = np.ma.filled(y_prod_masked, 0.0)
            z_prod = np.ma.filled(z_prod_masked, 0.0)
            t_prod = np.ma.filled(t_prod_masked, 0.0)

            # Track which entries were None (masked)
            vProd_is_none = x_prod_masked.mask if hasattr(x_prod_masked, 'mask') else np.zeros(n_particles, dtype=bool)

            # Mother relationships - need to reconstruct motherList from mother1/mother2
            mother1 = np.array(prt_data['mother1'])
            mother2 = np.array(prt_data['mother2'])
            
            # Color information
            col = np.array(prt_data['col'])
            acol = np.array(prt_data['acol'])
            
        # Convert Pythia8 status to HepMC3 status
        with profile_block('Convert status codes'):
            hepmc_status = np.array([self._pythia_to_hepmc_status(stat) for stat in pythia_status])
        
        # Create particles in batch
        with profile_block('Create particles'):
            hepevt_particles = []
            for i in range(n_particles):
                particle = hm.GenParticle(
                    hm.FourVector(px[i], py[i], pz[i], e[i]),
                    int(pid[i]), 
                    int(hepmc_status[i])
                )
                particle.set_generated_mass(mass[i])
                hepevt_particles.append(particle)
        
        # Build vertex structure
        with profile_block('Build vertices'):
            vertex_cache = hm.GenEvent().vertices()
            beam_particles = hm.GenEvent().particles()
            
            for i in range(n_particles):
                # Reconstruct mother list from mother1/mother2 (Pythia8 uses 1-based indexing)
                mothers_list = []
                if mother1[i] > 0:  # mother1 > 0 indicates valid mother
                    mothers_list.append(mother1[i] - 1)  # Convert to 0-based indexing
                if mother2[i] > 0 and mother2[i] != mother1[i]:  # mother2 different from mother1
                    mothers_list.append(mother2[i] - 1)
                # Handle range case: if mother2 > mother1, all particles in range are mothers
                elif mother2[i] > mother1[i]:
                    for m in range(mother1[i], mother2[i] + 1):
                        if m - 1 not in mothers_list:  # Avoid duplicates
                            mothers_list.append(m - 1)
                
                if len(mothers_list) > 0:
                    # Find or create production vertex
                    prod_vtx = None
                    for mother_idx in mothers_list:
                        if 0 <= mother_idx < len(hepevt_particles) and hepevt_particles[mother_idx].end_vertex():
                            prod_vtx = hepevt_particles[mother_idx].end_vertex()
                            break
                    
                    if prod_vtx is None:
                        prod_vtx = hm.GenVertex()
                        vertex_cache.append(prod_vtx)
                        for mother_idx in mothers_list:
                            if 0 <= mother_idx < len(hepevt_particles):
                                prod_vtx.add_particle_in(hepevt_particles[mother_idx])
                    
                    # Set vertex position if available and vertex position is currently zero
                    # This follows the logic from the original converter
                    if prod_vtx.position().is_zero():
                        # Only try to set position if vProd is not None for this particle
                        if not vProd_is_none[i]:
                            # Check if the production vertex is non-zero
                            if not (x_prod[i] == 0 and y_prod[i] == 0 and z_prod[i] == 0 and t_prod[i] == 0):
                                prod_pos = hm.FourVector(x_prod[i], y_prod[i], z_prod[i], t_prod[i])
                                prod_vtx.set_position(prod_pos)
                    
                    prod_vtx.add_particle_out(hepevt_particles[i])
                else:
                    beam_particles.append(hepevt_particles[i])
        
        # Add particles to event
        with profile_block('Add to event'):
            if len(beam_particles) < 2:
                print(f"Warning: {len(beam_particles)} != 2 beam particles")
                if self.m_crash_on_problem:
                    return False
            
            evt.add_tree(beam_particles)
        
        # Set particle attributes (color flow)
        with profile_block('Set attributes'):
            for i in range(n_particles):
                particle = hepevt_particles[i]
                
                # Determine color type from col/acol values
                if col[i] == 0 and acol[i] == 0:
                    # No color
                    pass
                elif col[i] == 0 and acol[i] > 0:
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
        
        return True

    def _pythia_to_hepmc_status(self, pythia_status):
        """
        Convert Pythia8 status code to HepMC3 status code.
        Based on the standard conversion used in Pythia8ToHepMC3.
        
        This is a simplified version - you might want to make it more comprehensive
        based on your specific needs.
        """
        # Common Pythia8 to HepMC3 status conversions:
        # Pythia8 -> HepMC3
        # -11, -12: incoming beam particles -> 4 (beam particle)
        # 11, 12: outgoing beam particles after hard interaction -> 4
        # 21-29: particles from hard process -> 3 (outgoing)
        # 31-39, 41-49, etc: particles from showering/hadronization -> 1 (final state)
        # 51-59: particles from hadron/tau decays -> 1 (final state) or 2 (intermediate)
        # 61-69: particles from beam remnants -> 1 (final state)
        # 71-79, 81-89: particles from multiple interactions -> varies
        
        if pythia_status in [-11, -12, 11, 12]:
            return 4  # Beam particle
        elif 21 <= pythia_status <= 29:
            return 3  # Outgoing from hard process
        elif pythia_status > 0 and (
            (31 <= pythia_status <= 39) or 
            (41 <= pythia_status <= 49) or
            (51 <= pythia_status <= 59) or
            (61 <= pythia_status <= 69) or
            (71 <= pythia_status <= 79) or
            (81 <= pythia_status <= 89)
        ):
            return 1  # Final state particle
        elif pythia_status < 0:
            return 2  # Intermediate state
        else:
            # Default mapping for other cases
            if pythia_status > 0:
                return 1  # Final state
            else:
                return 2  # Intermediate
                
    @profile_method('fill_batch_events_with_info')
    def fill_batch_events(self, awkward_batch, start_event_num=None):
        """
        Enhanced version that also handles PDF info, cross-sections, etc.
        from the awkward batch structure.
        
        Args:
            awkward_batch: Batch of events with 'prt' and 'info' fields
            start_event_num: Starting event number (optional)
        """
        from pyHepMC3 import HepMC3 as hm
        
        hepmc_events = self._fill_batch_events_no_info(awkward_batch, start_event_num)
        
        # Add event info from the 'info' field
        with profile_block('Add event info'):
            for i, evt in enumerate(hepmc_events):
                info_data = awkward_batch[i]['info']
                self._add_event_info(evt, info_data)
        
        return hepmc_events

    def _add_event_info(self, evt, info_data):
        """Add PDF, cross-section, and other event information from awkward array."""
        from pyHepMC3 import HepMC3 as hm
        
        # You'll need to check what fields are available in info_data
        # This is a template based on typical Pythia8 info content
        
        try:
            # PDF information - adjust field names as needed
            if self.m_store_pdf and 'pdf' in info_data.fields:
                pdf_info = info_data['pdf']
                id1pdf = int(pdf_info['id1'])
                id2pdf = int(pdf_info['id2'])
                
                if self.m_convert_gluon_to_0:
                    if id1pdf == 21:
                        id1pdf = 0
                    if id2pdf == 21:
                        id2pdf = 0
                
                pdfinfo = hm.GenPdfInfo()
                pdfinfo.set(id1pdf, id2pdf, 
                           float(pdf_info['x1']), float(pdf_info['x2']), 
                           float(pdf_info['QFac']), 
                           float(pdf_info['pdf1']), float(pdf_info['pdf2']))
                evt.set_pdf_info(pdfinfo)
            
            # Process information - adjust field names as needed  
            if self.m_store_proc:
                if 'nMPI' in info_data.fields:
                    evt.add_attribute("mpi", hm.IntAttribute(int(info_data['nMPI'])))
                if 'code' in info_data.fields:
                    evt.add_attribute("signal_process_id", hm.IntAttribute(int(info_data['code'])))
                if 'QRen' in info_data.fields:
                    evt.add_attribute("event_scale", hm.DoubleAttribute(float(info_data['QRen'])))
                if 'alphaS' in info_data.fields:
                    evt.add_attribute("alphaQCD", hm.DoubleAttribute(float(info_data['alphaS'])))
                if 'alphaEM' in info_data.fields:
                    evt.add_attribute("alphaQED", hm.DoubleAttribute(float(info_data['alphaEM'])))
            
            # Cross-section - adjust field names as needed
            if self.m_store_xsec and 'sigmaGen' in info_data.fields:
                xsec = hm.GenCrossSection()
                sigma_gen = float(info_data['sigmaGen'])
                sigma_err = float(info_data['sigmaErr']) if 'sigmaErr' in info_data.fields else 0.0
                xsec.set_cross_section(sigma_gen * 1e9, sigma_err * 1e9)  # Convert to pb
                evt.set_cross_section(xsec)
            
            # Weights - this might need adjustment based on how weights are stored
            if self.m_store_weights and 'weights' in info_data.fields:
                evt.weights().clear()
                weights = info_data['weights']
                if hasattr(weights, '__iter__'):
                    for weight in weights:
                        evt.weights().append(float(weight))
                else:
                    evt.weights().append(float(weights))
                    
        except Exception as e:
            print(f"Warning: Could not add event info: {e}")
            # Continue without crashing


######################################

class PythiaToHepMCBatchV3:
    def __init__(self, hepmc_dir=None):
        self.m_internal_event_number = 0
        self.m_free_parton_warnings = False
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True
        self.m_hanging_particle_check = False

        self.setup = HepMCSetup(hepmc_dir, verbose=False)
        python_dir = self.setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)

    @profile_method('_fill_batch_events_no_info')
    def _fill_batch_events_no_info(self, awkward_batch, start_event_num=None):
        """
        Vectorized Pythia8->HepMC3 - process all events simultaneously.
        
        Args:
            awkward_batch: Awkward array containing batch of Pythia8 events
            start_event_num: Starting event number (optional)
            
        Returns:
            List of HepMC3 GenEvent objects
        """
        from pyHepMC3 import HepMC3 as hm
        
        batch_size = len(awkward_batch)
        
        # Set up event numbering
        if start_event_num is not None:
            event_numbers = list(range(start_event_num, start_event_num + batch_size))
        else:
            event_numbers = list(range(self.m_internal_event_number, 
                                     self.m_internal_event_number + batch_size))
            self.m_internal_event_number += batch_size
        
        with profile_block('Extract all batch data'):
            # Extract data for all events
            all_prt_data = awkward_batch['prt']
            
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
            
            # Handle vProd with None values
            all_vProd_flat = ak.flatten(all_prt_data['vProd'])
            all_x_prod_masked = ak.to_numpy(all_vProd_flat['px'], allow_missing=True)
            all_y_prod_masked = ak.to_numpy(all_vProd_flat['py'], allow_missing=True)
            all_z_prod_masked = ak.to_numpy(all_vProd_flat['pz'], allow_missing=True)
            all_t_prod_masked = ak.to_numpy(all_vProd_flat['e'], allow_missing=True)
            
            all_x_prod = np.ma.filled(all_x_prod_masked, 0.0)
            all_y_prod = np.ma.filled(all_y_prod_masked, 0.0)
            all_z_prod = np.ma.filled(all_z_prod_masked, 0.0)
            all_t_prod = np.ma.filled(all_t_prod_masked, 0.0)
            all_vProd_is_none = all_x_prod_masked.mask if hasattr(all_x_prod_masked, 'mask') else np.zeros(len(all_x_prod), dtype=bool)
            
            # Get event boundaries (number of particles per event)
            particles_per_event = ak.num(all_prt_data)
            event_starts = np.concatenate([[0], np.cumsum(particles_per_event)[:-1]])
            event_ends = np.cumsum(particles_per_event)
        
        with profile_block('Convert all status codes'):
            # Pythia8->HepMC3 status conversion
            # TODO: Is there a more robust way to do it than this function?
            all_hepmc_status = self._pythia_to_hepmc_status(all_status)
        
        with profile_block('Create all particles'):
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
        
        with profile_block('Build all events'):
            # Split particles back into individual events and build event structures
            # hepmc_events = []
            hepmc_events = [hm.GenEvent() for i in range(batch_size)]
            for i in range(batch_size):
                event = hepmc_events[i] # event = hm.GenEvent()
                event.set_event_number(event_numbers[i])
                event.set_units(hm.Units.GEV, hm.Units.MM)
                
                # Get particles for this event
                start_idx = event_starts[i]
                end_idx = event_ends[i]
                event_particles = all_particles[start_idx:end_idx]
                
                # Build vertices for this event
                self._build_event_vertices(
                    event, event_particles,
                    all_mother1[start_idx:end_idx], all_mother2[start_idx:end_idx],
                    all_x_prod[start_idx:end_idx], all_y_prod[start_idx:end_idx],
                    all_z_prod[start_idx:end_idx], all_t_prod[start_idx:end_idx],
                    all_vProd_is_none[start_idx:end_idx],
                    start_idx  # offset for mother indexing
                )
                
                # Set color attributes for this event
                self._set_color_attributes(
                    event_particles,
                    all_col[start_idx:end_idx], all_acol[start_idx:end_idx]
                )
                
                # hepmc_events.append(event)
                
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
        
        hepmc_events = self._fill_batch_events_no_info(awkward_batch, start_event_num)
        
        # Add event info from the 'info' field
        with profile_block('Add event info'):
            for i, evt in enumerate(hepmc_events):
                info_data = awkward_batch[i]['info']
                self._add_event_info(evt, info_data)
        
        return hepmc_events

    def _pythia_to_hepmc_status(self, all_status):
        """
        Vectorized status code conversion for entire batch.
        Much faster than converting one by one.
        """
        # Create lookup table for common status codes
        status_map = {
            -11: 4, -12: 4, 11: 4, 12: 4,  # beam particles
        }
        
        # Start with default mappings
        hepmc_status = np.ones_like(all_status)  # default to 1 (final state)
        hepmc_status[all_status < 0] = 2  # intermediate state for negative status
        
        # Apply specific mappings vectorially
        for pythia_code, hepmc_code in status_map.items():
            hepmc_status[all_status == pythia_code] = hepmc_code
        
        # Handle ranges vectorially
        hard_process_mask = (all_status >= 21) & (all_status <= 29)
        hepmc_status[hard_process_mask] = 3  # outgoing from hard process
        
        # Final state particles from various processes
        final_state_ranges = [
            (31, 39), (41, 49), (51, 59), (61, 69), (71, 79), (81, 89)
        ]
        for start, end in final_state_ranges:
            range_mask = (all_status >= start) & (all_status <= end)
            hepmc_status[range_mask] = 1
        
        return hepmc_status

    def _build_event_vertices(self, event, event_particles, mother1, mother2, 
                                       x_prod, y_prod, z_prod, t_prod, vProd_is_none, offset):
        """
        Build vertices for a single event using vectorized data.
        """
        from pyHepMC3 import HepMC3 as hm
        
        n_particles = len(event_particles)
        vertex_cache = hm.GenEvent().vertices()
        beam_particles = hm.GenEvent().particles()
        
        for i in range(n_particles):
            # Reconstruct mother list from mother1/mother2 (adjust for global indexing)
            mothers_list = []
            if mother1[i] > 0:  # mother1 > 0 indicates valid mother
                local_mother1 = mother1[i] - 1 - offset  # Convert to local 0-based indexing
                if 0 <= local_mother1 < n_particles:
                    mothers_list.append(local_mother1)
                    
            if mother2[i] > 0 and mother2[i] != mother1[i]:
                local_mother2 = mother2[i] - 1 - offset  # Convert to local 0-based indexing
                if 0 <= local_mother2 < n_particles:
                    mothers_list.append(local_mother2)
            elif mother2[i] > mother1[i]:
                # Handle range case
                for m in range(mother1[i], mother2[i] + 1):
                    local_m = m - 1 - offset
                    if 0 <= local_m < n_particles and local_m not in mothers_list:
                        mothers_list.append(local_m)
            
            if len(mothers_list) > 0:
                # Find or create production vertex
                prod_vtx = None
                for mother_idx in mothers_list:
                    if event_particles[mother_idx].end_vertex():
                        prod_vtx = event_particles[mother_idx].end_vertex()
                        break
                
                if prod_vtx is None:
                    prod_vtx = hm.GenVertex()
                    vertex_cache.append(prod_vtx)
                    for mother_idx in mothers_list:
                        prod_vtx.add_particle_in(event_particles[mother_idx])
                
                # Set vertex position if available
                if prod_vtx.position().is_zero():
                    if not vProd_is_none[i]:
                        if not (x_prod[i] == 0 and y_prod[i] == 0 and z_prod[i] == 0 and t_prod[i] == 0):
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
            if col[i] == 0 and acol[i] == 0:
                # No color
                pass
            elif col[i] == 0 and acol[i] > 0:
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
        
        # You'll need to check what fields are available in info_data
        # This is a template based on typical Pythia8 info content
        
        try:
            # PDF information - adjust field names as needed
            if self.m_store_pdf and 'pdf' in info_data.fields:
                pdf_info = info_data['pdf']
                id1pdf = int(pdf_info['id1'])
                id2pdf = int(pdf_info['id2'])
                
                if self.m_convert_gluon_to_0:
                    if id1pdf == 21:
                        id1pdf = 0
                    if id2pdf == 21:
                        id2pdf = 0
                
                pdfinfo = hm.GenPdfInfo()
                pdfinfo.set(id1pdf, id2pdf, 
                           float(pdf_info['x1']), float(pdf_info['x2']), 
                           float(pdf_info['QFac']), 
                           float(pdf_info['pdf1']), float(pdf_info['pdf2']))
                evt.set_pdf_info(pdfinfo)
            
            # Process information - adjust field names as needed  
            if self.m_store_proc:
                if 'nMPI' in info_data.fields:
                    evt.add_attribute("mpi", hm.IntAttribute(int(info_data['nMPI'])))
                if 'code' in info_data.fields:
                    evt.add_attribute("signal_process_id", hm.IntAttribute(int(info_data['code'])))
                if 'QRen' in info_data.fields:
                    evt.add_attribute("event_scale", hm.DoubleAttribute(float(info_data['QRen'])))
                if 'alphaS' in info_data.fields:
                    evt.add_attribute("alphaQCD", hm.DoubleAttribute(float(info_data['alphaS'])))
                if 'alphaEM' in info_data.fields:
                    evt.add_attribute("alphaQED", hm.DoubleAttribute(float(info_data['alphaEM'])))
            
            # Cross-section - adjust field names as needed
            if self.m_store_xsec and 'sigmaGen' in info_data.fields:
                xsec = hm.GenCrossSection()
                sigma_gen = float(info_data['sigmaGen'])
                sigma_err = float(info_data['sigmaErr']) if 'sigmaErr' in info_data.fields else 0.0
                xsec.set_cross_section(sigma_gen * 1e9, sigma_err * 1e9)  # Convert to pb
                evt.set_cross_section(xsec)
            
            # Weights - this might need adjustment based on how weights are stored
            if self.m_store_weights and 'weights' in info_data.fields:
                evt.weights().clear()
                weights = info_data['weights']
                if hasattr(weights, '__iter__'):
                    for weight in weights:
                        evt.weights().append(float(weight))
                else:
                    evt.weights().append(float(weights))
                    
        except Exception as e:
            print(f"Warning: Could not add event info: {e}")
            # Continue without crashing