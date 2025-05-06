import os
import numpy as np
import ROOT as rt
import h5py as h5
import subprocess as sub
from util.pythia.utils import PythiaWrapper
from util.qol_utils.pdg import pdg_names, pdg_plotcodes, FillPdgHist
from util.calcs import Calculator
from util.particle_selection.algos import IsNeutrino
from util.hepmc import CreateHepMCEvent, CreateFullHepMCEvent, HepMCOutput
import util.qol_utils.qol_util as qu

class Generator:
    def __init__(self, pt_min, pt_max, configurator, pythia_rng=None, pythia_config_file=None, verbose=False):
        self.configurator = configurator
        self.pt_min = pt_min
        self.pt_max = pt_max

        # Create our Pythia wrapper.
        self.pythia = PythiaWrapper(verbose=verbose)
        self.pythia_rng = pythia_rng
        self.ConfigPythia(config_file=pythia_config_file)

        # Particle selectors and event filters.
        self.event_selection = None
        self.truth_selection = None
        self.final_state_selection = None
        self.event_filter = None
        self.jet_config = None
        self.n_truth = None

        # Things for the progress bar.
        self.prefix = 'Generating events for pT bin [{},{}]:'.format(self.pt_min,self.pt_max)
        self.suffix = 'Complete'
        self.bl = 50

        self.SetOutputDirectory()

        self.SetHistFilename('hists.root')
        self.SetFilename('events.hepmc')
        self.stats_filename = 'stats.h5'
        self.SetIndexOverlapFilename() # will give a default name
        self.SetEventFilterFlagFilename() # will give a default name
        self.filename_fullpath = None
        self.filename_truth_fullpath = None
        self.filename_full_fullpath = None # TODO: fix naming scheme

        self.diagnostic_plots = True
        self.InitializeHistograms()

        # Containers for event-level information.
        self.weights = None
        self.process_codes = None
        self.xsecs = None

        self.progress_bar = True

        self.loop_number = 0 # used for keeping track of successful generation loops
        self.nevents = None # number of events requested, will be set in generation function

        self.hepev_buffer_fs = [] # list containing HepMC events, to be written to a file
        self.hepev_buffer_truth = []
        self.hepev_buffer_full = []
        self.SetBufferSize(100) # TODO: Should this be configurable? Could be too much detail.
        self.buffername = None
        self.buffername_truth = None
        self.nevents_success = 0 # number of events successfully generated
        self.nevents_failed = 0 # number of events that failed (failure to pass basic cuts)

        self.header_status = False
        self.footer_status = False

        self.full_hepmc = self.configurator.GetMakeFullHepMCFile() # whether or not to create full HepMC3 files, that contain all the event's particles and all indices

        # Variables related to handling of pileup.
        # TODO: Might want to eventually reorganize this, or move it into a separate class?
        self.pileup_handler = self.configurator.GetPileupHandling()

        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetVerbose(self,val):
        self.pythia.SetVerbose(val)

    def SetEventSelection(self,selection):
        self.event_selection = selection
        try: self.event_selection.Initialize(self.configurator)
        except: pass

    def SetTruthSelection(self,selection):
        self.truth_selection = selection

    def SetFinalStateSelection(self,selection):
        self.final_state_selection = selection

    def SetEventFilter(self,filter):
        self.event_filter = filter
        if(self.event_filter is not None):
            self.event_filter.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import

    def SetEventFilterFlag(self,filter):
        self.event_filter_flag = filter
        if(self.event_filter_flag is not None):
            self.event_filter_flag.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import
        else:
            self.event_filter_flag_filename = None

    def SetJetConfig(self,config):
        self.jet_config = config

    def SetNTruth(self,n):
        self.n_truth = n

    def SetPythiaConfigFile(self,file=None):
        self.pythia_config_file = file
        if(file is None):
            self.pythia_config_file = self.configurator.GetPythiaConfigFile() # picked up from dictionary in config/config.py

    def SetDefaultFilenames(self,outdir=None):
        self.SetOutputDirectory(dir)
        self.SetFilename('events.hepmc')
        self.SetHistFilename('hists.root')
        self.stats_filename = 'stats.h5'
        self.SetIndexOverlapFilename()
        self.SetEventFilterFlagFilename()

    def SetProgressBar(self,flag):
        self.progress_bar = flag

    def SetOutputDirectory(self,dir=None):
        if(dir is None): dir = os.getcwd()
        self.outdir = dir

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetFilename(self,name, rename_extra_files=True):
        if('.hepmc' not in name):
            name = '{}.hepmc'.format(name)
        self.filename = name
        self.SetTruthFilename()
        if(rename_extra_files):
            self.SetIndexOverlapFilename(None)
            self.SetEventFilterFlagFilename(None)
        return

    def SetTruthFilename(self,name=None):
        if(name is None):
            self.truth_filename = self.filename.replace('.hepmc','_truth.hepmc')
            return

        if('.hepmc' not in name):
            name = '{}.hepmc'.format(name)
        self.truth_filename = name
        return

    def SetStatsFilename(self,name):
        if('.h5' not in name):
            name = '{}.h5'.format(name)
        self.stats_filename = name
        return

    def GetFilename(self):
        return self.filename

    def GetTruthFilename(self):
        return self.truth_filename

    def GetStatsFilename(self):
        return self.stats_filename

    def GetHistFilename(self):
        return self.hist_filename

    def SetIndexOverlapFilename(self,name=None):
        if(name is None): self.fs_truth_overlap_filename = self.filename.replace('.hepmc','_final-state_truth_overlap_indices.h5')
        else: self.fs_truth_overlap_filename = name
        return

    def SetEventFilterFlagFilename(self,name=None):
        if(name is None): self.event_filter_flag_filename = self.filename.replace('.hepmc','_event_filter_flag.h5')
        else: self.event_filter_flag_filename = name
        return

    def GetIndexOverlapFilename(self):
        return self.fs_truth_overlap_filename

    def GetEventFilterFlagFilename(self):
        return self.event_filter_flag_filename

    def ConfigPythia(self, config_file=None):
        """
        Prepare and apply Pythia configuration. This turns our settings (from our config file)
        into a list of strings ready to be input to Pythia8.
        """
        self.pythia_config = self.configurator.GetPythiaConfig(self.pt_min,self.pt_max)
        self.SetPythiaConfigFile(file=config_file)

        # Optionally set the Pythia RNG seed to something other than what's in the config.
        # TODO: Can we make this more tidy?
        if(self.pythia_rng is not None):
            self.pythia_config['Random:seed'] = self.pythia_rng

        # Now apply these configurations to our Generator's instance of PythiaWrapper.
        self.pythia.ClearConfigDict()
        self.pythia.AddToConfigDict(self.pythia_config)
        # self.pythia.ReadConfigDict()
        self.pythia.ReadStringsFromFile(self.pythia_config_file)
        self.pythia.InitializePythia()

    def InitializeHistograms(self):
        self.hists = []
        self.hist_fs_pdgid = rt.TH1D(qu.RN(),'Final-state particles;Particle;Count',47,0.,47)
        self.hists.append(self.hist_fs_pdgid)

        # Initialize histograms for sum of energy, E_T and p_T for final-state particles.
        self.kinematic_hists = {}
        hist_binning = (1000,0.,1000.)

        for key,key_insert in zip(('all','invisible','visible'),('',' (invisibles)',' (visibles)')):
            d = {}
            d['e'] = rt.TH1D(qu.RN(),'#sum Energy {};E [GeV];Count'.format(key_insert),*hist_binning)
            d['et'] = rt.TH1D(qu.RN(),'#sum ET {};ET [GeV];Count'.format(key_insert).replace('ET','E_{T}'),*hist_binning)
            d['pt'] = rt.TH1D(qu.RN(),'#sum pT {};pT [GeV];Count'.format(key_insert).replace('pT','p_{T}'),*hist_binning)
            self.kinematic_hists[key] = d
            for h in d.values():
                self.hists.append(h)

    def OutputHistograms(self):
        hist_filename = '{}/{}'.format(self.outdir,self.hist_filename)
        f = rt.TFile(hist_filename,'UPDATE')
        for i,h in enumerate(self.hists):
            canvas_name = 'c_{}'.format(i)
            hist_name = 'h_{}'.format(i)
            c = rt.TCanvas(canvas_name,'',800,600)

            if(i == 0): # pdg_hist
                pdg_names_inv = {}
                for key,val in pdg_plotcodes.items():
                    code = val
                    name = pdg_names[key]
                    pdg_names_inv[code] = name

                n = h.GetNbinsX()
                xaxis = h.GetXaxis()
                for j in range(n):
                    if(j == 0): name = ''
                    else:
                        try: name = pdg_names_inv[j]
                        except: name = ''
                    xaxis.SetBinLabel(j+1,name)
                xaxis.SetLabelSize(2.0e-2)
                rt.gPad.SetGrid(1,0)
            else: rt.gPad.SetGrid()

            h.Draw('HIST')
            rt.gPad.SetLogy()
            c.Draw()
            f.cd() # unnecessary
            # c.Write(canvas_name)
            h.Write(hist_name)
        f.Close()
        return

    def SetBufferSize(self,size=100):
        self.buffer_size = size

    def GetCurrentBufferSize(self):
        return len(self.hepev_buffer_fs)

    def ClearEventBuffers(self):
        self.hepev_buffer_fs.clear()
        self.hepev_buffer_truth.clear()
        self.hepev_buffer_full.clear()

    def FillEventBuffers(self,hepev_fs,hepev_truth=None,hepev_full=None):
        self.hepev_buffer_fs.append(hepev_fs)
        if(hepev_truth is not None): self.hepev_buffer_truth.append(hepev_truth)
        if(hepev_full is not None): self.hepev_buffer_full.append(hepev_full)

    def WriteEventBuffersToFile(self,header=False,footer=False):
        if(header): self.header_status = True
        if(footer): self.footer_status = True
        HepMCOutput(self.hepev_buffer_fs,self.buffername,self.filename_fullpath,header,footer)
        HepMCOutput(self.hepev_buffer_truth,self.buffername_truth,self.filename_truth_fullpath,header,footer)
        if(self.full_hepmc):
            HepMCOutput(self.hepev_buffer_full,self.buffername_full,self.filename_full_fullpath,header,footer)
        self.ClearEventBuffers()

    def GenerationLoop(self, nevents,i_real = 1, nevents_disp=None):
        """
        This is the function where Pythia8 matrix element generation + showering/hadronization happens.
        The results are filtered for the selected "truth" and "final state" particles, which are placed
        into HepMC3 events. These events are periodically written to a buffer file (which is then merged
        into the "master" HepMC3 file).
        """
        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that pyhepmc's WriterAscii works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality.

        self.filename_fullpath = '{}/{}'.format(self.outdir,self.filename)
        self.filename_truth_fullpath = self.filename_fullpath.replace('.hepmc','_truth.hepmc')
        self.filename_full_fullpath = self.filename_fullpath.replace('.hepmc','_full.hepmc')
        self.buffername = self.filename_fullpath.replace('.hepmc','_buffer.hepmc')
        self.buffername_truth = self.buffername.replace('.hepmc','_truth.hepmc')
        self.buffername_full = self.buffername.replace('.hepmc','_full.hepmc')

        for i in range(nevents):

            self.pythia.Generate() # generate an event!

            # Get the truth-level particles, using our truth selector.
            if(self.truth_selection is not None):
                truth_indices = np.atleast_1d(np.array(self.truth_selection(self.pythia),dtype=int))
                truth_selection_status = self.truth_selection.GetSelectionStatus()
            else:
                truth_indices = np.zeros(0)
                truth_selection_status = True # no truth particles, but it is okay

            # Some checks -- if we expect a fixed # of truth-level particles from our selector,
            # and it didn't give this, it means that something went wrong. In that case, we should
            # skip this event.
            if(not truth_selection_status):
                n_fail += 1
                continue

            # Here, we optionally apply pileup (or more generally, merging in some other event(s)).
            if(self.pileup_handler is not None):
                self.pileup_handler(self.pythia)

            # Get the final-state particles that we're saving (however we've defined our final state!).
            final_state_indices = np.atleast_1d(np.array(self.final_state_selection(self.pythia),dtype=int))
            final_state_selection_status = self.final_state_selection.GetSelectionStatus()
            if(not final_state_selection_status):
                n_fail += 1
                continue

            # ===================================+======
            # Now we apply an (optional) "event selection". This selects only a subset of the final-state particles
            # to save to our output HepMC file. In practice this can be useful for reducing these files' sizes.
            # For example, one may choose to only save final-state particles within some distance of one of
            # the selected truth-level particles.
            # ==========================================
            if(self.event_selection is not None):
                try: # TODO Running into some very rare/stochastic (?) error with certain event selections (using VectorCalcs.DeltaR2Vectorized()).
                    final_state_indices = self.event_selection(self.pythia, final_state_indices, truth_indices)
                except: # easiest thing to do for this kind of exception is to just throw out the event and try again
                    n_fail += 1
                    continue

            # ==========================================
            # Now we apply an (optional) "event filter". This applies some condition to the set
            # of truth and final-state particles selected above, and requires that the event pass
            # this condition. If not we will count the event as failed, and generate another.
            # ==========================================
            if(self.event_filter is not None):
                passed_filter = self.event_filter(self.pythia, final_state_indices)
                if(not passed_filter):
                    n_fail += 1
                    continue

            # ==========================================
            # Now we apply an (optional) "event filter" flag. This applies some condition to the set
            # of truth and final-state particles selected above, and checks if the event passes
            # this condition. The result is recorded in an HDF5 file that will be merged into the
            # final dataset.
            # ==========================================
            if(self.event_filter_flag is not None):
                event_filter_flag_filename_full = '{}/{}'.format(self.outdir,self.event_filter_flag_filename)
                f = h5.File(event_filter_flag_filename_full,'a')
                key = self.event_filter_flag.GetName()
                filter_flag = self.event_filter_flag(self.pythia, final_state_indices)
                if(i_real == 1):
                    filter_flag = np.expand_dims(filter_flag,axis=0)
                    f.create_dataset(key,data=filter_flag,compression='gzip',chunks=True,maxshape=(None,))
                else:
                    f[key].resize(f[key].shape[0] + 1,axis=0)
                    f[key][-1] = filter_flag
                f.close()

            # ==========================================
            # In some cases we may be interested in particles which are counted as both truth particles and final-state particles.
            # For example, we may be using our "truth particles" list to store the stable daughter particles of a particular
            # particle (like the stable daughters of a W boson decay).
            # These particles may be interesting to keep track of if we're using Delphes -- the Delphes output contains information
            # on the indices of particles that hit each detector element, and we may need to map between the detector
            # hits and any (stable) particles within among our truth particle array. To do this, we will ultimately
            # need to determine these truth particles' indices w.r.t. the final_state_hepev that we create further below.
            # ==========================================
            self.RecordFinalStateTruthOverlap(final_state_indices,truth_indices,i_real,'indices',self.configurator.GetNPars()['jet_n_par'])

            # ==========================================
            # Now lets create HepMC events -- one holding just the truth-level particles
            # and one holding just the final-state particles (however we've defined our
            # "truth selection" and "final-state selection" in our configuration).
            # By keeping these separate, we make it easy to (later) determine which
            # particles are passed to jet clustering, at the cost of each event being
            # split across two files.
            #
            # Keep in mind that these are also not totally complete HepMC files as we
            # have optionally applied an "event selection" above that may have thrown out
            # some particles (that are hopefully not relevant to whatever jets we're studying),
            # and partly as a consequence of this event reduction (as well as the fact that we
            # split things across two files) we're also not saving vertex information to
            # these HepMC files.
            #
            # TODO: Might want to include vertex info! This is being done with CreateFullHepMCEvent,
            #       but there might be some issues with it (events produced with it cannot be
            #       be visualized by pyhepmc, suggests some vertices might be missing from output.)
            # ==========================================

            final_state_hepev = CreateHepMCEvent(self.pythia,final_state_indices,i_real)

            truth_hepev = None
            if(len(truth_indices) > 0): # TODO: If this condition isn't met for whatever reason, will final-state and truth files' events not line up?
                truth_hepev = CreateHepMCEvent(self.pythia,truth_indices, i_real)

            full_hepev = None
            if(self.full_hepmc):
                full_hepev = CreateFullHepMCEvent(self.pythia,i_real) # TODO: Adding full event recording, with vertices. Currently unused in pipeline.

            # Fill the memory buffer with this event.
            self.FillEventBuffers(final_state_hepev, truth_hepev,full_hepev)

            # If buffer is full (i.e. has reached max size), write it to file & flush.
            if(self.GetCurrentBufferSize() == self.buffer_size):
                header = (self.loop_number == 0) and (not self.header_status)
                footer = False
                self.WriteEventBuffersToFile(header,footer)

            # Diagnostic plots.
            # TODO: Probably a lot of optimization to do for these, some new options thanks to PythiaWrapper.
            # if(self.diagnostic_plots):
                # Record the PDG codes of all selected final-state particles (i.e. all those that made it into the HepMC output).
                # HistogramFinalStateCodes(arr,self.hist_fs_pdgid)
                # HistogramFinalStateKinematics(arr,self.kinematic_hists)

            if(self.progress_bar): qu.printProgressBarColor(i_real,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)

            # Record this event's weight.
            self.weights[i_real-1] = self.pythia.GetEventWeight() # convert from 1-indexing to 0-indxing

            # Record this event's process code.
            self.process_codes[i_real-1] = self.pythia.GetProcessCode() # convert from 1-indexing to 0-indxing

            i_real += 1 # If success, increase i_real -- this is a counter for the number of successful events

        # Buffer gets written to the file in the loop above whenever it's full, but after exiting the loop
        # we should flush it once more -- even if empty, since this is still where we may write the footer.
        header = (self.loop_number == 0) and (not self.header_status) # for writing HepMC header
        footer = n_fail == 0
        self.WriteEventBuffersToFile(header,footer)

        # Delete the buffer files.
        for fname in [self.buffername,self.buffername_truth]:
            comm = ['rm', fname]
            try: sub.check_call(comm,stderr=sub.DEVNULL)
            except: pass

        return i_real-1, n_fail # note that i_real is using 1-indexing, which is what HepMC events use

    # Generate a bunch of events in the given pT range,
    # and save them to a HepMC file.
    # We do perform event selection: Only certain particles are saved to the file to begin with.
    def Generate(self,nevents):
        self.nevents = nevents # total number of events we request

        # File for keeping track of any particle indices that correspond with particles saved in *both*
        # our "final-state" and "truth" selections. Indices are stored in two ways:
        #   1) with respect to how they appear in the final-state HepMC3 files, and
        #   2) with respect to how they appear in the truth selection HepMC3 files.
        #
        # This may be useful for studying Delphes output, e.g. keeping track of which calorimeter towers were
        # hit by stable daughters of a W-boson (which were selected by "final-state" and "truth" selectors).
        fs_truth_overlap_filename_full = '{}/{}'.format(self.outdir,self.fs_truth_overlap_filename)
        try: sub.check_call(['rm',fs_truth_overlap_filename_full],stderr=sub.DEVNULL)
        except: pass

        # # Get the Fastjet banner out of the way # TODO: This should be done elsewhere.
        # tmp = InitFastJet()
        # del tmp

        if(self.progress_bar): qu.printProgressBarColor(0,nevents, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # Loop in such a way as to guarantee that we get as many events as requested.
        # This logic is required as events could technically fail selections, e.g. not have the
        # requested truth particles (depends on requested truth particles & processes).
        self.nevents_success = 0
        n_fail = nevents
        self.loop_number = 0
        self.weights       = np.zeros(nevents)
        self.process_codes = np.zeros(nevents, dtype=int)
        while(n_fail > 0):
            self.nevents_success, n_fail = self.GenerationLoop(
                nevents-self.nevents_success,
                i_real=self.nevents_success+1,
                nevents_disp = nevents
            )
            self.loop_number += 1

        # Fill out an array with each event's cross section (and the uncertainty on that cross-section).
        self.xsecs = np.zeros((nevents,2))
        xsec_dictionary = self.GetSigmaDictionary()
        for i,pcode in enumerate(self.process_codes):
            self.xsecs[i,:] = xsec_dictionary[pcode]

        # Create/append a stats file, and put in information on event weights & cross-sections.
        f = h5.File('{}/{}'.format(self.outdir,self.stats_filename),'a')
        compression = 'gzip'
        copts = 9
        f.create_dataset('mc_weight',data=self.weights,compression=compression,compression_opts=copts)
        f.create_dataset('process_code',data=self.process_codes,compression=compression,compression_opts=copts)
        f.create_dataset('cross_section',data=self.xsecs[:,0],compression=compression,compression_opts=copts)
        f.create_dataset('cross_section_uncertainty',data=self.xsecs[:,1],compression=compression,compression_opts=copts)
        f.close()

        # if(self.diagnostic_plots): self.OutputHistograms()
        return

    def RecordFinalStateTruthOverlap(self,final_state_indices,truth_indices,i_real,key='indices',indices_max=200):

        # Get indices (w.r.t. full Pythia event) of particles that are listed in both the final_state and truth.
        final_state_truth_overlap = np.intersect1d(final_state_indices,truth_indices)

        # We want to make an (n,2) array, where for n listings, we give
        # - the final_state_index in position 0,
        # - and the truth_index in position 1.
        # This is complicated by the possibility of a particle appearing multiple times in either list of indices.
        # (e.g. the truth selection selects decay products of W, but also all stable W daughters -- and the W decays to lepton + neutrino).
        # TODO: We currently deal with this complication via the "[0][0]" indexing in fs_truth_overlap_wrt_fs and fs_truth_overlap_wrt_truth.
        #       This means that we just count the first instance of the particle's listing. Maybe causes potential confusion in the case where
        #       it actually appears multiple times.

        # When final_state_indices is written to final_state_hepev, the particles will be re-indexed in the order they're given.
        # Note that HepMC3 uses 1-indexing, we will save these with 1-indexing too.
        # TODO: Seems to cause crash if a final-state particle shows up multiple times in the truth listing.
        fs_truth_overlap_wrt_fs    = np.array([np.where(final_state_indices==x)[0][0] + 1 for x in final_state_truth_overlap]).flatten()
        fs_truth_overlap_wrt_truth = np.array([np.where(truth_indices==x)[0][0] + 1 for x in final_state_truth_overlap]).flatten()

        indices = np.zeros((indices_max,2),dtype=int) # order will be (index w.r.t. final-state HepMC file, index w.r.t. truth-selection HepMC file)
        l = np.minimum(indices_max,len(fs_truth_overlap_wrt_fs))
        indices[:l,0] = fs_truth_overlap_wrt_fs[:l]
        indices[:l,1] = fs_truth_overlap_wrt_truth[:l]

        fs_truth_overlap_filename_full = '{}/{}'.format(self.outdir,self.fs_truth_overlap_filename)
        f = h5.File(fs_truth_overlap_filename_full,'a')
        if(i_real == 1):
            indices = np.expand_dims(indices,axis=0)
            f.create_dataset(key,data=indices,compression='gzip',chunks=True,maxshape=(None,indices_max,2))
        else:
            f[key].resize(f[key].shape[0] + 1,axis=0)
            f[key][-1] = indices
        f.close()


    # Get the MC event weights (will typically just be 1 for each event).
    # This reads from a file, which is only written to when an event is written to disk.
    # Thus these weights will line up with the events we've saved, i.e. we don't have to worry
    # about events that Pythia8 successfully generated but which we threw out because they failed
    # one of our selectors.
    def GetEventWeights(self):
        return self.weights

    def GetProcessCodes(self):
        return self.process_codes

    # This returns a list of all unique process codes encountered,
    # not a list of per-event process codes.
    def GetUniqueProcessCodes(self):
        return self.pythia.GetProcessCodes()

    # Get a dictionary containing cross sections (and their uncertainties)
    # for every process that was run, organized by process code.
    # Note that turning on a single Pythia process flag can in principle
    # turn on multiple processes, e.g. HardQCD will provide for many different
    # processes and the codes will distinguish between them.
    # Cross-sections are given in mb.
    # TODO: The estimates will possibly be inaccurate, as we throw out certain
    # events that Pythia8 has generated when they fail our particle selectors, but
    # those thrown out events are still included in computing the cross-section estimates.
    def GetSigmaDictionary(self):
        return self.pythia.GetSigmaDictionary()

    # def HistogramFinalStateCodes(self,arr):
    #     codes = []
    #     for entry in arr:
    #         codes.append(entry[-2])
    #     FillPdgHist(self.hist_fs_pdgid,codes)
    #     return

    # Some utility functions, which were previously in a separate sub-library.
    # -- plotting functionality below --
    def HistogramFinalStateCodes(self,arr,hist):
        codes = []
        for entry in arr:
            codes.append(entry[-2])
        FillPdgHist(hist,codes)
        return

    def HistogramFinalStateKinematics(self,arr,kinematic_hist_dict):
        # We will make histograms of sum E, E_T and p_T. For each, we consider a sum over all particles,
        # as well as separate sums over invisibles (i.e. neutrinos) and visibles.
        full_sum = np.zeros(4)
        inv_sum = np.zeros(4)
        vis_sum = np.zeros(4)

        for par in arr:
            par_array = np.array([par[0],par[1],par[2],par[3]] )
            full_sum += par_array
            invisible = IsNeutrino(p=par) # TODO: Not sure if this will still work based on other updates to the code?
            if(invisible): inv_sum += par_array
            else: vis_sum += par_array

        for vec,key in zip((full_sum,inv_sum,vis_sum),('all','invisible','visible')):
            e = vec[0]
            px = vec[1]
            py = vec[2]
            pz = vec[3]
            vec_cyl = self.calculator.PxPyPzEToPtEtaPhiM([px],[py],[pz],[e])[0]
            pt = vec_cyl[0]
            et = e / np.cosh(vec_cyl[1])
            kinematic_hist_dict[key]['e'].Fill(e)
            kinematic_hist_dict[key]['et'].Fill(et)
            kinematic_hist_dict[key]['pt'].Fill(pt)