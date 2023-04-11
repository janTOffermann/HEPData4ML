import os
import numpy as np
import ROOT as rt
import h5py as h5
import pyhepmc as hep
import subprocess as sub

from util.pythia.utils import PythiaWrapper
from util.qol_utils.pdg import pdg_names, pdg_plotcodes, FillPdgHist
from util.calcs import PxPyPzEToPtEtaPhiM
from util.particle_selection.algos import IsNeutrino
from util.hepmc import CreateHepMCEvent, HepMCOutput
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
        self.filename_full = None
        self.filename_truth_full = None

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
        self.SetBufferSize(100)
        self.buffername = None
        self.buffername_truth = None
        self.nevents_success = 0 # number of events successfully generated
        self.nevents_failed = 0 # number of events that failed (failure to pass basic cuts)

    def SetEventSelection(self,selection):
        self.event_selection = selection

    def SetTruthSelection(self,selection):
        self.truth_selection = selection

    def SetFinalStateSelection(self,selection):
        self.final_state_selection = selection

    def SetEventFilter(self,filter):
        self.event_filter = filter
        self.event_filter.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import

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

    def SetProgressBar(self,flag):
        self.progress_bar = flag

    def SetOutputDirectory(self,dir=None):
        if(dir is None): dir = os.getcwd()
        self.outdir = dir

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetFilename(self,name):
        if('.hepmc' not in name):
            name = '{}.hepmc'.format(name)
        self.filename = name
        self.SetTruthFilename()
        self.SetIndexOverlapFilename(None)
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

    def SetIndexOverlapFilename(self,name=None):
        if(name is None): self.fs_truth_overlap_filename = self.filename.replace('.hepmc','_final-state_truth_overlap_indices.h5')
        else: self.fs_truth_overlap_filename = name
        return

    def GetIndexOverlapFilename(self):
        return self.fs_truth_overlap_filename

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

    def ClearEventBuffers(self):
        self.hepev_buffer_fs.clear()
        self.hepev_buffer_truth.clear()

    def FillEventBuffers(self,hepev_fs,hepev_truth=None):
        self.hepev_buffer_fs.append(hepev_fs)
        if(hepev_truth is not None): self.hepev_buffer_truth.append(hepev_truth)

    def WriteEventBuffersToFile(self,header=False,footer=False):
        HepMCOutput(self.hepev_buffer_fs,self.buffername,self.filename_full,header,footer)
        HepMCOutput(self.hepev_buffer_truth,self.buffername_truth,self.filename_truth_full,header,footer)
        self.ClearEventBuffers()

    def GenerationLoop(self, nevents,i_real = 1, nevents_disp=None):
        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that pyhepmc's WriterAscii works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality.

        self.filename_full = '{}/{}'.format(self.outdir,self.filename)
        self.filename_truth_full = self.filename_full.replace('.hepmc','_truth.hepmc')
        self.buffername = self.filename_full.replace('.hepmc','_buffer.hepmc')
        self.buffername_truth = self.buffername.replace('.hepmc','_truth.hepmc')

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
                final_state_indices = self.event_selection(self.pythia, final_state_indices, truth_indices)

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
            # In some cases we may be interested in particles which are counted as both truth particles and final-state particles.
            # For example, we may be using our "truth particles" list to store the stable daughter particles of a particular
            # particle (like the stable daughters of a W boson decay).
            # These particles may be interesting to keep track of if we're using Delphes -- the Delphes output contains information
            # on the indices of particles that hit each detector element, and in practice we may want to be able to map between
            # the detector hits and any (stable) particles we within among our truth particle array. To do this, we will ultimately
            # need to determine these truth particles' indices w.r.t. the final_state_hepev that we create further below.
            # ==========================================
            final_state_truth_overlap = np.intersect1d(final_state_indices,truth_indices) # still gives indices w.r.t. the full Pythia event, not what we want to save.

            # When final_state_indices is written to final_state_hepev, the particles will be re-indexed in the order they're given.
            # Note that HepMC3 uses 1-indexing, we will save these with 1-indexing too.
            fs_truth_overlap_wrt_fs    = np.array([np.where(final_state_indices==x)[0] + 1 for x in final_state_truth_overlap]).flatten()
            fs_truth_overlap_wrt_truth = np.array([np.where(truth_indices==x)[0] + 1 for x in final_state_truth_overlap]).flatten()

            indices_max = 200 # TODO: Make this configurable!
            indices = np.zeros((indices_max,2),dtype=int) # order will be (index w.r.t. final-state HepMC file, index w.r.t. truth-selection HepMC file)
            l = np.minimum(indices_max,len(fs_truth_overlap_wrt_fs))
            indices[:l,0] = fs_truth_overlap_wrt_fs[:l]
            indices[:l,1] = fs_truth_overlap_wrt_truth[:l]

            fs_truth_overlap_filename_full = '{}/{}'.format(self.outdir,self.fs_truth_overlap_filename)
            f = h5.File(fs_truth_overlap_filename_full,'a')
            key = 'indices'
            if(i_real == 1):
                indices = np.expand_dims(indices,axis=0)
                f.create_dataset(key,data=indices,compression='gzip',chunks=True,maxshape=(None,indices_max,2))
            else:
                f[key].resize(f[key].shape[0] + 1,axis=0)
                f[key][-1] = indices
            f.close()

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
            # TODO: We may want to totally rework this part of the code, so that the HepMC files
            #       contain the *full* event (incl. vertex information), and the selections
            #       (event selection, truth selection, final-state selection) are applied later
            #       when in the HepMC -> Delphes/HDF5 pipeline. The current method has the advantage
            #       of reducing the HepMC filesize but what we get aren't really "full" HepMC files,
            #       plus redesigning things to use more standard HepMC files (with vertex info) would
            #       make it easier to utilize HepMC files from outside sources (e.g. events generated
            #       externally with a different generator like MadGraph).
            # ==========================================
            final_state_hepev = CreateHepMCEvent(self.pythia,final_state_indices,i_real)
            # HepMCOutput(final_state_hepev,self.buffername,filename_full,self.loop_number,i,i_real,nevents,n_fail) # writes event to disk

            truth_hepev = None
            if(len(truth_indices) > 0): # TODO: If this condition isn't met for whatever reason, will final-state and truth files' events not line up?
                truth_hepev = CreateHepMCEvent(self.pythia,truth_indices, i_real)
                # HepMCOutput(truth_hepev,self.buffername_truth,filename_truth,self.loop_number,i,i_real,nevents,n_fail) # writes event to disk

            # Fill the memory buffer with this event.
            self.FillEventBuffers(final_state_hepev, truth_hepev)

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

        # Write the event buffer to the event buffer files. Then, copy info from the buffer files to the full HepMC files.
        header = self.loop_number == 0 # for writing HepMC header
        footer = i_real >= self.nevents
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

        # Create a stats file, and put in information on event weights & cross-sections.
        f = h5.File('{}/{}'.format(self.outdir,self.stats_filename),'w')
        compression = 'gzip'
        copts = 7
        f.create_dataset('mc_weight',data=self.weights,compression=compression,compression_opts=copts)
        f.create_dataset('process_code',data=self.process_codes,compression=compression,compression_opts=copts)
        f.create_dataset('cross_section',data=self.xsecs[:,0],compression=compression,compression_opts=copts)
        f.create_dataset('cross_section_uncertainty',data=self.xsecs[:,1],compression=compression,compression_opts=copts)
        f.close()

        # if(self.diagnostic_plots): self.OutputHistograms()
        return

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
            invisible = IsNeutrino(par,use_hepmc=False)
            if(invisible): inv_sum += par_array
            else: vis_sum += par_array

        for vec,key in zip((full_sum,inv_sum,vis_sum),('all','invisible','visible')):
            e = vec[0]
            px = vec[1]
            py = vec[2]
            pz = vec[3]
            vec_cyl = PxPyPzEToPtEtaPhiM([px],[py],[pz],[e])[0]
            pt = vec_cyl[0]
            et = e / np.cosh(vec_cyl[1])
            kinematic_hist_dict[key]['e'].Fill(e)
            kinematic_hist_dict[key]['et'].Fill(et)
            kinematic_hist_dict[key]['pt'].Fill(pt)