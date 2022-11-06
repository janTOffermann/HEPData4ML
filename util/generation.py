import sys, glob
import numpy as np
import ROOT as rt
import subprocess as sub

# import numpythia as npyth # Pythia, hepmc_write
from util.pythia.utils import PythiaWrapper

from util.config import GetEventSelection, GetTruthSelection, GetFinalStateSelection, GetPythiaConfig, GetPythiaConfigFile, GetJetConfig, GetNPars
from util.qol_utils.pdg import pdg_names, pdg_plotcodes
from util.gen_utils.utils import CreateHepMCEvent, HepMCOutput, HistogramFinalStateKinematics, HistogramFinalStateCodes
from util.conv_utils.utils import InitFastJet
import util.qol_utils.qol_util as qu

class Generator:
    def __init__(self, pt_min, pt_max):
        self.pt_min = pt_min
        self.pt_max = pt_max

        # Create our Pythia wrapper.
        self.pythia = PythiaWrapper()
        self.ConfigPythia()

        # Configure our particle selectors and event filters.
        self.event_selection = GetEventSelection()
        self.truth_selection = GetTruthSelection()
        self.final_state_selection = GetFinalStateSelection()
        self.jet_config = GetJetConfig()
        self.n_truth = GetNPars()['n_truth']

        # Things for the progress bar.
        self.prefix = 'Generating events for pT bin [{},{}]:'.format(self.pt_min,self.pt_max)
        self.suffix = 'Complete'
        self.bl = 50

        self.outdir = None # TODO: Use this with filenames.

        self.hist_filename = 'hists.root'
        self.SetFilename('events.hepmc')

        self.diagnostic_plots = True
        self.InitializeHistograms()

    def SetOutputDirectory(self,dir):
        self.outdir = dir

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

    def SetFilename(self,name, truth=True):
        if('.hepmc' not in name):
            name = '{}.hepmc'.format(name)
        self.filename = name
        if(truth):
            self.SetTruthFilename(name.replace('.hepmc','_truth.hepmc'))
        return

    def SetTruthFilename(self,name):
        if('.hepmc' not in name):
            name = '{}.hepmc'.format(name)
        self.truth_filename = name
        return

    def GetFilename(self):
        return self.filename

    def GetTruthFilename(self):
        return self.truth_filename

    def ConfigPythia(self):
        """
        Prepare and apply Pythia configuration. This turns our settings (from our config file)
        into a list of strings ready to be input to Pythia8.
        """
        self.pythia_config = GetPythiaConfig(self.pt_min,self.pt_max)
        self.pythia_config_file = GetPythiaConfigFile()

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

    def GenerationLoop(self, nevents,filename,i_real = 1, nevents_disp=None,loop_number=0):
        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that pyhepmc_ng's WriterAscii works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality.

        filename_truth = filename.replace('.hepmc','_truth.hepmc')
        buffername = filename.replace('.hepmc','_buffer.hepmc')
        buffername_truth = buffername.replace('.hepmc','_truth.hepmc')

        for i in range(nevents):
            success = True
            self.pythia.Generate() # generate an event!

            # Get the truth-level particles, using our truth selector.
            truth_indices = np.atleast_1d(np.array(self.truth_selection(self.pythia),dtype=int))
            truth_selection_status = self.truth_selection.GetSelectionStatus()

            # Some checks -- if we expect a fixed # of truth-level particles from our selector,
            # and it didn't give this, it means that something went wrong. In that case, we should
            # skip this event.
            if(not truth_selection_status):
                success = False
                n_fail += 1
                continue

            # Get the final-state particles that we're saving (however we've defined our final state!).
            final_state_indices = np.atleast_1d(np.array(self.final_state_selection(self.pythia),dtype=int))
            final_state_selection_status = self.final_state_selection.GetSelectionStatus()
            if(not final_state_selection_status):
                success = False
                n_fail += 1
                continue

            # ==========================================
            # Now we apply an (optional) "event filter". This selects only a subset of the final-state particles
            # to save to our output HepMC file. In practice this can be useful for reducing these files' sizes.
            # For example, one may choose to only save final-state particles within some distance of one of
            # the selected truth-level particles.
            # ==========================================

            # print(final_state_indices)
            if(self.event_selection is not None):
                final_state_indices = self.event_selection(self.pythia, final_state_indices, truth_indices)
            # print(final_state_indices)

            # ==========================================
            # Now lets create HepMC events -- one holding just the truth-level particles,
            # and one holding just the final-state particles. By keeping these separate,
            # we make it easy to (later) determine which particles are passed to jet clustering.
            # ==========================================
            final_state_hepev = CreateHepMCEvent(self.pythia,final_state_indices,i_real) # internally will fetch momenta
            truth_hepev       = CreateHepMCEvent(self.pythia,truth_indices,      i_real)

            # With the HepMC events created, we write them to disk.
            # TODO: This is probably a bottleneck since it involves I/O, is there a clever way
            #       for us to chunk this part and write a bunch of HepMC events to a memory buffer first?

            # ----- File I/O -----
            # Write this event to the HepMC buffer file, then copy the buffer contents to the full HepMC file.
            HepMCOutput(final_state_hepev,buffername,filename,loop_number,i,i_real,nevents,n_fail)
            HepMCOutput(truth_hepev,buffername_truth,filename_truth,loop_number,i,i_real,nevents,n_fail)
            # ----- End File I/O -----

            # Diagnostic plots.
            # TODO: Probably a lot of optimization to do for these, some new options thanks to PythiaWrapper.
            # if(self.diagnostic_plots):
                # Record the PDG codes of all selected final-state particles (i.e. all those that made it into the HepMC output).
                # HistogramFinalStateCodes(arr,self.hist_fs_pdgid)
                # HistogramFinalStateKinematics(arr,self.kinematic_hists)

            qu.printProgressBarColor(i_real,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)
            i_real += 1 # If success, increase i_real -- this is a counter for the number of successful events
            continue

        # Delete the buffer files.
        for fname in [buffername,buffername_truth]:
            comm = ['rm', fname]
            try: sub.check_call(comm,stderr=sub.DEVNULL)
            except: pass

        return i_real-1, n_fail

    # Generate a bunch of events in the given pT range,
    # and save them to a HepMC file.
    # We do perform event selection: Only certain particles are saved to the file to begin with.
    def Generate(self,nevents): # TODO: implement file chunking
        # pythia = npyth.Pythia(config=self.pythia_config_file,params=self.pythia_config)
        filename = '{}/{}'.format(self.outdir,self.filename)

        # # Get the Fastjet banner out of the way
        # tmp = InitFastJet()
        # del tmp

        qu.printProgressBarColor(0,nevents, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # Loop in such a way as to guarantee that we get as many events as requested.
        # This logic is required as events could technically fail selections, e.g. not have the
        # requested truth particles (depends on requested truth particles & processes).
        n_success = 0
        n_fail = nevents
        nloops = 0
        while(n_fail > 0):
            n_success, n_fail = self.GenerationLoop(nevents-n_success, filename,
                                            i_real=n_success+1, nevents_disp = nevents, loop_number = nloops)
            nloops = nloops + 1

        # if(self.diagnostic_plots): self.OutputHistograms()
        return

    # def HistogramFinalStateCodes(self,arr):
    #     codes = []
    #     for entry in arr:
    #         codes.append(entry[-2])
    #     FillPdgHist(self.hist_fs_pdgid,codes)
    #     return