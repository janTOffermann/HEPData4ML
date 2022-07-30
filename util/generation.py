import sys, glob
import numpy as np
import ROOT as rt
import subprocess as sub
import pyhepmc_ng as hep
import numpythia as npyth # Pythia, hepmc_write
from util.config import GetEventSelection, GetTruthSelection, GetFinalStateSelection, GetPythiaConfig, GetPythiaConfigFile, GetJetConfig
# from util.fastjet import BuildFastjet
from util.gen_utils.utils import CreateHepMCEvent, HepMCOutput
from util.conv_utils.utils import InitFastJet
from util.hepmc import RestructureParticleArray
import util.qol_utils.qol_util as qu

# # --- FASTJET IMPORT ---
# # TODO: Can this be done more nicely?
# fastjet_dir = BuildFastjet(j=8)
# fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
# if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
# import fastjet as fj
# # ----------------------

class Generator:
    def __init__(self, pt_min, pt_max):
        self.pt_min = pt_min
        self.pt_max = pt_max

        self.ConfigPythia()
        self.event_selection = GetEventSelection()
        self.truth_selection = GetTruthSelection()
        self.final_state_selection = GetFinalStateSelection()
        self.jet_config = GetJetConfig()

        # Things for the progress bar.
        self.prefix = 'Generating events for pT bin [{},{}]:'.format(self.pt_min,self.pt_max)
        self.suffix = 'Complete'
        self.bl = 50

        self.outdir = None # TODO: Use this with filenames.

        # TODO: Add some histograms.
        self.InitializeHistograms()

    def SetOutputDirectory(self,dir):
        self.outdir = dir

    def ConfigPythia(self):
        self.pythia_config = GetPythiaConfig(self.pt_min,self.pt_max)
        self.pythia_config = self.PythiaCompatibilityCheck(self.pythia_config)
        self.pythia_config_file = GetPythiaConfigFile()

    # TODO: This is a compatibility check for the pythia configuration.
    # Numpythia v. 1.2 uses Pythia 8.244, but Numpythia v. 1.1 uses Pythia 8.226.
    # There are some noticeable differences in available process configs, and the latest numpythia
    # may not be available on macOS/arm yet.
    # The check is not very efficient, it will make many nympythia objects.
    def PythiaCompatibilityCheck(self,pythia_config):
        pythia_config_safe = {}
        safe_keys = ['Beam','Random','Stat','Print','PhaseSpace','PartonLevel','HadronLevel']
        for key,val in pythia_config.items():
            safe = False
            for safe_key in safe_keys:
                if(safe_key in key):
                    safe = True
                    break
            if(safe):
                pythia_config_safe[key] = val
                continue

            tmp_config = {key_safe:val_safe for key_safe,val_safe in pythia_config_safe.items()}
            tmp_config[key] = val

            with qu.stdout_redirected():
                try:
                    pythia = npyth.Pythia(params=tmp_config)
                    pythia_config_safe[key] = val
                    del pythia # TODO: Is this helpful?
                except: pass
        return pythia_config_safe

    def InitializeHistograms(self):
        self.hists = []
        self.hist_fs_pdgid = rt.TH1D(qu.RN(),'PDG codes (all final-state particles);PDG code;Count',1000,-500,500)
        self.hists.append(self.hist_fs_pdgid)

    def OutputHistograms(self,filetype='root'):
        for i,h in enumerate(self.hists):
            c = rt.TCanvas(qu.RN(),'',800,600)
            h.Draw('HIST')
            rt.gPad.SetLogy()
            c.Draw()
            c.SaveAs('{}/hist_{}.{}'.format(self.outdir,str(i).zfill(2),filetype))
            del c
        return

    def GenerationLoop(self,pythia, nevents,filename,i_real = 1, nevents_disp=None,loop_number=0):

        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that pyhepmc_ng's WriterAscii works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality.

        filename_truth = filename.replace('.hepmc','_truth.hepmc')
        buffername = filename.replace('.hepmc','_buffer.hepmc')
        buffername_truth = buffername.replace('.hepmc','_truth.hepmc')

        for i,event in enumerate(pythia(events=nevents)):
            success = True

            # Get the truth-level particles.
            arr_truth = self.truth_selection(event=event)

            for truth in arr_truth:
                if(len(truth) == 0): # missing some truth particle -> potential trouble?
                    n_fail += 1
                    success = False
                    break

            if(not success): continue

            # Get the final-state particles, as a numpy array.
            status, arr = self.final_state_selection(event=event)

            # If we didn't pick up any final-state particles, discard this event.
            if(not status):
                n_fail += 1
                success = False
                continue

            # Optionally filter down events (e.g. throw out particles too far from selected truth particles).
            # This is typically some sort of filtering done just to reduce the HepMC file size (which can be rather large).
            if(self.event_selection is not None):
                status, arr, arr_truth = self.event_selection(arr,arr_truth)

            if(not status):
                n_fail += 1
                success = False
                continue

            # self.HistogramFinalStateCodes(arr)

            # Create a GenEvent (pyhepmc_ng) containing these selected particles.
            # Note that the event comes from pyhepmc_ng, *not* numpythia.
            # Thus we cannot just use GenParticles from numpythia (e.g. using return_hepmc=True above).
            hepev       = CreateHepMCEvent(arr      ,i_real)
            hepev_truth = CreateHepMCEvent(arr_truth,i_real)

            # ----- File I/O -----
            # Write this event to the HepMC buffer file, then copy the buffer contents to the full HepMC file.
            HepMCOutput(hepev,buffername,filename,loop_number,i,i_real,nevents,n_fail)
            HepMCOutput(hepev_truth,buffername_truth,filename_truth,loop_number,i,i_real,nevents,n_fail)
            # ----- End File I/O -----

            qu.printProgressBarColor(i_real,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)
            i_real += 1

        # Delete the buffer files.
        for fname in [buffername,buffername_truth]:
            comm = ['rm', fname]
            try: sub.check_call(comm,stderr=sub.DEVNULL)
            except: pass

        return i_real-1, n_fail

    # Generate a bunch of events in the given pT range,
    # and save them to a HepMC file.
    # We do perform event selection: Only certain particles are saved to the file to begin with.
    def Generate(self,nevents, filename = 'events.hepmc'): # TODO: implement file chunking
        pythia = npyth.Pythia(config=self.pythia_config_file,params=self.pythia_config)

        # # Get our (Pythonic) Fastjet. # TODO: Make this optional (only need if *not* using Delphes)
        # fastjet_dir = BuildFastjet(j=8)
        # fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]

        # Get the Fastjet banner out of the way
        tmp = InitFastJet()
        del tmp

        qu.printProgressBarColor(0,nevents, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # Loop in such a way as to guarantee that we get as many events as requested.
        # This logic is required as events could technically fail selections, e.g. not have the
        # requested truth particles (depends on requested truth particles & processes).
        n_success = 0
        n_fail = nevents
        nloops = 0
        while(n_fail > 0):
            n_success, n_fail = self.GenerationLoop(pythia, nevents-n_success, filename,
                                            i_real=n_success+1, nevents_disp = nevents, loop_number = nloops)
            nloops = nloops + 1

        # self.OutputHistograms()
        return

    def HistogramFinalStateCodes(self,arr):
        for entry in arr:
            pdgid = entry[-2]
            self.hist_fs_pdgid.Fill(pdgid)
        return