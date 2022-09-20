import sys, glob
import numpy as np
import ROOT as rt
import subprocess as sub
import pyhepmc_ng as hep
import numpythia as npyth # Pythia, hepmc_write
from util.config import GetEventSelection, GetTruthSelection, GetFinalStateSelection, GetPythiaConfig, GetPythiaConfigFile, GetJetConfig, GetNPars
from util.qol_utils.pdg import pdg_names, pdg_plotcodes
from util.gen_utils.utils import CreateHepMCEvent, HepMCOutput, HistogramFinalStateKinematics, HistogramFinalStateCodes
from util.conv_utils.utils import InitFastJet
import util.qol_utils.qol_util as qu

class Generator:
    def __init__(self, pt_min, pt_max):
        self.pt_min = pt_min
        self.pt_max = pt_max

        self.ConfigPythia()
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
        self.diagnostic_plots = True
        self.InitializeHistograms()

    def SetOutputDirectory(self,dir):
        self.outdir = dir

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

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

            # # Debug: Print full event record. #TODO: This is useful to keep here, is there a nicer way to implement this?
            # all_particles = event.all()
            # print('\n---START---')
            # for par in all_particles:
            #     print('\t',par)
            # print('\n--- END ---')

            # Get the truth-level particles.
            try: arr_truth = self.truth_selection(event=event)
            except: arr_truth = np.zeros((0,4))

            if(len(arr_truth) != self.n_truth and self.truth_selection is not None): # missing some desired truth particle -> potential trouble, discard this event.
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

            if(self.diagnostic_plots):
                # Record the PDG codes of all selected final-state particles (i.e. all those that made it into the HepMC output).
                HistogramFinalStateCodes(arr,self.hist_fs_pdgid)
                HistogramFinalStateKinematics(arr,self.kinematic_hists)

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

        filename = '{}/{}'.format(self.outdir,filename)

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

        if(self.diagnostic_plots): self.OutputHistograms()
        return

    # def HistogramFinalStateCodes(self,arr):
    #     codes = []
    #     for entry in arr:
    #         codes.append(entry[-2])
    #     FillPdgHist(self.hist_fs_pdgid,codes)
    #     return