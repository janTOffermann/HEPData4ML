import sys,os,glob
import ROOT as rt
import numpy as np
import pyhepmc as hep
import pythia8 as pyth8

class SimplePileup:

    def __init__(self, pileup_hepmc3,rng_seed=1,mu_file=None):
        self.files = sorted(glob.glob(pileup_hepmc3))

        # determine the total number of events in these files
        self.n_dict = {}
        self.first_idx = {}
        self.n_total = 0
        for file in self.files:
            with open(file,'r') as f:
                n =  len([x for x in f.readlines() if 'E ' in x])
                self.first_idx[file] = self.n_total
                self.n_dict[file] = n
                self.n_total += n

        self.event_indices = np.arange(self.n_total)
        self.event_mask = np.full(self.n_total,True,dtype=bool)
        self.selected_indices = None # temporary storage for indices selected for a particular event

        self.rng_seed = rng_seed
        self.rng = np.random.default_rng(self.rng_seed)

        # Set up the distribution for # of interactions per crossing
        self.mu_file = mu_file
        self.InitializeMuDistibution()
        self.mu = None # store the current value of mu

    def InitializeMuDistibution(self):
        hist_name = 'SimplePileup_mu'

        if(self.mu_file is None):
            # default
            nbins = 80
            self.mu_distribution = rt.TH1F(hist_name,'',nbins,0,nbins)
            # very approximate for Run 2,
            # see https://atlas.web.cern.ch/Atlas/GROUPS/DATAPREPARATION/PublicPlots/2018/DataSummary/figs/mu_2015_2018.png
            mu = 33.7
            sigma = 11.5
            for i in range(nbins):
                self.mu_distribution.SetBinContent(i+1,self._gaussian(i,mu,sigma))

        else:
            try:
                f = rt.TFile(self.mu_file,"READ")
                self.mu_distribution = f.Get(hist_name).Clone()
                self.mu_distribution.SetDirectory(0)
                f.Close()
            except:
                print("Error: SimplePileup failed to open file {}.".format(self.mu_file))
                print("Falling back on default mu distribution.")
                self.mu_file = None
                self.InitializeMuDistibution()

    def SampleMuDistribution(self):
        n = self.mu_distribution.GetNbinsX()
        mu_min = self.mu_distribution.GetXaxis().GetBinLowEdge(1)
        mu_max = self.mu_distribution.GetXaxis().GetBinLowEdge(n+1)
        mu_values = np.arange(mu_min,mu_max)
        mu_probs = np.array([self.mu_distribution.GetBinContent(i+1) for i in range(n)],dtype=float)
        mu_probs /= np.sum(mu_probs)
        self.mu = int(self.rng.choice(mu_values,p=mu_probs))

    def PickEventIndices(self, update_mask=False):
        self.SampleMuDistribution()
        print('mu = {}'.format(self.mu))
        self.selected_indices = self.rng.choice(self.event_indices[self.event_mask],size=self.mu)
        if(update_mask): self.UpdateMask()

    def UpdateMask(self,indices=None):
        """
        Remove "indices" from pileup event indices that can be sampled from in the future.
        Used to achieve sampling without replacement, but having this function separated
        allows us to reuse the indices in some cases (e.g. pileup events were selected
        but that event was then thrown out by some generation-level cut).
        """
        if(indices is None): indices = self.selected_indices # use currently selected indices
        self.event_mask[self.selected_indices] = False

    def FetchEvent(self,idx):
        file_idx = np.where(idx < np.array(list(self.first_idx.values()),dtype=int))[0][0]
        filename = self.files[file_idx]
        idx_in_file = idx - self.first_idx[filename]

        # TODO: Would be nice to find alternative, but for now we have to iterate through
        #       the given HepMC3 file since pyhepmc treats files as iterable but does
        #       not have any way to just fetch the i'th event. This will scale poorly
        #       with large HepMC3 files.
        event = 0
        with hep.io.ReaderAscii(filename) as f:
            for i,evt in enumerate(f):
                while(i != idx_in_file): continue
                event = evt

        print(event)

    def __call__(self,pythia_wrapper):
        print("Called SimplePileup.")
        # determine indices of events to add
        self.PickEventIndices()

        for i,idx in enumerate(self.selected_indices):
            if(i > 0): break
            self.FetchEvent(idx)


    def _gaussian(self,x,mu,sig,A=None):
        if(A is None): A = 1. / (np.sqrt(2.0 * np.pi))
        return A * np.exp(-np.power((x - mu) / sig, 2.0) / 2)





