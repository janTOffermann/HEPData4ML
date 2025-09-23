# Some functions for filtering out events during generation -- if an event fails to pass some condition
# it can be thrown out instead of being written to the HepMC file -- and it is not counted towards the
# total number of events we have generated.
import numpy as np
from util.reconstruction.post_processing.jets import TruthJetFinder
from typing import TYPE_CHECKING

if(TYPE_CHECKING):
    from util.config.config import Configurator
    from util.pythia.pythia import PythiaPythonWrapper

class DefaultFilter:
    """
    A class for making an event filter
    that simply defaults to True/False.
    """
    def __init__(self,value=True):
        self.value = value

    def __call__(self,pythia_wrapper:'PythiaPythonWrapper'):
        return self.value

    def Initialize(self,configurator):
        pass

class NotFilter:
    """
    A class for making an inverse of an event filter.
    """
    def __init__(self,filter):
        self.filter = filter

    def __call__(self,pythia_wrapper:'PythiaPythonWrapper'):
        return (not filter(pythia_wrapper))

    def Initialize(self,configurator):
        self.filter.Initialize(configurator)
        return

class MultiFilter:
    """
    A class for combining multiple event filters; this checks that
    all filters are passed. Note that in practice, the ordering of
    the sequence can affect runtime, possibly if one is quicker to
    evaluate than another.
    """
    def __init__(self,filters=[]):
        self.filters = filters

    def __call__(self,pythia_wrapper:'PythiaPythonWrapper'):
        for filter in self.filters:
            individual_status = filter(pythia_wrapper)
            if(not individual_status): return False
        return True

    def Initialize(self,configurator):
        for filter in self.filters:
            filter.Initialize(configurator)
        return

class PtJetFilter:
    """
    A filter that requires an event to have at least a certain
    number of truth-level anti-kt jets with pT within some range.
    Defaults to just using a lower threshold.
    """
    def __init__(self,jet_radius,njet,pt_min_jet=35., pt_max_jet=None, eta_max_jet=None):
        self.SetJetRadius(jet_radius)
        self.pt_min_jet = pt_min_jet
        self.pt_max_jet = pt_max_jet
        self.eta_max_jet = eta_max_jet
        self.njet = njet
        self.configurator = None

        self.jet_finder = TruthJetFinder(radius=jet_radius)

    def Initialize(self,configurator:'Configurator'):
        self.configurator = configurator
        self.jet_finder.SetConfigurator(configurator)
        self.jet_finder.Initialize()

    def SetJetRadius(self,radius):
            self.radius = radius
            self.jet_finder.SetRadius(self.radius)

    def __call__(self,pythia_wrapper:'PythiaPythonWrapper'):

        # Get the four-vectors of the visible final-state particles.
        # Using the signature (E, px, py, pz).
        stable_indices = np.where(pythia_wrapper.GetStatus(hepmc=True)==1)[0]
        visible_indices =  np.where(~np.isin(np.abs(pythia_wrapper.GetPdgId()), [12, 14, 16]))[0]
        stable_indices = np.intersect1d(stable_indices,visible_indices,assume_unique=True)
        input_vecs = pythia_wrapper.GetEPxPyPz(stable_indices)

        self.jet_finder.Process(input_vecs)
        jet_vecs_cyl = self.jet_finder.jet_vectors_cyl

        mask = np.full(len(jet_vecs_cyl),True)

        # Apply optional eta window cut to the jets we consider.
        if(self.eta_max_jet is not None):
            mask *= np.where(np.abs(jet_vecs_cyl[:,1]) < self.eta_max_jet)[0]
        if(np.sum(mask) < self.njet): return False

        if(self.pt_min_jet is not None):
            mask *= np.where(jet_vecs_cyl[:,0] > self.pt_min_jet)[0]
        if(np.sum(mask) < self.njet): return False

        if(self.pt_max_jet is not None):
            mask *= np.where(jet_vecs_cyl[:,0] < self.pt_max_jet)[0]
        if(np.sum(mask) < self.njet): return False
        return True