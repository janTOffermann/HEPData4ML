# Some functions for filtering out events during generation -- if an event fails to pass some condition
# on its selected truth and final-state particles (the inputs to jet clustering), the event is thrown out
# and not written to the HepMC file -- is not counted towards the total number of events we have generated.
import sys
import numpy as np
from util.calcs import Calculator
from util.fastjet import FastJetSetup

class FilterFlag:
    """
    A wrapper class for a filter. This is used for the "event_filter_flag" configuration option,
    which runs an event filter -- but instead of literally filtering out events that fail the filter,
    it creates a boolean flag that will be added to the final dataset. This is achieved by producing
    a special HDF5 file during the generation step, as that is the step where the event filter algo is run.
    """
    def __init__(self,filter,name='event_filter'):
        self.filter = filter
        self.name = name

    def __call__(self,pythia_wrapper,final_state_indices):
        return self.filter(pythia_wrapper,final_state_indices)

    def GetName(self):
        return self.name

    def Initialize(self,configurator):
        self.filter.Initialize(configurator)
        return

class NotFilter:
    """
    A class for making an inverse of an event filter.
    """
    def __init__(self,filter):
        self.filter = filter

    def __call__(self,pythia_wrapper,final_state_indices):
        return (not filter(pythia_wrapper,final_state_indices))

    def Initialize(self,configurator):
        self.filter.Initialize(configurator)
        return

class MultiFilter:
    """
    A class for combining multiple event filters.
    """
    def __init__(self,filters=[]):
        self.filters = filters

    def __call__(self,pythia_wrapper,final_state_indices):
        for filter in self.filters:
            individual_status = filter(pythia_wrapper,final_state_indices)
            if(not individual_status): return False
        return True

    def Initialize(self,configurator):
        for filter in self.filters:
            filter.Initialize(configurator)
        return

class DefaultFilter:
    """
    A class for making an event filter
    that simply defaults to True/False.
    """
    def __init__(self,value=True):
        self.value = value

    def __call__(self,pythia_wrapper,final_state_indices):
        return self.value

    def Initialize(self,configurator):
        pass

class PtJetFilter:
    """
    A filter that requires an event to have at least
    n jets (truth-level) with pT within some range.
    Defaults to just using a lower threshold.
    """
    def __init__(self,jet_radius,njet,pt_min_jet=35., pt_max_jet=None, eta_max_jet=None):
        self.SetJetRadius(jet_radius)
        self.pt_min_jet = pt_min_jet
        self.pt_max_jet = pt_max_jet
        self.eta_max_jet = eta_max_jet
        self.njet = njet

        self.configurator = None

    def Initialize(self,configurator):
        self.configurator = configurator
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # this import should work, will be others peppered throughout for scope reasons but Python caches imports!
        self.jetdef = fj.JetDefinition(fj.antikt_algorithm, self.radius)
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetJetRadius(self,radius):
            self.radius = radius

    def __call__(self,pythia_wrapper,final_state_indices):
        import fastjet as fj
        # Get the four-vectors of the final-state particles.
        # Using the signature (px, py, pz, E) since this is what fastjet uses
        fs_particles = pythia_wrapper.GetPxPyPzE(final_state_indices)

        # Perform jet clustering.
        pseudojets = [fj.PseudoJet(*x) for x in fs_particles]
        cluster_sequence = fj.ClusterSequence(pseudojets, self.jetdef)
        jets = cluster_sequence.inclusive_jets()

        # Apply optional eta window cut to the jets we consider.
        if(self.eta_max_jet is not None):
            jets_eta = np.array([jet.eta() for jet in jets])
            jets = jets[np.where(np.abs(jets_eta) < self.eta_max_jet)[0]]

        if(len(jets) < self.njet): return False # not enough jets -> already failed

        # Check if there are sufficient jets with pT above/below thresholds.
        jets_pt = np.array([jet.pt() for jet in jets],dtype=float)
        if(self.pt_min_jet is not None): jets = jets[np.where(jets_pt > self.pt_min_jet)[0]]
        if(self.pt_max_jet is not None): jets = jets[np.where(jets_pt < self.pt_max_jet)[0]]
        if(len(jets) < self.njet): return False
        return True

class MassJetFilter:
    """
    A filter that requires an event to have at least
    n jets (truth-level) with mass within some range.
    Defaults to just using a lower threshold.
    """
    def __init__(self,jet_radius,njet,mass_min_jet=35., mass_max_jet=None, eta_max_jet=None):
        self.SetJetRadius(jet_radius)
        self.mass_min_jet = mass_min_jet
        self.mass_max_jet = mass_max_jet
        self.eta_max_jet = eta_max_jet
        self.njet = njet

        self.configurator = None

    def Initialize(self,configurator):
        self.configurator = configurator
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # this import should work, will be others peppered throughout for scope reasons but Python caches imports!
        self.jetdef = fj.JetDefinition(fj.antikt_algorithm, self.radius)
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetJetRadius(self,radius):
            self.radius = radius

    def __call__(self,pythia_wrapper,final_state_indices):
        import fastjet as fj
        # Get the four-vectors of the final-state particles.
        # Using the signature (px, py, pz, E) since this is what fastjet uses
        fs_particles = pythia_wrapper.GetPxPyPzE(final_state_indices)

        # Perform jet clustering.
        pseudojets = [fj.PseudoJet(*x) for x in fs_particles]
        cluster_sequence = fj.ClusterSequence(pseudojets, self.jetdef)
        jets = cluster_sequence.inclusive_jets()

        # Apply optional eta window cut to the jets we consider.
        if(self.eta_max_jet is not None):
            jets_eta = np.array([jet.eta() for jet in jets])
            jets = jets[np.where(np.abs(jets_eta) < self.eta_max_jet)[0]]

        if(len(jets) < self.njet): return False # not enough jets -> already failed

        # Check if there are sufficient jets with mass above/below thresholds.
        jets_mass = np.array([jet.m() for jet in jets],dtype=float)
        if(self.mass_min_jet is not None): jets = jets[np.where(jets_mass > self.mass_min_jet)[0]]
        if(self.mass_max_jet is not None): jets = jets[np.where(jets_mass < self.mass_max_jet)[0]]
        if(len(jets) < self.njet): return False
        return True

class PtMatchedJetFilter:
    """
    A multi-step filter:
    - Cluster jets in the event (truth-level, i.e. no DELPHES).
    - Require there to be a jet within some dR of a particular truth particle (identified by its index in the truth selection).
    - Require that jet's pT to be within some range of that truth particle's pT.
    """
    def __init__(self,jet_radius,truth_selector,matching_radius=None, pt_window=None, pt_window_frac=None, pt_min_jet=None, eta_max_jet=None):
        self.SetJetRadius(jet_radius)
        self.SetTruthSelector(truth_selector) # index of truth particle we want the jet to be close to
        if(matching_radius is None): matching_radius = jet_radius
        if(pt_window is None): pt_window = 0.
        self.pt_window = pt_window
        self.pt_window_frac_mode = False
        self.SetPtWindowFraction(pt_window_frac)
        self.configurator = None

        self.SetMatchingRadius(matching_radius)
        if(pt_min_jet) is None: pt_min_jet = 0.
        self.pt_min_jet = pt_min_jet
        self.eta_max_jet = eta_max_jet

    def Initialize(self,configurator):
        self.configurator = configurator
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # this import should work, will be others peppered throughout for scope reasons but Python caches imports!
        self.jetdef = fj.JetDefinition(fj.antikt_algorithm, self.radius)
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetJetRadius(self,radius):
        self.radius = radius

    def SetMatchingRadius(self,radius):
        self.matching_radius = radius

    def SetTruthSelector(self,selector):
        self.truth_selector = selector

    def SetPtWindowFraction(self,val):
        self.pt_window_frac = val
        if(val is not None): self.pt_window_frac_mode = True

    def __call__(self,pythia_wrapper,final_state_indices):
        import fastjet as fj
        # Get the truth particle from our truth selector.
        truth_indices = np.atleast_1d(np.array(self.truth_selector(pythia_wrapper),dtype=int))
        truth_selection_status = self.truth_selector.GetSelectionStatus()
        if(not truth_selection_status): return False
        self.truth_index = truth_indices[0] # assume we are only using a selector giving 1 output, otherwise just take the first one

        # Get the four-vectors of the final-state particles.
        # Using the signature (px, py, pz, E) since this is what fastjet uses
        fs_particles = pythia_wrapper.GetPxPyPzE(final_state_indices)

        # Perform jet clustering.
        pseudojets = [fj.PseudoJet(*x) for x in fs_particles]
        cluster_sequence = fj.ClusterSequence(pseudojets, self.jetdef)
        jets = cluster_sequence.inclusive_jets()

        if(len(jets) == 0): return False # no jets -> already failed

        # Select a jet within matching_radius of the truth particle with truth_index
        truth_particle_eta_phi = np.array([pythia_wrapper.GetEta(self.truth_index), pythia_wrapper.GetPhi(self.truth_index)])
        truth_pt = pythia_wrapper.GetPt(self.truth_index)
        jets_eta_phi = np.array([[jet.eta(),jet.phi()] for jet in jets])

        d2 = np.array([self.calculator.DeltaR2(*truth_particle_eta_phi,*j) for j in jets_eta_phi])
        amin = np.argmin(d2)
        d2_min = d2[amin]
        if(d2_min > self.matching_radius * self.matching_radius): return False
        jet = jets[amin]
        jet_pt = jet.pt()
        jet_eta = jet.eta()

        if(jet_pt < self.pt_min_jet): return False
        if(self.eta_max_jet is not None):
            if(np.abs(jet_eta) > np.abs(self.eta_max_jet)): return False

        # Now, compare the selected jet's pT to the truth particle's pT.
        if(self.pt_window_frac_mode):
            if(jet_pt > (1. - self.pt_window_frac) * truth_pt and jet_pt < (1. + self.pt_window_frac) * truth_pt): return True
            return False

        if(np.abs(jet_pt - truth_pt) > self.pt_window): return False
        return True

class ContainedJetFilter:
    """
    A multi-step filter:
    - Cluster jets in the event (truth-level, i.e. no DELPHES).
    - Require there to be a jet within some dR of a particular truth particle (identified by its index in the truth selection).
    - Require all the stable daughters of a particular truth particle -- above some pT threshold -- to be within that jet's radius.
    """
    def __init__(self,jet_radius,truth_selector, daughter_selector,matching_radius=None, pt_threshold=None, pt_threshold_sum=None, pt_threshold_frac=None, pt_min_jet=None, eta_max_jet=None):
        self.SetJetRadius(jet_radius)
        self.SetTruthSelector(truth_selector) # index of truth particle we want the jet to be close to
        if(matching_radius is None): matching_radius = jet_radius
        if(pt_threshold is None): pt_threshold = 0.
        self.SetDaughterPtThreshold(pt_threshold) # daughter particles with pT below this value (GeV) are allowed to not be contained
        self.pt_threshold_frac_mode = False
        self.pt_threshold_sum_mode = False
        self.SetDaughterPtThresholdFraction(pt_threshold_frac)
        self.SetDaughterPtThresholdSum(pt_threshold_sum)
        self.configurator = None

        if(self.pt_threshold != 0. and self.pt_threshold_frac_mode):
            print('\n\t: Warning: ContainedJetFilter has pt_min != 0 and pt_min_frac != 0. Using fractional pT mode.\n')

        self.SetMatchingRadius(matching_radius)
        self.SetDaughterSelector(daughter_selector)
        self.jetdef = None
        if(pt_min_jet) is None: pt_min_jet = 0.
        self.pt_min_jet = pt_min_jet
        self.eta_max_jet = eta_max_jet

    def Initialize(self,configurator):
        self.configurator = configurator
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # this import should work, will be others peppered throughout for scope reasons but Python caches imports!
        self.jetdef = fj.JetDefinition(fj.antikt_algorithm, self.radius)
        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetJetRadius(self,radius):
        self.radius = radius

    def SetMatchingRadius(self,radius):
        self.matching_radius = radius

    def SetTruthSelector(self,selector):
        self.truth_selector = selector

    def SetDaughterSelector(self,selector):
        self.daughter_selector = selector

    def SetDaughterPtThreshold(self,val):
        self.pt_threshold = val

    def SetDaughterPtThresholdFraction(self,val):
        self.pt_threshold_frac = val
        if(val is not None):
            self.pt_threshold_frac_mode = True
            self.pt_threshold_sum_mode = False

    def SetDaughterPtThresholdSum(self,val):
        self.pt_threshold_sum = val
        if(val is not None):
            self.pt_threshold_sum_mode = True
            self.pt_threshold_frac_mode = False

    def __call__(self,pythia_wrapper,final_state_indices):
        import fastjet as fj
        # Get the truth particle from our truth selector.
        truth_indices = np.atleast_1d(np.array(self.truth_selector(pythia_wrapper),dtype=int))
        truth_selection_status = self.truth_selector.GetSelectionStatus()
        if(not truth_selection_status): return False
        self.truth_index = truth_indices[0] # assume we are only using a selector giving 1 output, otherwise just take the first one

        # Get the four-vectors of the final-state particles.
        # Using the signature (px, py, pz, E) since this is what fastjet uses
        fs_particles = pythia_wrapper.GetPxPyPzE(final_state_indices)

        # Perform jet clustering.
        pseudojets = [fj.PseudoJet(*x) for x in fs_particles]
        cluster_sequence = fj.ClusterSequence(pseudojets, self.jetdef)
        jets = cluster_sequence.inclusive_jets()

        if(len(jets) == 0): return False # no jets -> already failed

        # Select a jet within matching_radius of the truth particle with truth_index
        truth_particle_eta_phi = np.array([pythia_wrapper.GetEta(self.truth_index), pythia_wrapper.GetPhi(self.truth_index)])
        jets_eta_phi = np.array([[jet.eta(),jet.phi()] for jet in jets])

        d2 = np.array([self.calculator.DeltaR2(*truth_particle_eta_phi,*j) for j in jets_eta_phi])
        amin = np.argmin(d2)
        d2_min = d2[amin]
        if(d2_min > self.matching_radius * self.matching_radius): return False
        jet = jets[amin]
        jet_pt = jet.pt()
        jet_eta = jet.eta()
        jet_phi = jet.phi()

        if(jet_pt < self.pt_min_jet): return False
        if(self.eta_max_jet is not None):
            if(np.abs(jet_eta) > np.abs(self.eta_max_jet)): return False

        # Now get the particles from our daughter selector. It's up to the user to make sure that this selector's definition
        # makes sense for their use-case!
        daughter_particles = np.atleast_1d(np.array(self.daughter_selector(pythia_wrapper),dtype=int))
        if(len(daughter_particles) == 0): return False

        # Get the distances between the selected jet and each daughter particle.
        daughters_eta_phi = np.array([[pythia_wrapper.GetEta(x),pythia_wrapper.GetPhi(x)] for x in daughter_particles])
        daughters_pt = np.array([pythia_wrapper.GetPt(x) for x in daughter_particles])

        # If we use the "sum mode", we check the pT of the sum of uncontained daughters versus the threshold.
        if(self.pt_threshold_sum_mode):
            d2 = np.array([self.calculator.DeltaR2(*daughter,jet_eta,jet_phi) for daughter in daughters_eta_phi])
            mask = d2 > np.square(self.matching_radius)
            uncontained_daughters = daughter_particles[mask]

            if(len(uncontained_daughters) == 0): return True

            else:
                # do a 4-vector sum -- just need to get the pT so we only need total px and py
                uncontained_sum = np.sum(np.array([pythia_wrapper.GetPxPyPzE(x) for x in uncontained_daughters]),axis=0)
                uncontained_pt = np.linalg.norm(uncontained_sum[:2],ord=2)

            if(uncontained_pt > self.pt_threshold_sum): return False
            return True

        if(self.pt_threshold_frac_mode):
            daughter_sum = np.sum(np.array([pythia_wrapper.GetPxPyPzE(x) for x in daughter_particles]),axis=0)
            self.pt_threshold = np.sqrt(daughter_sum[0] * daughter_sum[0] + daughter_sum[1] * daughter_sum[1]) * self.pt_threshold_frac

        mask = daughters_pt >= self.pt_threshold
        daughters_eta_phi = daughters_eta_phi[mask]
        if(len(daughters_eta_phi) == 0): return False

        d2 = np.array([self.calculator.DeltaR2(*daughter,jet_eta,jet_phi) for daughter in daughters_eta_phi])
        d2_max = np.max(d2)
        if(d2_max > self.matching_radius * self.matching_radius): return False
        return True