import ROOT as rt
import numpy as np
from numpy.typing import NDArray
from typing import TYPE_CHECKING, Optional,Union,List,Annotated,Any

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class RhoAreaSubtraction:
    """
    Performs "rho-area subtraction", to remove pileup from jets.
    Rho is an estimate of the pileup activity/density; it can be
    provided as a direct input (e.g. Delphes has a module for
    computing it), or it can be computed on-the-fly.

    Depending on the mode,

    Note that jet vectors are affected, but jet constituents are not.
    """

    def __init__(self,rho_input:Optional[str]=None,pt_min:Annotated[float,"GeV"]=20., rho_eta_edges:Optional[Union[List[Any],NDArray]]=None, eta_bins:Optional[Union[list,NDArray]]=None, phi_bins:Optional[Union[list,NDArray]]=None,mode='decorate'):

        self.mode = mode
        assert mode in ['overwrite','decorate']

        # Print a warning if the mode is "overwrite", since this
        # is sort of non-standard; it'll actually modify jet momenta,
        # but not the constituent momenta.
        if(mode=='overwrite'):
            self._print('Warning: Using \'overwrite\' mode, will replace')
            self._print('raw jet momenta with rho-area-corrected momenta,')
            self._print('but will *not* affect jet constituents.')
            self._print('Constituent momenta may no longer sum up to jet momentum!')
            print()

        # Variables for reading in rho input,
        # for pre-computed rho (from Delphes).
        self.rho_input = rho_input
        self.pt_min = pt_min
        self.rho_name = None
        self.rho_edge_name = None
        self.compute = False

        # Variables for computing rho.
        self.eta_bins = eta_bins
        self.phi_bins = phi_bins
        self.estimators = None

        # Variables for saving to new branches
        self.pmu = None
        self.pmu_cyl = None
        self.pmu_name = None
        self.pmu_cyl_name = None

        # transient storage
        self.rho = None
        self.rho_eta_edges = rho_eta_edges # edges in pseudo-rapidity -- transient only in case of compute==False

        self.name = 'RhoAreaSubtraction'
        self.branch_name = None
        self.print_prefix = '\n\t{}'.format(self.name)
        self.citations = {
            "Rho-Area Subtraction":
            """
@article{Cacciari:2007fd,
    author = "Cacciari, Matteo and Salam, Gavin P.",
    title = "{Pileup subtraction using jet areas}",
    eprint = "0707.1378",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "LPTHE-07-01",
    doi = "10.1016/j.physletb.2007.09.077",
    journal = "Phys. Lett. B",
    volume = "659",
    pages = "119--126",
    year = "2008"
}
            """
        }
    def GetCitations(self):
        return self.citations

    def ModifyInitialization(self,obj):

        if(self.rho_input is not None):
            self._set_rho_input(obj)
            self.compute = False
        else:
            self._init_estimator()
            self.compute = True # will do on-the-fly computation
        return

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        Performs rho-area subtraction on the jets.

        If mode=='overwrite', modifies them in-place
        (adjusts jet momentum, *not* constituents)!

        If mode=='decorate', writes the rho-area-corrected
        momentum to a new set of branches.
        """

        import fastjet as fj # can do this since JetFinder ensures this is set up and cached correctly

        # fetch/compute rho
        if(self.compute):
            self._compute_rho() # fills self.rho (self.rho_eta_edges must be already filled in this case)
        else:
            self._fetch_rho() # fills self.rho, self.rho_eta_edges

        delete_indices = [] # keep track of jets to remove entirely (or assign 0 corrected momentum)

        self.pmu = {i:np.zeros(4) for i in obj.jets_dict.keys()}
        self.pmu_cyl = {i:np.zeros(4) for i in obj.jets_dict.keys()}

        for i,jet in obj.jets_dict.items():

            # get the jet area
            area = None
            if(jet.has_area()):

                # Get the area 4-vector
                area = jet.area_4vector() # pseudojet

            if(area is None):
                continue

            # determine the appropriate value of rho to use, based on jet's location
            rho = None
            rho_idx = np.searchsorted(self.rho_eta_edges, jet.eta(), side='right') - 1
            if(0 <= rho_idx < len(self.rho)):
                rho = self.rho[rho_idx]

            if(rho is None):
                continue

            rho_area = rho * area # Fastjet::PseudoJet

            if(rho_area.pt() >= jet.pt()):
                delete_indices.append(i)
                continue

            if(self.mode=='overwrite'):
                obj.jets_dict[i] -= rho_area
                if(obj.jets_dict[i].pt() < self.pt_min):
                    delete_indices.append(i)
                    continue

            else:
                corr_jet = jet - rho_area
                if(corr_jet.pt() < self.pt_min):
                    delete_indices.append(i)
                    continue
                self.pmu[i]     = np.array([jet.e(),jet.px(),jet.py(),jet.pz()])
                self.pmu_cyl[i] = np.array([jet.pt(),jet.eta(),jet.phi(),jet.m()])

        if(self.mode=='overwrite'):
            # Now, drop any jets that failed the minimum pt check
            # or were totally removed by rho-area subtraction
            for idx in delete_indices:
                del obj.jets_dict[idx]
            self._jetsToVectors() # force computation of jet vectors again (pseudojets -> numpy arrays)

        return

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='decorate'):
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addToBuffer(obj)
        return

    def ModifyConstituents(self, obj : 'JetFinder'):
        # TODO figure out what, if anything, to do here
        return

    def _set_rho_input(self,obj : 'JetFinder'):
        self.rho_name = '{}.Rho'.format(self.rho_input)
        self.rho_edge_name = '{}.Edges'.format(self.rho_input)

        for name in [self.rho_name, self.rho_edge_name]:
            if(name not in obj.input_collection_arrays.keys()):
                obj.input_buffer.read_branch(name)
        return

    def _fetch_rho(self,obj : 'JetFinder'):
        self.rho = np.array(obj.input_buffer[self.rho_name])
        self.rho_eta_edges = np.array(obj.input_buffer[self.rho_eta_edges])

    def _compute_rho(self,obj):
        import fastjet as fj # can do this since JetFinder has imported for us

        # Turn the JetFinder inputs into pseudojets, and pass them to the rho estimators.
        # Each estimator covers a different rapidity range, but it's OK to hand each all the inputs.
        # We'll actually fetch the pseudojet objects from JetFinderBase; we shouldn't need
        # to worry about indexing (it doesn't matter if these ended up in jets or not).
        pseudojets = obj.pseudojets[:len(obj.input_vecs)] # list of Fastjet::PseudoJet -- it's a buffer so truncate to correct length (could have unused extra pseudojets)
        self.rho = np.zeros(len(self.estimators))
        for i,estimator in enumerate(self.estimators):
            estimator.set_particles(pseudojets)
            self.rho[i] = estimator.rho()
        return

    def _init_estimator(self):
        import fastjet as fj # can do this since JetFinder has imported for us

        n = len(self.rho_eta_edges) - 1
        self.eta_bins = np.array(self.eta_bins)
        self.phi_bins = np.array(self.phi_bins)
        # TODO: Explicitly handle broadcasting of eta/phi bins if of length 1

        self.estimators = []
        for i in range(n):
            self.estimators.append(
                fj.GridMedianBackgroundEstimator(
                    self.rho_eta_edges[i],
                    self.rho_eta_edges[i+1],
                    self.eta_bins[i],
                    self.phi_bins[i]
                )
            )
        return

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Adds branches to the buffer corresponding with rho-area-corrected
        jet four-momentum.
        """
        self._createBranchNames(obj)

        if(self.pmu_name not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.pmu_name,ndim=2,dtype=np.dtype('f8'))

        if(self.pmu_cyl_name not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.pmu_cyl_name,ndim=2,dtype=np.dtype('f8'))
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        self.branch_name_prefixname = '{}.{}'.format(obj.jet_name,self.name)

        self.pmu_name = '{}.Pmu'.format(self.branch_name_prefixname)
        self.pmu_cyl_name = '{}.Pmu_cyl'.format(self.branch_name_prefixname)
        return

    def _addToBuffer(self,obj : 'JetFinder'):
        """
        Adds the corrected jet momentum to the buffer.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.output_buffer.set(self.pmu_name,obj._i,self.pmu)
        obj.output_buffer.set(self.pmu_cyl_name,obj._i,self.pmu_cyl)

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return