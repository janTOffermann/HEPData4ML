import ROOT as rt
import numpy as np
from numpy.typing import NDArray
from typing import TYPE_CHECKING, Optional,Union,List,Annotated,Any
from util.reconstruction.post_processing.utils.postprocessor_base import PostProcessorBase

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class RhoAreaSubtraction(PostProcessorBase):
    """
    Performs "rho-area subtraction", to remove pileup from jets.
    Rho is an estimate of the pileup activity/density; it can be
    provided as a direct input (e.g. Delphes has a module for
    computing it), or it can be computed on-the-fly.

    Depending on the mode,

    Note that jet vectors are affected, but jet constituents are not.
    """

    def __init__(self,rho_input:Optional[str]=None,pt_min:Annotated[float,"GeV"]=20., rho_eta_edges:Optional[Union[List[Any],NDArray]]=None, eta_bins:Optional[Union[list,NDArray]]=None, phi_bins:Optional[Union[list,NDArray]]=None, save_rho_area:bool=False):
        super().__init__()

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

        self.save_rho_area = save_rho_area

        # Variables for saving to new branches
        self.buffers = {
            'Pmu':None,
            'Pmu_cyl':None,
            'RhoArea.Pmu':None,
            'RhoArea.Pmu_cyl':None,
            'Passed':None
        }
        self.branch_names = {key:None for key in self.buffers.keys()}

        # transient storage
        self.rho = None
        self.rho_eta_edges = rho_eta_edges # edges in pseudo-rapidity -- transient only in case of compute==False

        self.name = 'RhoAreaSubtraction'
        self.branch_name = None
        self.print_prefix = '{}'.format(self.name)
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

    def ModifyInitialization(self,obj):
        if(self.obj_name_input is None):
            self.SetInputObjectName(obj.jet_name)
        self.obj_name_output = '.'.join([self.obj_name_input,'RhoPileupSubtracted'])
        self.obj_name_constituents_output = self.obj_name_output

        obj.SetUseArea(True) # will need jet areas

        if(self.rho_input is not None):
            self._set_rho_input(obj)
            self.compute = False
        else:
            self._init_estimator()
            self.compute = True # will do on-the-fly computation

        self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
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
            self._compute_rho(obj) # fills self.rho (self.rho_eta_edges must be already filled in this case)
        else:
            self._fetch_rho(obj) # fills self.rho, self.rho_eta_edges

        pmu_key = '{}.Pmu'.format(self.obj_name_input)
        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)

        pmu_dict = obj.output_buffer.get(pmu_key,filter=obj.jet_ordering)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        # NOTE: This post-processor directly accesses the original jets, because it needs information
        #       on jet areas. We implicitly assume these are unchanged by any prior post-processors.
        # TODO: Consider storing jet area as a branch; then we can use it here w/out fastjet!

        delete_indices = [] # keep track of jets to remove entirely (or assign 0 corrected momentum)
        self.buffers['Pmu'] = {i:np.zeros(4) for i in pmu_cyl_dict.keys()}
        self.buffers['Pmu_cyl'] = {i:np.zeros(4) for i in pmu_cyl_dict.keys()}
        self.buffers['RhoArea.Pmu'] = {i:np.zeros(4) for i in pmu_cyl_dict.keys()}
        self.buffers['RhoArea.Pmu_cyl'] = {i:np.zeros(4) for i in pmu_cyl_dict.keys()}
        self.buffers['Passed'] = {i:False for i in pmu_cyl_dict.keys()}

        for i,pmu_cyl in pmu_cyl_dict.items():
            raw_jet = obj.jets_dict[i] # NOTE: we use the original jet from jet clustering here -- not affected by things like jet energy scale calibration! OK for accessing area.
            area = None
            if(raw_jet.has_area()):
                area = raw_jet.area_4vector() # 4-vector area, Fastjet::PseudoJet
            if(area is None):
                continue

            # determine the appropriate value of rho to use, based on jet's location
            rho = None
            rho_idx = np.searchsorted(self.rho_eta_edges[:,0], pmu_cyl[1], side='right') - 1
            if(0 <= rho_idx < len(self.rho)):
                rho = self.rho[rho_idx]

            if(rho is None):
                continue

            rho_area = rho * area # Fastjet::PseudoJet
            self.buffers['RhoArea.Pmu'][i]     = np.array([rho_area.e(),rho_area.px(),rho_area.py(),rho_area.pz()])
            self.buffers['RhoArea.Pmu_cyl'][i] = np.array([rho_area.pt(),rho_area.eta(),rho_area.phi(),rho_area.m()])

            if(rho_area.pt() >= pmu_cyl[0]):
                delete_indices.append(i)
                continue

            # Now, find the corrected jet 4-momentum.
            # Here, we're careful to use the input jet object,
            # *not* the raw jet (that we used to access the area).
            pmu = pmu_dict[i]
            jet = fj.PseudoJet(pmu[1],pmu[2],pmu[3],pmu[0]) # NOTE: won't have any clustering history, unlike raw_jet, but that's OK here

            corr_jet = jet - rho_area
            if(corr_jet.pt() < self.pt_min):
                delete_indices.append(i)
                continue
            self.buffers['Pmu'][i]     = np.array([corr_jet.e(),corr_jet.px(),corr_jet.py(),corr_jet.pz()])
            self.buffers['Pmu_cyl'][i] = np.array([corr_jet.pt(),corr_jet.eta(),corr_jet.phi(),corr_jet.m()])
            self.buffers['Passed'][i] = True
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        self._addToBuffer(obj)
        return

    def _set_rho_input(self,obj : 'JetFinder'):
        self.rho_name = '{}.Rho'.format(self.rho_input)
        self.rho_edge_name = '{}.Edges.Eta'.format(self.rho_input)

        for name in [self.rho_name, self.rho_edge_name]:
            if(name not in obj.input_collection_arrays.keys()):
                obj.input_buffer.read_branch(name)
        return

    def _fetch_rho(self,obj : 'JetFinder'):
        self.rho = np.array(obj.input_buffer[self.rho_name])
        self.rho_eta_edges = np.array(obj.input_buffer[self.rho_edge_name])

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

        for key,value in self.branch_names.items():
            if('Pmu' not in key):
                continue

            if('RhoArea.Pmu' in key and not self.save_rho_area):
                continue

            if(value not in obj.output_buffer.keys()):
                obj.output_buffer.create_array(value,ndim=2,dtype=np.dtype('f8'))
        obj.output_buffer.create_array(self.branch_names['Passed'],ndim=1,dtype=np.dtype('f8'))
        return

    def _createBranchNames(self,obj:'JetFinder'):
        self.branch_names = {key:'{}.{}'.format(self.obj_name_output,key) for key in self.branch_names}
        return

    def _addToBuffer(self,obj : 'JetFinder'):
        """
        Adds the corrected jet momentum to the buffer.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        for key,value in self.branch_names.items():
            if('RhoArea.Pmu' in key and not self.save_rho_area):
                continue
            obj.output_buffer.set(value,obj._i,self.buffers[key])
