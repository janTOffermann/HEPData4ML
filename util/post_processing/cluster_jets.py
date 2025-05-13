# The purpose of this code is to apply the Johns Hopkins top tagger (arXiv:0806.0848 [hep-ph])
# to the jets in the dataset.
import sys,itertools
import subprocess as sub
import numpy as np
import h5py as h5
# from util.post_processing.utils.jhtagger import JHTaggerSetup, JHTopTagger
from util.fastjet import FastJetSetup
from util.calcs import Calculator
import util.qol_utils.qol_util as qu

class JetFinder:
    """
    This class uses the Fastjet library to perform jet clustering, on some
    (arbitrary) set of inputs representing four-momenta of some objects.
    """

    def __init__(self, jet_algorithm='anti_kt',radius=0.4, jet_name='AK04Jets', save_constituents=True, fastjet_dir=None,verbose=False):

        self.status = False

        self.jet_algorithm_name = jet_algorithm
        self.jet_name = jet_name
        self.jet_radius = radius
        self.min_pt = 15.
        self.max_eta = 4.
        self.constituents_flag = save_constituents
        self.n_jets_max = 10 # max number of jets to save per event (will be pt-ordered)
        self.n_constituents_max = 200 # max number of constituents to save per jet

        self.fastjet_dir = fastjet_dir
        self.fastjet_init_flag = False

        # some fastjet-specific vars, for internal usage
        self.jet_algorithm = None
        self.jetdef = None
        self.cluster_sequence = None
        self.jets = None
        self.jet_vectors = None
        self.jet_vectors_cyl = None
        self.constituent_vectors = None
        self.constituent_vectors_cyl = None
        self.constituent_indices = None
        self.buffer = {}

        self.input_collection_arrays = None

        self.SetVerbosity(verbose)
        self.configurator = None
        self.calculator = None

        self.print_prefix = '\n\tJetFinder'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tRunning JetFinder:'
        self.progress_bar_suffix = 'Complete'
        self.setup = None
        self.tagger = None

        self.copts = 9

        self.error = False

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

    def _initialize_fastjet(self):

        if(self.fastjet_init_flag):
            return

        if(self.fastjet_dir is None):
            if(self.configurator is None):
                return
            else: # Fetch fastjet directory from configurator. This is foreseen as the "typical" usage.
                self.fastjet_dir = self.configurator.GetFastjetDirectory()

        # Now, if possible we initialize fastjet
        if(self.fastjet_dir is not None):
            self.SetupFastjet()
            sys.path.append(self.fastjet_dir)
            import fastjet as fj # This is where fastjet is really imported & cached by Python, but there may be other import statements peppered throughout since this has limited scope.
            self.fastjet_init_flag = True
        return

    def RequiresIndexing(self):
        return False

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetInputCollections(self,collections):
        self.input_collections = collections

    def SetMaxEta(self,eta):
        self.max_eta = np.abs(eta)

    def SetMinPt(self,pt):
        self.min_pt = pt # GeV

    def SetNConstituentsMax(self,n):
        self.n_constituents_max = n

    def SetRadius(self,radius):
        self.radius = radius

    def SetConfigurator(self,configurator):
        self.configurator = configurator


    def _input_consistency_check(self):
        f = h5.File(self.h5_file,'r')
        keys = list(f.keys())
        f.close()

        cleaned_collections = []
        for collection in self.input_collections:
            # we will be using the Cartesian versions of each collection for clustering
            branch = '{}.Pmu'.format(collection)
            found = branch in keys
            if(not found):
                self._print('Warning: Did not find key {} in file {}. Disabling as input...'.format(branch,self.h5_file))
            else:
                cleaned_collections.append(collection)
        self.input_collections = cleaned_collections
        if(len(self.input_collections)==0):
            self._print('Error: No input collections.')
            return False
        return True

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in the input HDF5 file, and places the required arrays in memory.
        """

        # If already initialized, no need to do it again.
        if(self.status):
            return

        # Input consistency check.
        if(not self._input_consistency_check()):
            self.status = False
            return


        # Initialize fastjet.
        self._initialize_fastjet()
        if(not self.fastjet_init_flag):
            self.status = False
            return

        # Initialize the fastjet jet definition
        self._initialize_jet_definition()

        # Read in the input 4-momenta from the input file.
        # TODO: Currently reading things into memory, which is not ideal for large files.
        #       Should ultimately move towards batching things, which will require
        #       keeping the input file open the whole time.
        f = h5.File(self.h5_file,'r')
        self.input_collection_arrays = {
            key:f[key][:] for key in self.input_collections
        }
        self.nevents = f[self.input_collections[0]].shape[0]

        # With nevents defined, we can initialize the buffer.
        self._initializeBuffer()

        self.status = True
        f.close()

    def Process(self):
        self.Initialize()

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        for i in range(self.nevents):

            # Gather the different input collections together, into one array of four-momenta.
            input_vecs = np.vstack([self.input_collection_arrays[key][i] for key in self.input_collection_arrays.keys()])
            self._clusterJets(input_vecs) # fills self.jets

            # optionally extract information on jet constituents
            if(self.constituents_flag):
                self._fetchJetConstituents() # fills self.constituent_vectors, self.constituent_indices

            # now write to buffer
            self._writeToBuffer()

            if(self.verbose):
                qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        return

    def Write(self,output_file):
        if(output_file is None): output_file = self.h5_file
        else: sub.check_call(['cp',self.h5_file,output_file])
        if(self.verbose):
            self._print('Writing output to {}.'.format(output_file))
        f = h5.File(output_file,'a')

        for key,val in self.buffer.items():
            d = f.create_dataset(key, data=val, compression='gzip',compression_opts=self.copts)
        f.close()
        return output_file

    # Using a generic signature -- should consider making the various post-processors inherit from a single parent class!
    def __call__(self,hepmc_file,h5_file,output_file=None,verbose=None, copts=9, key=None):
        if(verbose is not None): self.SetVerbosity(verbose)
        self.SetH5EventFile(h5_file)
        self.Process()
        output_file = self.Write(output_file)
        return output_file

    # --- Utility functions for the class ---
    def SetupFastJet(self):
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        return

    def _parse_jet_algorithm(self):
        self._initialize_fastjet()
        import fastjet as fj

        self.jet_algorithm = None
        for key in ['anti_kt','anti kt','anti-kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.antikt_algorithm
                return True

        for key in ['kt']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.kt_algorithm
                return True

        # TODO: C/A algorithm
        return False

    def _initialize_jet_definition(self):
        self._initialize_fastjet()
        import fastjet as fj
        # Print the FastJet banner -- it's unavoidable (other packages don't do this!).
        fj.ClusterSequence.print_banner()

        # determine what jet algorithm to use
        self._parse_jet_algorithm()

        self.jetdef = fj.JetDefinition(self.jet_algorithm, self.radius)
        return

    def _clusterJets(self,vecs):
        self._initialize_fastjet()
        import fastjet as fj # hacky, but will work because SetupFastJet() was run in __init__()

        # vecs has format (E,px,py,pz) -- FastJet uses (px,py,pz,E) so we must modify it. Using np.roll.
        pj = [fj.PseudoJet(*x) for x in np.roll(vecs,-1,axis=-1)]

        # Attach indices to the pseudojet objects, so that we can trace them through jet clustering.
        # Indices will correspond to the order they were input (with zero-indexing).
        for i,pseudojet in enumerate(pj):
            pseudojet.set_user_index(i)

        # selector = fj.SelectorPtMin(jet_config['jet_min_pt']) & fj.SelectorAbsEtaMax(jet_config['jet_max_eta'])
        # Note: Switched from the old method, this is more verbose but seems to do the same thing anyway.
        self.cluster_sequence = fj.ClusterSequence(pj, self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
        self.jets = self.cluster_sequence.inclusive_jets()

        # Now we will apply our (optional) pt and eta cuts to the jets.
        # TODO: Should be achievable with Fastjet selector classes too.
        jet_eta = np.array([jet.eta() for jet in self.jets])
        jet_pt  = np.array([jet.pt()  for jet in self.jets])

        selected_eta = np.where(np.abs(jet_eta) <= np.abs(self.max_eta))[0]
        selected_pt  = np.where(jet_pt          >= self.min_pt )[0]
        selected = np.sort(np.intersect1d(selected_eta, selected_pt)) # TODO: Is the sort needed?
        self.jets = [self.jets[x] for x in selected]

        # sort by pT, truncate if necessary
        jet_pt  = np.array([jet.pt()  for jet in self.jets])
        pt_sorting = np.argsort(-jet_pt)
        if(len(pt_sorting) > self.n_jets_max):
            pt_sorting = pt_sorting[:self.n_jets_max]
        self.jets = list(itertools.compress(self.jets,pt_sorting))

        # also record the jet four-momenta explicitly
        self.jet_vectors = np.array([[x.e(), x.px(), x.py(),x.pz()] for x in self.jets])
        self.jet_vectors_cyl = np.array([[x.pt(), x.eta(), x.phi(),x.m()] for x in self.jets])

    def _fetchJetConstituents(self):
        results = [self._fetchJetConstituentsSingle(jet, self.n_constituents_max) for jet in self.jets]
        self.constituent_vectors = [x[0] for x in results]
        self.constituent_vectors_cyl = [x[1] for x in results]
        self.constituent_indices = [x[2] for x in results]

    def _fetchJetConstituentsSingle(self,jet,n_constituents):
        pt,eta,phi,m,e,px,py,pz = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m(), x.e(), x.px(),x.py(),x.pz()] for x in jet.constituents()]),8)
        pt  = pt.flatten() # use this for sorting
        eta = eta.flatten()
        phi = phi.flatten()
        m   = m.flatten()
        e   = e.flatten()
        px  = px.flatten()
        py  = py.flatten()
        pz  = pz.flatten()

        # The indices of the jet constituents, corresponding with the order in which they
        # were passed to jet clustering.
        indices = np.array([x.user_index() for x in jet.constituents()],dtype=np.dtype('i4')).flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt)
        l = int(np.minimum(n_constituents,len(pt)))

        pt  =  pt[sorting][:l]
        eta = eta[sorting][:l]
        phi = phi[sorting][:l]
        m   = m[sorting][:l]

        e  =  e[sorting][:l]
        px = px[sorting][:l]
        py = py[sorting][:l]
        pz = pz[sorting][:l]
        indices = indices[sorting][:l]

        vecs = np.array([e,px,py,pz],dtype=np.dtype('f8')).T
        vecs_cyl = np.array([pt,eta,phi,m],dtype=np.dtype('f8')).T

        return vecs, vecs_cyl, indices

    def _initializeBuffer(self):
        self.buffer['{}.N'.format(self.jet_name)] = np.zeros((self.nevents),dtype=np.dtype('i4'))
        self.buffer['{}.Pmu'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,4),dtype=np.dtype('f8'))
        self.buffer['{}.Pmu_cyl'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,4),dtype=np.dtype('f8'))
        if(self.constituents_flag):
            self.buffer['{}.Constituents.N'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max),dtype=np.dtype('i4'))
            self.buffer['{}.Constituents.Pmu'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,self.n_constituents_max, 4),dtype=np.dtype('f8'))
            self.buffer['{}.Constituents.Pmu_cyl'.format(self.jet_name)] = np.zeros((self.nevents,self.n_jets_max,self.n_constituents_max,4),dtype=np.dtype('f8'))
        return


    def _writeToBuffer(self,event_index):

        # fill jets
        self.buffer['{}.N'.format(self.jet_name)][event_index] = len(self.jet_vectors)
        self.buffer['{}.Pmu'.format(self.jet_name)][event_index] = self.jet_vectors
        self.buffer['{}.Pmu_cyl'.format(self.jet_name)][event_index] = self.jet_vectors_cyl

        # fill jet constituents
        if(self.constituents_flag):
            self.buffer['{}.Constituents.N'.format(self.jet_name)][event_index] = len(self.constituent_vectors)
            self.buffer['{}.Constituents.Pmu'.format(self.jet_name)][event_index] = self.constituent_vectors
            self.buffer['{}.Constituents.Pmu_cyl'.format(self.jet_name)][event_index] = self.constituent_vectors_cyl
        return

