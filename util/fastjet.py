import sys, os, glob,pathlib, operator
import numpy as np
import subprocess as sub
from util.qol_utils.misc import stdout_redirected

class FastJetSetup:
    def __init__(self,fastjet_dir=None, full_setup=False,verbose=True):
        self.fastjet_version = '3.4.2'
        self.SetDirectory(fastjet_dir)
        if(full_setup): self.PrepFastjet(verbose=verbose)

    def SetDirectory(self,fastjet_dir=None):
        self.fastjet_dir = fastjet_dir
        if(self.fastjet_dir is None):
                self.fastjet_dir = os.path.dirname(os.path.abspath(__file__)) + '/../fastjet'
        self.fastjet_dir = os.path.normpath(self.fastjet_dir)
        os.makedirs(self.fastjet_dir,exist_ok=True)

        self.SetPythonDirectory() # attempts to get the Python library subdir, only sets it if found

        self.logfile = '{}/log.stdout'.format(self.fastjet_dir)
        self.errfile = '{}/log.stderr'.format(self.fastjet_dir)

        self.source_dir = '{}/fastjet-{}'.format(self.fastjet_dir,self.fastjet_version)
        self.install_dir = '{}/fastjet-install'.format(self.fastjet_dir)
        return

    def GetDirectory(self):
        return self.fastjet_dir

    def SetPythonDirectory(self):
        try:
            self.python_dir = os.path.normpath(glob.glob('{}/**/site-packages'.format(self.fastjet_dir),recursive=True)[0])
        except:
            self.python_dir = None

    def GetPythonDirectory(self):
        return self.python_dir

    def PrepFastjet(self, j=4, force=False, verbose=True):
        # Check if Fastjet is already built at destination.
        # Specifically, we will look for some Python-related files.
        if(not force):

            files_to_find = [
                '{}/**/site-packages/fastjet.py'.format(self.fastjet_dir),
                '{}/**/site-packages/_fastjet.a'.format(self.fastjet_dir),
                '{}/**/site-packages/_fastjet.so.0'.format(self.fastjet_dir)
            ]

            files_to_find = [glob.glob(x,recursive=True) for x in files_to_find]
            files_to_find = [len(x) for x in files_to_find]
            found_fastjet = False
            for entry in files_to_find:
                if(entry > 0):
                    found_fastjet = True # Loose condition since not all these files may exist (e.g. no .so on macos). Should be okay.
                    break
            if(found_fastjet):
                if(verbose): print('Found existing Fastjet installation with Python extension @ {}.'.format(self.fastjet_dir))
                # TODO: Do something else here?
                return

        # Fastjet was not found -> download and build.
        self.DownloadFastjet(verbose=verbose)
        self.BuildFastjet(verbose=verbose)

    def DownloadFastjet(self,force=False,verbose=True):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            if(not force and pathlib.Path(self.source_dir).exists()): # check if fastjet source is already present, if so we do not need to download it again
                print('Found local fastjet source code, but fastjet is not locally built yet.')
                return

            # Fetch the Fastjet source
            fastjet_download = 'http://fastjet.fr/repo/fastjet-{}.tar.gz'.format(self.fastjet_version)
            fastjet_file = fastjet_download.split('/')[-1]
            if(verbose): print('Downloading fastjet from {}.'.format(fastjet_download))

            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget): sub.check_call(['wget', fastjet_download], shell=False, cwd=self.fastjet_dir, stdout=f, stderr=g)
            else: sub.check_call(['curl',fastjet_download,'-o',fastjet_file], shell=False, cwd=self.fastjet_dir, stdout=f, stderr=g)
            sub.check_call(['tar', 'zxvf', fastjet_file], shell=False, cwd=self.fastjet_dir, stdout=f, stderr=g)
            sub.check_call(['rm', fastjet_file], shell=False, cwd=self.fastjet_dir, stdout=f, stderr=g)

    def BuildFastjet(self,j=4,verbose=True):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            # Now configure. We create the python bindings.
            if(verbose): print('Configuring fastjet.')
            sub.check_call(['./configure', '--prefix={}'.format(self.install_dir), '--enable-pyext'],
                        shell=False, cwd=self.source_dir, stdout=f, stderr=g)

            # Now make and install. Will skip "make check".
            if(verbose): print('Making fastjet.')
            sub.check_call(['make', '-j{}'.format(j)],
                        shell=False, cwd=self.source_dir, stdout=f, stderr=g)
            if(verbose): print('Installing fastjet.')
            sub.check_call(['make', 'install'],
                        shell=False, cwd=self.source_dir, stdout=f, stderr=g)
        self.SetPythonDirectory()
        return


class JetFinderBase:
    """
    This is a base class, for using the Fastjet library to perform jet clustering.
    It is leveraged by the JetFinder class in utils/post_processing/jets.py .
    This base class exists to avoid some possible import loops.
    It lacks some useful functionality that JetFinder includes.
    """

    def __init__(self,fastjet_dir=None):

        self.jet_algorithm_name = ''
        self.jet_name = ''
        self.radius = 0.4
        self.min_pt = 0.
        self.max_eta = 10.

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
        self.user_info = None
        self.pt_sorting = None # allows access to the sorting array

        self.input_vecs = None
        self.configurator = None

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
            self._setupFastJet()
            sys.path.append(self.fastjet_dir)
            import fastjet as fj # This is where fastjet is really imported & cached by Python, but there may be other import statements peppered throughout since this has limited scope.
            self.fastjet_init_flag = True
        return

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

    def SetInputs(self,vecs):
        self.input_vecs = vecs

    def _setupFastJet(self):
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

        for key in ['c/a','cambridge','aachen']:
            if(key in self.jet_algorithm_name.lower()):
                self.jet_algorithm = fj.cambridge_aachen_algorithm
                return True

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

    def _clusterJets(self):
        self._initialize_fastjet()
        import fastjet as fj # hacky, but will work because _setupFastJet() was run in __init__()

        # vecs has format (E,px,py,pz) -- FastJet uses (px,py,pz,E) so we must modify it. Using np.roll.
        pj = [fj.PseudoJet(*x) for x in np.roll(self.input_vecs,-1,axis=-1)]

        # Attach indices to the pseudojet objects, so that we can trace them through jet clustering.
        # Indices will correspond to the order they were input (with zero-indexing).
        for i,pseudojet in enumerate(pj):
            pseudojet.set_user_index(i)

        # Attach any optional information to the pseudojet objects. This can be leveraged by other classes
        # or extensions.
        if(self.user_info is not None):
            for idx, val in self.user_info.items(): # in practice, val will be a dictionary itself -- allows for attaching multiple things
                pj[idx].set_python_info(val)


        # selector = fj.SelectorPtMin(jet_config['jet_min_pt']) & fj.SelectorAbsEtaMax(jet_config['jet_max_eta'])
        # Note: Switched from the old method, this is more verbose but seems to do the same thing anyway.
        self.cluster_sequence = fj.ClusterSequence(pj, self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
        self.jets = self.cluster_sequence.inclusive_jets()

        # Now we will apply our (optional) pt and eta cuts to the jets.
        # TODO: Make this a post-processor of the jets.
        jet_eta = np.array([jet.eta() for jet in self.jets])
        jet_pt  = np.array([jet.pt()  for jet in self.jets])

        selected_eta = np.where(np.abs(jet_eta) <= np.abs(self.max_eta))[0]
        selected_pt  = np.where(jet_pt          >= self.min_pt )[0]
        selected = np.sort(np.intersect1d(selected_eta, selected_pt)) # TODO: Is the sort needed?
        self.jets = [self.jets[x] for x in selected]
        self._jetsToVectors()

    def _jetsToVectors(self):
        self.jet_vectors = np.array([[x.e(), x.px(), x.py(),x.pz()] for x in self.jets])
        self.jet_vectors_cyl = np.array([[x.pt(), x.eta(), x.phi(),x.m()] for x in self.jets])

    def _ptSort(self):
        """
        Sorts jets by decreasing pT, and truncates to
        take only the first self.n_jets_max jets.
        """
        if len(self.jets) <= 1:
            return

        jet_pt  = np.array([jet.pt() for jet in self.jets])
        self.pt_sorting = np.argsort(-jet_pt)
        if(len(self.pt_sorting) > self.n_jets_max):
            self.pt_sorting = self.pt_sorting[:self.n_jets_max]

        self.jets = list(operator.itemgetter(*self.pt_sorting)(self.jets)) #NOTE: Possibly a little obscure, but maybe more efficient than list comprehension? -Jan
        self._jetsToVectors()

    def _fetchJetConstituents(self):
        results = [self._fetchJetConstituentsSingle(jet, self.n_constituents_max) for jet in self.jets]
        self.constituent_vectors = [x[0] for x in results]
        self.constituent_vectors_cyl = [x[1] for x in results]
        self.constituent_indices = [x[2] for x in results]

    def _fetchJetConstituentsSingle(self,jet,n_constituents=-1):
        pt,eta,phi,m,e,px,py,pz = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m(), x.e(), x.px(),x.py(),x.pz()] for x in jet.constituents()]),8)

        # The indices of the jet constituents, corresponding with the order in which they
        # were passed to jet clustering.
        indices = np.array([x.user_index() for x in jet.constituents()],dtype=np.dtype('i4')).flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt.flatten())
        l = len(pt)
        if(n_constituents > 0):
            l = int(np.minimum(n_constituents,l))
        vecs = np.vstack([x.flatten() for x in [e,px,py,pz]]).T[sorting][:l]
        vecs_cyl = np.vstack([x.flatten() for x in [pt,eta,phi,m]]).T[sorting][:l]
        indices = indices[sorting][:l]

        return vecs, vecs_cyl, indices


class ParticleInfo(object):
    """Illustrative class for use in assigning pythonic user information
    to a PseudoJet.
    """
    def __init__(self, particle_index, status, pdg_id=0):
        self.particle_index = particle_index
        self.status = status
        self.pdg_id = pdg_id

    def __str__(self):
        return "particle_index={0}, status={1}, pdg_id={2}".format(
            self.particle_index, self.status, self.pdg_id)