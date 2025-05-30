import os,glob,pathlib
import subprocess as sub
import numpy as np
import ROOT as rt
from util.calcs import embed_array_inplace
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder

class JHTaggerSetup:
    """
    This class is for setting up the Johns Hopkins top tagger (JHTagger) code.
    This is some custom C++/ROOT code, that leverages the JHTagger implementation
    in the Fastjet library.
    """
    def __init__(self,configurator,jhtagger_dir=None):
        self.SetDirectory(jhtagger_dir)
        self.SetConfigurator(configurator)
        self.executable = '/bin/bash'
        self.lib_extensions = ['so','dylib']
        self.cmake_template = self.dir + '/../cmake_templates/CMakeLists_jhtagger.txt'
        self.status = False

    def SetDirectory(self,dir):
        if(dir is None):
            dir = os.path.dirname(os.path.abspath(__file__)) + '/../../root/jhtagger'
        self.dir = dir

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def FullPreparation(self):
        self.status = False
        try: self.LoadJHTagger()
        except:
            self.PrepareCMake() # writes a CMakeLists.txt file with current ROOT version, based on template.
            self.BuildJHTagger()
            self.LoadJHTagger()
        self.status = True
        return

    def PrepareCMake(self):
        with open(self.cmake_template,'r') as f:
            lines = f.readlines()

        # Get the ROOT version
        root_version = rt.__version__.split('/')[0]

        for i,line in enumerate(lines):
            lines[i] = lines[i].replace('ROOT_VERSION',root_version)

        cmake_file = self.dir + '/CMakeLists.txt'
        with open(cmake_file,'w') as f:
            for line in lines:
                f.write(line)
        return

    def BuildJHTagger(self):
        """
        Builds the JHTagger C++/ROOT library.
        This runs a shell script, where certain
        environment variables have been set by
        this Python script (via our configurator,
        which tells us where Fastjet is installed).
        """
        self.build_script = 'build.sh'
        fastjet_dir = self.configurator.GetFastjetDirectory()
        # TODO: This might crash if fastjet is not built yet? Could be an issue for a first run, need to test!
        fastjet_lib = glob.glob('{}/**/lib'.format(fastjet_dir),recursive=True)[0]
        fastjet_inc = glob.glob('{}/**/include'.format(fastjet_dir),recursive=True)[0]
        env = os.environ.copy()
        env['CMAKE_PREFIX_PATH'] += ':{}'.format(fastjet_lib)
        env['FASTJET_INCLUDE_DIR'] = fastjet_inc

        command = ['.',self.build_script]
        sub.check_call(command,shell=False,cwd=self.dir,env=env,executable=self.executable,stderr=sub.DEVNULL,stdout=sub.DEVNULL)
        return

    def LoadJHTagger(self,quiet=False):
        # Load our custom ROOT library.
        try:
            a = rt.JHTagger
            return
        except:
            pass

        # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
        # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
        build_dir = self.dir + '/jhtagger/build'
        custom_lib_paths = [os.path.realpath('{}/lib/libJHTagger.{}'.format(build_dir,x)) for x in self.lib_extensions]
        custom_inc_paths = os.path.realpath('{}/include/jhtagger'.format(build_dir))
        custom_inc_paths = glob.glob(custom_inc_paths + '/**/*.h',recursive=True)

        # Check for any of the libraries.
        found_libary = False
        for libpath in custom_lib_paths:
            if(pathlib.Path(libpath).exists()):
                custom_lib_path = libpath
                found_libary = True
                break

        if(not found_libary):
            if(not quiet): print('Error: The JHTagger lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The JHTagger lib did not load properly.'
        return

class JohnsHopkinsTagger:
    """
    The Johns Hopkins top tagger. This class leverages our custom
    ROOT::JHTagger::JohnnyTagger() class, that interfaces with the Fastjet
    JH tagger implementation.

    See: https://arxiv.org/abs/0806.0848 [Phys.Rev.Lett. 101 (2008) 142001]

    """
    def __init__(self,delta_p=0.1,delta_r=0.19,cos_theta_W_max=0.7,top_mass_range=(150.,200.),W_mass_range=(65.,95.), mode='filter',tag_name=None, n_w_constituents_max=100):
        self.mode = mode
        assert self.mode in ['tag','filter']
        self.tag_name = tag_name
        self.w_name = None
        self.w_nconst_name = None
        self.w_constituents_name = None
        self.w_constituents_name = None
        # Prepare the JHTagger library.
        self.jh_setup = None

        self.tagger = None
        self.SetDeltaP(delta_p,init=False)
        self.SetDeltaR(delta_r,init=False)
        self.SetCosThetaWMax(cos_theta_W_max,init=False)
        self.SetTopMassRange(*top_mass_range,init=False)
        self.SetWMassRange(*W_mass_range,init=False)
        # self.InitTagger()

        # Transient, per-jet variables
        self.tag_status = False
        self.w_candidate = None
        self.w_candidate_cyl = None
        self.w_constituent = None
        self.w_constituent_cyl = None

        # Event-level; one extra dim w.r.t. above
        # TODO: Rework this? Might be confusing naming scheme
        self.tags = None
        self.w_candidates = None
        self.w_candidates_cyl = None
        self.n_constituents = None
        self.w_constituents = None
        self.w_constituents_cyl = None

        self.n_w_constituents_max = n_w_constituents_max

        self.print_prefix = '\n\t\tJohnsHopkinsTagger'

    def SetDeltaP(self,p,init=True):
        self.delta_p = p
        if(init): self.InitTagger()

    def SetDeltaR(self,r,init=True):
        self.delta_r = r
        if(init): self.InitTagger()

    def SetCosThetaWMax(self,c,init=True):
        self.cos_theta_W_max = c
        if(init): self.InitTagger()

    def SetTopMassRange(self,mmin,mmax,init=True):
        self.top_mass_range = (mmin,mmax)
        if(init): self.InitTagger()

    def SetWMassRange(self,mmin,mmax,init=True):
        self.W_mass_range = (mmin,mmax)
        if(init): self.InitTagger()

    def InitTagger(self):
        if(self.tagger is not None):
            del self.tagger # TODO: Is this necessary?
        self.tagger = rt.JHTagger.JohnnyTagger(self.delta_p,self.delta_r,self.cos_theta_W_max,*self.top_mass_range,*self.W_mass_range)

    def _tag(self,vecs):
        #NOTE: This usage of TagJet() is a bit awkward, but I've had some issues when trying to pass a fastjet.PseudoJet object.
        #      I think this has to do with some weirdness around the Fastjet Python interface, or Python/C++ interfaces in general.
        self.tagger.TagJet( vecs[:,0].flatten(), vecs[:,1].flatten(), vecs[:,2].flatten(), vecs[:,3].flatten() )
        self.tag_status = self.tagger.GetStatus()
        self.w_candidate = None
        if(self.tag_status):
            try:
                w = self.tagger.GetWCandidate() # Fastjet::PseudoJet (C++ type, PyROOT seems to handle interface here!)
                self.w_candidate     = np.array([w.e(),w.px(),w.py(),w.pz()  ]) # TODO: Can this array return be handled by the C++/ROOT class?
                self.w_candidate_cyl = np.array([w.pt(),w.eta(),w.phi(),w.m()])
            except:
                self.w_candidate     = np.full(4,np.nan)
                self.w_candidate_cyl = np.full(4,np.nan)
        return

    def _getWConstituents(self):
        if(not self.tag_status): return None
        pt = np.array(self.tagger.GetWCandidateConstituentsProperty("pt"))
        n = pt.shape[0]
        if(n == 0):
            return
        ordering = np.argsort(-pt)
        self.w_constituent     = np.vstack([self.tagger.GetWCandidateConstituentsProperty(x) for x in ["E","px","py","pz"]]  ).T[ordering]
        self.w_constituent_cyl = np.vstack([self.tagger.GetWCandidateConstituentsProperty(x) for x in ["pt","eta","phi","m"]]).T[ordering]

    def ModifyInitialization(self, obj : 'JetFinder'):
        #NOTE: Might want to move this to constructor, which will need to take obj as input.
        self.jh_setup = JHTaggerSetup(obj.configurator)
        # self.jh_setup.SetConfigurator(obj.configurator)
        self.jh_setup.FullPreparation()
        self.InitTagger()
        return

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will tag jets with the JH tagger, and fill the corresponding branches.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        # Fetch the jet constituents, just to be safe -- this makes sure that they are up-to-date.
        obj._fetchJetConstituents()

        self.tags = {i:False for i in obj.jets_dict.keys()}
        self.w_candidates = {i:np.zeros(4) for i in obj.jets_dict.keys()}
        self.w_candidates_cyl = {i:np.zeros(4) for i in obj.jets_dict.keys()}
        self.n_constituents = {i:0 for i in obj.jets_dict.keys()}
        self.w_constituents = {i:np.zeros((self.n_w_constituents_max,4)) for i in obj.jets_dict.keys()}
        self.w_constituents_cyl = {i:np.zeros((self.n_w_constituents_max,4)) for i in obj.jets_dict.keys()}

        for key in obj.jets_dict.keys():
            self._tag(obj.constituent_vectors[key]) # fills self.tag_status, self.w_candidate
            self.tags[key] = self.tag_status
            if(self.tag_status): # if not tagged, there is no W -- so we can safely skip filling
                self._getWConstituents() # fills self.w_constituent, self.w_constituent_cyl
                self.n_constituents[key] = len(self.w_constituent)
                embed_array_inplace(self.w_constituent,self.w_constituents[key])
                embed_array_inplace(self.w_constituent_cyl,self.w_constituents_cyl[key])

        # for i in range(len(obj.constituent_vectors)):
        #     self._tag(obj.constituent_vectors[i]) # fills self.tag_status, self.w_candidate
        #     self.tags[i] = self.tag_status
        #     self.w_candidates[i] = self.w_candidate
        #     if(self.tag_status): # if not tagged, there is no W -- can safely skip
        #         self._getWConstituents() # fills self.w_constituent, self.w_constituent_cyl
        #         self.n_constituents[i] = len(self.w_constituent)
        #         embed_array_inplace(self.w_constituent,self.w_constituents[i])
        #         embed_array_inplace(self.w_constituent_cyl,self.w_constituents_cyl[i])

        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            obj._updateJetDictionary()
            # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
            obj._jetsToVectors()
            obj._fetchJetConstituents()

    def ModifyConstituents(self, obj : 'JetFinder'):
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addFlagToBuffer(obj)
            self._addWToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Used if self.mode=='tag', in which case we're writing all jets,
        and including a new branch to indicate whether or not a jet is
        JH-tagged.
        """
        self._createBranchNames(obj)

        if(self.tag_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.tag_name,(obj.n_jets_max,),dtype=bool)

        if(self.w_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.w_name,         (obj.n_jets_max,4),dtype=np.dtype('f8'))
            obj.buffer.create_array(self.w_name + '_cyl',(obj.n_jets_max,4),dtype=np.dtype('f8'))

        if(self.w_constituents_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.w_nconst_name, (obj.n_jets_max),dtype=np.dtype('i4'))
            obj.buffer.create_array(self.w_constituents_name,          (obj.n_jets_max,self.n_w_constituents_max,4),dtype=np.dtype('f8'))
            obj.buffer.create_array(self.w_constituents_name + '_cyl', (obj.n_jets_max,self.n_w_constituents_max,4),dtype=np.dtype('f8'))
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        if(self.tag_name is None):
            self.tag_name = '{}.JHTag'.format(obj.jet_name)

        # Also create keys corresponding to the W-boson candidate, and its constituents,
        # that are identified by the JH tagger.
        self.w_name = '{}.WBosonCandidate.Pmu'.format(self.tag_name)
        self.w_nconst_name = '{}.WBosonCandidate.Constituents.N'.format(self.tag_name)
        self.w_constituents_name = '{}.WBosonCandidate.Constituents.Pmu'.format(self.tag_name)

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the JH tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.tag_name,obj._i,[self.tags[i] for i in obj.jet_ordering])

    def _addWToBuffer(self,obj : 'JetFinder'):
        """
        Adds the JH W candidate info to the buffer.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.w_name,obj._i,[self.w_candidates[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_name + '_cyl',obj._i,[self.w_candidates_cyl[i] for i in obj.jet_ordering])

        obj.buffer.set(self.w_nconst_name,obj._i,[self.n_constituents[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_constituents_name,obj._i,[self.w_constituents[i] for i in obj.jet_ordering])
        obj.buffer.set(self.w_constituents_name + '_cyl',obj._i,[self.w_constituents_cyl[i] for i in obj.jet_ordering])

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
