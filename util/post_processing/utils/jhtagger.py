import os,glob,pathlib
import subprocess as sub
import numpy as np
import ROOT as rt

class JHTaggerSetup:
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

    def LoadJHTagger(self):
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
            print('Error: The JHTagger lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The JHTagger lib did not load properly.'
        return

class JHTopTagger:
    def __init__(self,delta_p=0.1,delta_r=0.19,cos_theta_W_max=0.7,top_mass_range=(150.,200.),W_mass_range=(65.,95.)):
        self.SetDeltaP(delta_p,init=False)
        self.SetDeltaR(delta_r,init=False)
        self.SetCosThetaWMax(cos_theta_W_max,init=False)
        self.SetTopMassRange(*top_mass_range,init=False)
        self.SetWMassRange(*W_mass_range,init=False)
        self.InitTagger()

        self.status = False
        self.w_candidate = None

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
        self.tagger = rt.JHTagger.JohnnyTagger(self.delta_p,self.delta_r,self.cos_theta_W_max,*self.top_mass_range,*self.W_mass_range)

    def TagEvent(self,pmu):
        # Have to do things like this since fastjet objects on the Python side of things aren't actually identical
        # to fastjet objects on the C++ side of things. I.e. even if our custom ROOT library's function takes
        # a fastjet::PseudoJet as an argument, we can't pass it a fastjet.PseudoJet from Python -- it won't
        # be recognized as the same type.
        self.tagger.TagJet( pmu[:,0].flatten(), pmu[:,1].flatten(), pmu[:,2].flatten(), pmu[:,3].flatten() )
        self.status = self.tagger.GetStatus()
        self.w_candidate = None
        if(self.status):
            w = self.tagger.GetWCandidate() # Fastjet::PseudoJet (C++ type, PyROOT seems to handle interface here!)
            self.w_candidate = np.array([w.e(),w.px(),w.py(),w.pz()]) # TODO: Can this array return be handled by the C++/ROOT class?
        return

    def GetStatus(self):
        return self.status

    def GetWCandidate(self):
        return self.w_candidate

    def GetWCandidateConstituents(self):
        if(not self.status): return None
        # n = self.tagger.GetWCandidateNConstituents() #TODO: Why does n seem to be 1 too small in some tests?
        # if(n == 0):
        #     return None
        pt = np.array(self.tagger.GetWCandidateConstituentsProperty("pt")) # used for ordering
        n = pt.shape[0]
        if(n == 0):
            return None
        ordering = np.argsort(-pt)
        vecs = np.zeros((n,4))
        vecs[:,0] = np.array(self.tagger.GetWCandidateConstituentsProperty("E"))[ordering]
        vecs[:,1] = np.array(self.tagger.GetWCandidateConstituentsProperty("px"))[ordering]
        vecs[:,2] = np.array(self.tagger.GetWCandidateConstituentsProperty("py"))[ordering]
        vecs[:,3] = np.array(self.tagger.GetWCandidateConstituentsProperty("pz"))[ordering]
        return vecs
