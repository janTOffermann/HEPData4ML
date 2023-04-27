# For loading our custom VectorCalcs library.
import ROOT as rt
import os,glob,pathlib
import subprocess as sub
import util.qol_utils.qol_util as qu

class VectorCalcsManager:
    def __init__(self,verbose=True):
        self.executable = '/bin/bash'
        self.verbose = verbose
        self.lib_extensions = ['so','dylib']
        self.build_script = 'build.sh'
        self.script_dir = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root/vcalcs')
        self.cmake_template = self.script_dir + '/../cmake_templates/CMakeLists_vectorcalcs.txt'
        self.status = False

    def FullPreparation(self,force=False):
        self.status = False
        try:
            self.PrepareCMake() # writes a CMakeLists.txt file with current ROOT version, based on template.
            self.BuildVectorCalcs(force) # builds VectorCalcs library, if it is not built yet. (force=True will make it always build)
            self.LoadVectorCalcs()
            self.status = True
        except:
            pass
        return

    def GetStatus(self):
        return self.status

    def PrepareCMake(self):
        with open(self.cmake_template,'r') as f:
            lines = f.readlines()

        # Get the ROOT version
        root_version = rt.__version__.split('/')[0]

        for i,line in enumerate(lines):
            lines[i] = lines[i].replace('ROOT_VERSION',root_version)

        cmake_file = self.script_dir + '/CMakeLists.txt'
        with open(cmake_file,'w') as f:
            for line in lines:
                f.write(line)
        return

    def BuildVectorCalcs(self,force=False):
        # Check if the VectorCalcs library is already built, by trying to load it.
        # We will disable printouts for this check, since it's okay if it fails.
        if(not force):
            with qu.stdout_redirected():
                try:
                    self.LoadVectorCalcs()
                    return
                except: pass

        comm = ['.',self.build_script]
        env = os.environ.copy()
        if(self.verbose): print('Building VectorCalcs library.')
        sub.check_call(comm,shell=False,cwd=self.script_dir,env=os.environ.copy(),executable=self.executable,stderr=sub.DEVNULL,stdout=sub.DEVNULL)
        return

    def LoadVectorCalcs(self):
        # Load our custom ROOT library.
        try:
            a = rt.VectorCalcs
            return
        except: pass

        # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
        # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
        # TODO: Is there a more elegant workaround for this?
        custom_lib_paths = [os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root/vcalcs/vectorcalcs/build/lib/libVectorCalcs.{}'.format(x)) for x in self.lib_extensions]
        custom_inc_paths = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root/vcalcs/vectorcalcs/build/include/vectorcalcs')
        custom_inc_paths = glob.glob(custom_inc_paths + '/**/*.h',recursive=True)

        # Check for any of the libraries.
        found_libary = False
        for libpath in custom_lib_paths:
            if(pathlib.Path(libpath).exists()):
                custom_lib_path = libpath
                found_libary = True
                break

        if(not found_libary):
            print('Error: The VectorCalcs lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The VectorCalcs lib did not load properly.'
        return