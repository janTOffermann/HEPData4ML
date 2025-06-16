
from util.qol_utils.misc import stdout_redirected
from util.qol_utils.progress_bar import printProgressBar
import subprocess as sub
import sys,os,glob,re
from typing import Optional,Union

class HepMCSetup:
    """
    The purpose of this class is to install HepMC3 if necessary.
    This is achieved via the PrepHepMC() function, which will
    attempt to access the HepMC3 Python interface -- and optionally
    try to access the HepMC3/ROOT Python interface.
    If it is unable to do so, it will download and install HepMC3
    into the external/hepmc subdirectory.
    """

    def __init__(self,verbose:bool=False):
        self.version = '3.3.1' # what version we'll try to install, if necessary
        self.SetVerbose(verbose)
        self.prefix = self.__class__.__name__
        self.python_dir = None

    def SetVerbose(self,verbose:bool):
        self.verbose = verbose

    def SetDirectory(self,hepmc_dir:Optional[str]=None):
        self.hepmc_dir = hepmc_dir
        if(self.hepmc_dir is None):
                self.hepmc_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/hepmc'
        self.hepmc_dir = os.path.normpath(self.hepmc_dir)
        os.makedirs(self.hepmc_dir,exist_ok=True)

        # self.SetPythonDirectory() # attempts to get the Python library subdir, only sets it if found

        self.logfile = '{}/log.stdout'.format(self.hepmc_dir)
        self.errfile = '{}/log.stderr'.format(self.hepmc_dir)

        self.source_dir = '{}/HepMC3-{}'.format(self.hepmc_dir,self.version)
        self.build_dir = '{}/hepmc3-build'.format(self.hepmc_dir)
        self.install_dir = '{}/hepmc3-install'.format(self.hepmc_dir)
        return

    def SetPythonDirectory(self):
        try:
            self.python_dir = os.path.normpath(glob.glob('{}/**/site-packages'.format(self.install_dir),recursive=True)[0])
        except:
            self.python_dir = None

    def GetPythonDirectory(self) -> Union[str,None]:
        """
        Returns the 'site-packages' directory, where the pyHepMC library is located.
        Adding this to the $PYTHONPATH (sys.path) will allow one to import it.
        """
        return self.python_dir

    def DownloadHepMC3(self):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:

            # Fetch the HepMC3 source
            hepmc_download = 'http://hepmc.web.cern.ch/hepmc/releases/HepMC3-{}.tar.gz'.format(self.version)
            hepmc_file = hepmc_download.split('/')[-1]
            if(self.verbose): self._print('Downloading HepMC3 from {}.'.format(hepmc_download))

            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget): sub.check_call(['wget', hepmc_download], cwd=self.hepmc_dir, stdout=f, stderr=g)
            else: sub.check_call(['curl',hepmc_download,'-o',hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)
            sub.check_call(['tar', 'xzf', hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)
            sub.check_call(['rm', hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)

    def MakeHepMC3(self,j:int=4):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:

            os.makedirs(self.build_dir,exist_ok=True)

            if(self.verbose): self._print('Configuring HepMC3 with CMake.')
            # CMake configuration
            command = [
                'cmake',
                '-DCMAKE_INSTALL_PREFIX={}'.format(self.install_dir),
                '-DHEPMC3_ENABLE_ROOTIO:BOOL=ON',
                '-DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF',
                '-DHEPMC3_ENABLE_TEST:BOOL=OFF',
                '-DHEPMC3_INSTALL_INTERFACES:BOOL=ON',
                '-DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF',
                '-DHEPMC3_BUILD_DOCS:BOOL=OFF',
                '-DHEPMC3_ENABLE_PYTHON:BOOL=ON',
                '-DHEPMC3_PYTHON_VERSIONS={}.{}'.format(sys.version_info.major,sys.version_info.minor),
                '-DHEPMC3_Python_SITEARCH{a}{b}={c}/lib/python{a}.{b}/site-packages'.format(a=sys.version_info.major,b=sys.version_info.minor,c=self.install_dir),
                self.source_dir
            ]
            for entry in command:
                if(entry[0] == '-'):
                    self._print('\t{}'.format(entry))
            sub.check_call(command,cwd=self.build_dir, stdout=f, stderr=g)

            # Now make
            if(self.verbose): self._print('Building HepMC3.')
            command = ['make', '-j{}'.format(j)]
            self.run_make_with_progress(command)
            # sub.check_call(command,cwd=self.build_dir, stdout=f, stderr=g)

            if(self.verbose): self._print('Installing HepMC3.')
            command = ['make','install']
            sub.check_call(command,cwd=self.build_dir, stdout=f, stderr=g)

    def PrepHepMC(self, j:int=4, require_root:bool=True, force=False):
        status = True

        # First, we check if HepMC Python bindings are available, and have what we need
        try:
            assert self._basic_check()
        except:
            if(self.verbose): self._print('pyHepMC3 not found.')
            status = False

        # If we need the ROOT interface, we'll explicitly check that too.
        # The HepMC3 available via CVMFS has Python bindings, but might specifically lack this.
        if(require_root):
            try:
                assert self._root_check()
            except:
                if(self.verbose): self._print('pyHepMC3 ROOT interface requested, but not found.')
                status = False

        if(status and not force):
            # NOTE: Would be best to fix this so that self.hepmc_dir points to existing install
            return

        self.SetDirectory()
        self.DownloadHepMC3()
        self.MakeHepMC3(j)
        self.SetPythonDirectory()

    def _basic_check(self):
        try:
            from pyHepMC3 import HepMC3 as hm
            evt = hm.GenEvent()
        except:
            return False
        return True

    def _root_check(self):
        try:
            import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hmroot
            reader_class = hmroot.ReaderRoot
        except:
            return False
        return True

    def _print(self,val:str):
        print('{}: {}'.format(self.prefix,val))

    def run_make_with_progress(self,command=['make'],prefix='Building HepMC'):
        # Pattern to match progress indicators
        progress_pattern = re.compile(r'\[\s*(\d+)%\]')

        with open(self.logfile, 'w') as f, open(self.errfile, 'w') as g:
            process = sub.Popen(
                command,
                cwd=self.build_dir,
                stdout=sub.PIPE,
                stderr=sub.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )

            # Read stdout line by line
            for line in iter(process.stdout.readline, ''):
                f.write(line)  # Write to log file
                f.flush()

                # Check for progress indicator
                match = progress_pattern.search(line)
                if match:
                    progress = int(match.group(1))
                    printProgressBar(progress, 100, prefix=prefix, suffix='Complete')

            # Read any remaining stderr
            stderr_output = process.stderr.read()
            if stderr_output:
                g.write(stderr_output)

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise sub.CalledProcessError(return_code, command)
