import os,sys,glob,pathlib
import subprocess as sub
from util.qol_utils.misc import stdout_redirected

class PythiaInstaller:
    """
    This class will simply fetch and install Pythia8.
    We do this because we need an installation with both
    Python3 bindings, and HepMC3 support. The Pythia8 available
    via LCG (on CVMFS) seems to lack the latter, unfortunately,
    and custom HepMC3 writing with pyhepmc is a mess with ISR/FSR.
    """

    def __init__(self):
        self.pythia_download = 'https://pythia.org/download/pythia83/pythia8315.tgz' # LCG_105 - LCG_107 have 3.5.1pre09
        self.pythia_file = 'pythia8315'

    def SetDirectory(self,directory=None):
        self.directory = directory
        if(self.directory is None):
            self.directory = os.path.dirname(os.path.abspath(__file__)) + '/../../pythia'
        self.directory = os.path.normpath(self.directory)

        # Make the Pythia dir if it does not exist.
        os.makedirs(self.directory,exist_ok=True)

        # Files that are used for logging progress with downloading/building Pythia.
        self.logfile = '{}/log.stdout'.format(self.directory)
        self.errfile = '{}/log.stderr'.format(self.directory)

    def PrepPythia(self, j=4, force=False, verbose=False):
        """
        Checks if Pythia8 is built. If it's not, or if it lacks the features we need,
        this will download and build it.
        """
        # Check if Pythia8 is already built at destination.
        if(not force):
            pass # TODO: Will be a bit involved, need to test Python3 bindings and HepMC3 support
            # ex = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)
            # if(len(ex) > 0):
            #     self.executable = ex[0]
            #     if(verbose): print('Found DelphesHepMC3 @ {}'.format(self.executable))
            #     return

        self.DownloadPythia() # TODO: Could check to see if source code is downloaded, though it's hard to make a perfect check.
        self.ConfigurePythia()
        self.BuildPythia(j=j)
        # self.executable = ex = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)[0]
        return

    def DownloadPythia(self):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            # Fetch Delphes source.
            print('Downloading Pythia.')
            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget):
                sub.check_call('wget --no-check-certificate --content-disposition {} -O {}'.format(self.pythia_download,self.pythia_file).split(' '), shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            else:
                sub.check_call('curl -LJ {} -o {}'.format(self.pythia_download,self.pythia_file).split(' '), shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['tar', '-zxf', self.pythia_file], shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['rm', self.pythia_file], shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
        return

    def ConfigurePythia(self):
        print('Configuring Pythia installation.')

    def BuildPythia(self, j=4):
        print('Making Pythia')

        # with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
        #     print('Making Delphes.')
        #     sub.check_call(['make', '-j{}'.format(j)],
        #                 shell=False, cwd = '{}/{}'.format(self.delphes_dir,self.pythia_file.replace('.tar.gz','')), stdout=f, stderr=g)
        return
