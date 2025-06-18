import os, glob, re
import subprocess as sub
import numpy as np
from util.qol_utils.misc import stdout_redirected
from util.qol_utils.progress_bar import printProgressBar
from typing import Optional

class DelphesSetup:
    """
    A class for running downloading/building (if necessary) the
    Delphes fast detector simulation framework.
    """
    def __init__(self,delphes_dir:Optional[str]=None):
        self.SetDirectory(delphes_dir)
        self.delphes_download = 'https://github.com/delphes/delphes/archive/refs/tags/3.5.1pre12.tar.gz' # LCG_105 - LCG_107 has 3.5.1pre09
        self.delphes_file = 'delphes-{}'.format(self.delphes_download.split('/')[-1]) # TODO: A bit hacky
        self.executable = None
        self.prefix = self.__class__.__name__

    def SetDirectory(self,delphes_dir:Optional[str]=None):
        self.delphes_dir = delphes_dir
        if(self.delphes_dir is None):
            self.delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/delphes'
        self.delphes_dir = os.path.normpath(self.delphes_dir)

        # Make the Delphes dir if it does not exist.
        os.makedirs(self.delphes_dir,exist_ok=True)

        # Files that are used for logging progress with downloading/building Delphes.
        self.logfile = '{}/log.stdout'.format(self.delphes_dir)
        self.errfile = '{}/log.stderr'.format(self.delphes_dir)

    def GetDirectory(self)->str:
        return self.delphes_dir

    def GetExecutable(self)->str:
        return self.executable

    def PrepDelphes(self, j:int=4, force:bool=False, verbose:bool=False):
        """
        Checks if Delphes is built. If it's not, it will download
        and build Delphes.
        """
        # Check if Delphes is already built at destination.
        if(not force):
            ex1 = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)
            if(len(ex1) > 0):
                self.executable = ex1[0]
                if(verbose): self._print('Found DelphesHepMC3 @ {}'.format(self.executable['hepmc']))
                return

        self.DownloadDelphes() # TODO: Could check to see if source code is downloaded, though it's hard to make a perfect check.
        self.BuildDelphes(j=j)
        self.executable = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)[0]
        return

    def DownloadDelphes(self):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            # Fetch Delphes source.
            self._print('Downloading Delphes from {} .'.format(self.delphes_download))
            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget):
                sub.check_call('wget --no-check-certificate --content-disposition {} -O {}'.format(self.delphes_download,self.delphes_file).split(' '),cwd=self.delphes_dir, stdout=f, stderr=g)
            else:
                sub.check_call('curl -LJ {} -o {}'.format(self.delphes_download,self.delphes_file).split(' '),cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['tar', '-zxf', self.delphes_file],cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['rm', self.delphes_file],cwd=self.delphes_dir, stdout=f, stderr=g)
        return

    def BuildDelphes(self, j:int=4):
        self.build_dir = '{}/{}'.format(self.delphes_dir,self.delphes_file.replace('.tar.gz','')) # TODO: clean this up
        self._print('Building Delphes @ {} .'.format(self.build_dir))
        command = ['make', '-j{}'.format(j)]
        self.run_make_with_progress(command)
        return

    def _print(self,val:str):
        print('{}: {}'.format(self.prefix,val))

    def run_make_with_progress(self,command=['make'],prefix='Building Delphes'):
        """
        Runs "make", with a progress bar.
        NOTE: This is quite hard-coded, due to how the make printouts
        work for Delphes (they don't have their own progress tracking).
        """

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
            counter = 0
            N = 281
            for line in iter(process.stdout.readline, ''):
                f.write(line)  # Write to log file
                f.flush()

                progress = (counter + 1)/N * 100.
                progress = np.maximum(progress,100) # just in case -- this is all a little hard-coded
                printProgressBar(progress, 100, prefix=prefix, suffix='Complete')
                counter +=1

            # Read any remaining stderr
            stderr_output = process.stderr.read()
            if stderr_output:
                g.write(stderr_output)

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise sub.CalledProcessError(return_code, command)

class DelphesROOTHepMC3Setup:
    """
    A class for running building the
    DelphesHepMC3ROOT executable -- our own custom one
    for handling HepMC3/ROOT files.
    """
    def __init__(self):
        self.SetDirectory()
        self.executable = None
        self.prefix = self.__class__.__name__

    def SetDirectory(self,delphes_dir:Optional[str]=None):
        self.delphes_dir = delphes_dir
        if(self.delphes_dir is None):
            self.delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/delphes_custom'
        self.delphes_dir = os.path.normpath(self.delphes_dir)

        # # Make the Delphes dir if it does not exist. # Commented out for now, while delphes_custom is part of this package directly
        # os.makedirs(self.delphes_dir,exist_ok=True)

        # Files that are used for logging progress with building DelphesHepMC3ROOT.
        self.logfile = '{}/log.stdout'.format(self.delphes_dir)
        self.errfile = '{}/log.stderr'.format(self.delphes_dir)

    def GetDirectory(self)->str:
        return self.delphes_dir

    def GetExecutable(self)->str:
        return self.executable

    def PrepDelphesHepMC3ROOT(self, j:int=4, force:bool=False, verbose:bool=False):
        """
        Checks if DelphesHepMC3ROOT is built. If it's not,
        it will build it.
        """
        # Check if DelphesHepMC3ROOT is already built at destination.
        if(not force):
            ex1 = glob.glob('{}/**/DelphesHepMC3ROOT'.format(self.delphes_dir),recursive=True)
            if(len(ex1) > 0):
                self.executable = ex1[0]
                if(verbose): self._print('Found DelphesHepMC3ROOT @ {}'.format(self.executable))
                return


        self.Build(j=j)
        self.executable = glob.glob('{}/**/DelphesHepMC3ROOT'.format(self.delphes_dir),recursive=True)[0]
        return

    def Build(self, j:int=4):
        self.build_dir = '{}/{}'.format(self.delphes_dir,'build')
        self._print('Building DelphesHepMC3ROOT @ {} .'.format(self.build_dir))
        command = ['cmake', '../']
        sub.check_call(command,cwd=self.build_dir)
        command = ['make', '-j{}'.format(j)]
        self.run_make_with_progress(command)
        return

    def _print(self,val:str):
        print('{}: {}'.format(self.prefix,val))

    def run_make_with_progress(self,command=['make'],prefix='Building DelphesHepMC3ROOT'):
        """
        Runs "make", with a progress bar.
        """
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
