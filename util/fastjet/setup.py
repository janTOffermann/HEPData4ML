import os, glob,pathlib
import numpy as np
import subprocess as sub
from util.qol_utils.misc import stdout_redirected
from util.qol_utils.progress_bar import printProgressBar

class FastJetSetup:
    def __init__(self,fastjet_dir=None, full_setup=False,verbose=True):
        self.fastjet_version = '3.4.2'
        self.SetDirectory(fastjet_dir)
        if(full_setup): self.PrepFastjet(verbose=verbose)

    def SetDirectory(self,fastjet_dir=None):
        self.fastjet_dir = fastjet_dir
        if(self.fastjet_dir is None):
                self.fastjet_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/fastjet'
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
            if(verbose): print('Building fastjet.')

            self.run_make_with_progress(command=['make', '-j{}'.format(j)])

            # sub.check_call(['make', '-j{}'.format(j)],
            #             shell=False, cwd=self.source_dir, stdout=f, stderr=g)

            if(verbose): print('Installing fastjet.')
            sub.check_call(['make', 'install'],
                        shell=False, cwd=self.source_dir, stdout=f, stderr=g)
        self.SetPythonDirectory()
        return

    def run_make_with_progress(self,command=['make'],prefix='Building Fasatjet'):
        """
        Runs "make", with a progress bar.
        NOTE: This is quite hard-coded, due to how the make printouts
        work for Fastjet (they don't have their own progress tracking).
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
            N = 947
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