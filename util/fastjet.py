import os, glob,pathlib
import subprocess as sub
import util.qol_utils.qol_util as qu

class FastJetSetup:
    def __init__(self,fastjet_dir=None, full_setup=False,verbose=True):
        self.fastjet_version = '3.4.0'
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
            with qu.stdout_redirected():
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