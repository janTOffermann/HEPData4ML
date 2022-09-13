import sys, os, glob
import subprocess as sub
import util.qol_utils.qol_util as qu

def BuildFastjet(fastjet_dir=None, j=4, force=False, verbose=False):
    if(fastjet_dir is None):
        fastjet_dir = os.path.dirname(os.path.abspath(__file__)) + '/../fastjet'

    # Check if Fastjet is already built at destination.
    # Specifically, we will look for some Python-related files.
    if(not force):

        files_to_find = [
            '{}/**/site-packages/fastjet.py'.format(fastjet_dir),
            '{}/**/site-packages/_fastjet.a'.format(fastjet_dir),
            '{}/**/site-packages/_fastjet.so.0'.format(fastjet_dir)
        ]

        files_to_find = [glob.glob(x,recursive=True) for x in files_to_find]
        files_to_find = [len(x) for x in files_to_find]
        found_fastjet = False
        for entry in files_to_find:
            if(entry > 0):
                found_fastjet = True # Loose condition since not all these files may exist (e.g. no .so on macos). Should be okay.
                break
        if(found_fastjet):
            if(verbose): print('Found existing Fastjet installation with Python extension @ {}.'.format(fastjet_dir))
            return fastjet_dir

    # Make the Fastjet dir if it does not exist.
    try: os.makedirs(fastjet_dir)
    except: pass

    # Put output into log files.
    logfile = '{}/log.stdout'.format(fastjet_dir)
    errfile = '{}/log.stderr'.format(fastjet_dir)

    with open(logfile,'w') as f, open(errfile,'w') as g:

        # Fetch the Fastjet source
        fastjet_download = 'http://fastjet.fr/repo/fastjet-3.4.0.tar.gz'
        fastjet_file = fastjet_download.split('/')[-1]
        print('Downloading fastjet from {}.'.format(fastjet_download))

        # Depending on Linux/macOS, we use wget or curl.
        has_wget = True
        with qu.stdout_redirected():
            try: sub.check_call('which wget'.split(' '))
            except:
                has_wget = False
                pass

        if(has_wget): sub.check_call(['wget', fastjet_download], shell=False, cwd=fastjet_dir, stdout=f, stderr=g)
        else: sub.check_call(['curl',fastjet_download,'-o',fastjet_file], shell=False, cwd=fastjet_dir, stdout=f, stderr=g)
        sub.check_call(['tar', 'zxvf', fastjet_file], shell=False, cwd=fastjet_dir, stdout=f, stderr=g)
        sub.check_call(['rm', fastjet_file], shell=False, cwd=fastjet_dir, stdout=f, stderr=g)

        source_dir  = '{}/fastjet-3.4.0'.format(fastjet_dir)
        install_dir = '{}/fastjet-install'.format(fastjet_dir)

        # Now configure. We create the python bindings.
        print('Configuring fastjet.')
        sub.check_call(['./configure', '--prefix={}'.format(install_dir), '--enable-pyext'],
                       shell=False, cwd=source_dir, stdout=f, stderr=g)

        # Now make and install. Will skip "make check".
        print('Making fastjet.')
        sub.check_call(['make', '-j{}'.format(j)],
                       shell=False, cwd=source_dir, stdout=f, stderr=g)
        print('Installing fastjet.')
        sub.check_call(['make', 'install'],
                       shell=False, cwd=source_dir, stdout=f, stderr=g)
    return install_dir

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