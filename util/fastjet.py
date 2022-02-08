import sys, os, glob
import subprocess as sub
from util.qol_util import printProgressBarColor

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
        if(0 not in files_to_find):
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
        print('Downloading fastjet from {}.'.format(fastjet_download))
        sub.check_call(['wget', fastjet_download], 
                       shell=False, cwd=fastjet_dir, stdout=f, stderr=g)
        sub.check_call(['tar', 'zxvf', 'fastjet-3.4.0.tar.gz'], 
                       shell=False, cwd=fastjet_dir, stdout=f, stderr=g)
        sub.check_call(['rm', 'fastjet-3.4.0.tar.gz'], 
                       shell=False, cwd=fastjet_dir, stdout=f, stderr=g)

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
    """illustrative class for use in assigning pythonic user information
    to a PseudoJet.
    """
    def __init__(self, particle_index, status, pdg_id=0):
        self.particle_index = particle_index
        self.status = status
        self.pdg_id = pdg_id

    def __str__(self):
        return "particle_index={0}, status={1}, pdg_id={2}".format(
            self.particle_index, self.status, self.pdg_id)