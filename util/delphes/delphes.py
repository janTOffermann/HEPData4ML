import os, glob, pathlib
import subprocess as sub
from util.qol_utils.misc import stdout_redirected
from util.delphes.setup import DelphesSetup, DelphesROOTHepMC3Setup
from typing import Union,Optional

class DelphesWrapper:
    """
    A class for running the
    Delphes fast detector simulation framework.
    """
    def __init__(self,delphes_dir=None, use_root=False):

        self.setup = DelphesSetup(delphes_dir=delphes_dir)
        self.setup.Prepare()
        self.delphes_dir = self.setup.GetDirectory()
        self.executable = {'hepmc': self.setup.GetExecutable()}
        self.default_card = "delphes_card_ATLAS.tcl" # Default to the ATLAS card. Will search for this within the Delphes directory.

        # Also prepare an instance of the setup class for our custom executable,
        # for handling HepMC3/ROOT input.
        self.setup_hepmc3root = DelphesROOTHepMC3Setup() # its delphes_dir is fixed
        if(use_root):
            self._setup_root()

    def _setup_root(self):
        self.setup_hepmc3root.Prepare()
        self.delphes_root_dir = self.setup_hepmc3root.GetDirectory()
        self.executable['root'] = self.setup_hepmc3root.GetExecutable()

    def GetExecutable(self)->str:
        return self.executable

    def HepMC3ToDelphes(self,hepmc_file, output_file=None, delphes_card=None, logfile=None, cwd=None, force=False):
        """
        This function runs the Delphes executable,
        to convert HepMC3(ROOT) files to DELPHES ROOT files.
        (It supports both plaintext and ROOT-based HepMC3).
        """
        if(cwd is not None):
            hepmc_file = '{}/{}'.format(cwd,hepmc_file)

        if(output_file is not None and cwd is not None):
            output_file = '{}/{}'.format(cwd,output_file)

        # default to using HepMC filename for ROOT output
        if(output_file is None):
            output_file = hepmc_file.replace('.hepmc','.root')

        # if(cwd is not None): output_file_nodir = hepmc_file_nodir.replace('.hepmc','.root')
        output_file_nodir = output_file.split('/')[-1]

        # check if the output file already exists (i.e. look for file with same name)
        if(pathlib.Path(output_file).exists() and not force):
            print('\t\tDelphes ROOT file {} already found, skipping its generation.'.format(output_file))
            if(cwd is not None): return output_file_nodir
            return output_file


        delphes_ex = self.executable['hepmc']
        if(hepmc_file.split('.')[-1].lower() == 'root'):
            delphes_ex = self.executable['root']

        # default card
        if(delphes_card is None): delphes_card = glob.glob('{}/**/{}'.format(self.delphes_dir,self.default_card),recursive=True)[0]

        # Delphes will crash if output file already exists, so we need to remove it.
        try: os.remove(output_file)
        except: pass

        if(logfile is not None):
            with open(logfile,'a') as f:
                sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                            shell=False, stdout=f, stderr=f)

        else:
            sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                        shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)

        if(cwd is not None): return output_file_nodir
        return output_file # return the name (esp. useful if none was provided)