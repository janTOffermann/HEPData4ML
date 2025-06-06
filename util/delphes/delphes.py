import os, glob, pathlib
import subprocess as sub
from util.qol_utils.misc import stdout_redirected
from util.delphes.setup import DelphesSetup
from typing import Union,Optional

class DelphesWrapper:
    """
    A class for running the
    Delphes fast detector simulation framework.
    """
    def __init__(self,delphes_dir=None):

        self.setup = DelphesSetup(delphes_dir=delphes_dir)
        self.setup.PrepDelphes()
        self.delphes_dir = self.setup.GetDirectory()
        self.executable = self.setup.GetExecutable()
        self.default_card = "delphes_card_ATLAS.tcl" # Default to the ATLAS card. Will search for this within the Delphes directory.

    def GetExecutable(self)->str:
        return self.executable

    def HepMC3ToDelphes(self,hepmc_file, output_file=None, delphes_card=None, logfile=None, cwd=None, force=False):
        hepmc_file_nodir = hepmc_file
        if(cwd is not None):
            hepmc_file = '{}/{}'.format(cwd,hepmc_file)

        if(output_file is not None and cwd is not None):
            output_file = '{}/{}'.format(cwd,output_file)

        # default to using HepMC filename for ROOT output
        if(output_file is None):
            output_file = hepmc_file.replace('.hepmc','.root')

        if(cwd is not None): output_file_nodir = hepmc_file_nodir.replace('.hepmc','.root')
        else: output_file_nodir = output_file.split('/')[-1]

        # check if the output file already exists (i.e. look for file with same name)
        if(pathlib.Path(output_file).exists() and not force):
            print('\t\tDelphes ROOT file {} already found, skipping its generation.'.format(output_file))
            if(cwd is not None): return output_file_nodir
            return output_file

        # delphes_ex = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)
        # if(len(delphes_ex) == 0):
        #     print('Error: Delphes executable not found at {}/**?DelphesHepMC3 .'.format(self.delphes_dir))
        #     assert(False)
        # delphes_ex = delphes_ex[0]

        delphes_ex = self.executable

        # default card
        if(delphes_card is None): delphes_card = glob.glob('{}/**/{}'.format(self.delphes_dir,self.default_card),recursive=True)[0]

        # Delphes will crash if output file already exists, so we need to remove it.
        try: os.remove(output_file)
        except: pass

        if(logfile is not None):
            with open(logfile,'w') as f:
                sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                            shell=False, stdout=f, stderr=f)

        else:
            sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                        shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)

        if(cwd is not None): return output_file_nodir
        return output_file # return the name (esp. useful if none was provided)