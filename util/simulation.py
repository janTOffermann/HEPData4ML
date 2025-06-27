
import numpy as np
from util.config import Configurator
from util.delphes.delphes import DelphesWrapper

from typing import TYPE_CHECKING
if TYPE_CHECKING: # Only imported during type checking -- limits the risk of circular imports, as the code is further developped
    from util.config import Configurator

class DetectorSimulator:
    def __init__(self):
        self.input_files = []
        self.output_files = []

    def SetConfigurator(self,configurator : 'Configurator'):
        self.configurator = configurator

    def SetInputs(self,files:list):
        self.input_files = files

    def SetOutputDirectory(self,directory:str):
        self.outdir = directory

    def GetOutputFiles(self):
        return self.output_files

    def SetLogfile(self,file:str):
        self.logfile = file

class DelphesSimulator(DetectorSimulator):

    def __init__(self,configurator,output_directory,logfile=None):
        super().__init__()
        self.SetConfigurator(configurator)
        self.delphes_card = self.configurator.GetDelphesCard()

        use_root = self.configurator.GetHepMCFormat().lower() == 'root'
        self.delphes_wrapper = DelphesWrapper(self.configurator.GetDelphesDirectory(),use_root=use_root)
        self.SetOutputDirectory(output_directory)
        self.SetLogfile(logfile)
        self.mode = None

    def SetMode(self,mode:str):
        self.mode = mode.lower()

        # if(mode == 'root'): # this requires special setup TODO: Move upstream to initialization?
        #     self.delphes_wrapper._setup_root()

    def Process(self,files=None):
        if(files is not None):
            self.SetInputs(files)

        for i,hep_file in enumerate(self.input_files):

            if(hep_file.split('.')[-1].lower() == 'root'):
                self.SetMode('root')
            else:
                self.SetMode('hepmc') # or 'ascii'?

            delphes_file = '.'.join(hep_file.split('.')[:-1]) + '_delphes.root'

            if(self.mode=='root'):
                print('Running DelphesHepMC3ROOT: {} -> {}.'.format(hep_file,delphes_file))
            else:
                print('Running DelphesHepMC3: {} -> {}.'.format(hep_file,delphes_file))
            if(i == 0):
                print('\tDelphes executable: {}'.format(self.delphes_wrapper.GetExecutable()[self.mode]))
                print('\tDelphes card: {}'.format(self.delphes_card))
            delphes_file = self.delphes_wrapper.HepMC3ToDelphes(
                hepmc_file=hep_file,
                output_file=delphes_file,
                cwd=self.outdir,
                delphes_card=self.delphes_card,
                logfile=self.logfile
            )
            self.output_files.append(delphes_file)
        return
