
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

class DelphesSimulator(DetectorSimulator):

    def __init__(self,configurator,output_directory):
        super().__init__()
        self.SetConfigurator(configurator)
        self.delphes_card = self.configurator.GetDelphesCard()
        self.delphes_wrapper = DelphesWrapper(self.configurator.GetDelphesDirectory())
        self.SetOutputDirectory(output_directory)
        self.Initialize()

    def Initialize(self):
        self.delphes_wrapper.PrepDelphes() # this will download/build Delphes if necessary

    def Process(self,files=None):
        if(files is not None):
            self.SetInputs(files)

        for i,hep_file in enumerate(self.input_files):
            delphes_file = hep_file.replace('.hepmc','.root')
            print('Running DelphesHepMC3: {} -> {}.'.format(hep_file,delphes_file))
            if(i == 0):
                print('\tDelphes executable: {}'.format(self.delphes_wrapper.GetExecutable()))
                print('\tDelphes card: {}'.format(self.delphes_card))
                delphes_file = self.delphes_wrapper.HepMC3ToDelphes(
                    hepmc_file=hep_file,
                    output_file=delphes_file,
                    cwd=self.outdir,
                    delphes_card=self.delphes_card
                )
                self.output_files.append(delphes_file)
        return
