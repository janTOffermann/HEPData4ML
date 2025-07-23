
import re, json
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

        self.delphes_card_text = {}

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

    def FetchDelphesCard(self):
        if(self.delphes_card is not None):
            with open(self.delphes_card,'r') as f:
                lines = f.readlines()
            lines = ''.join(lines)
            self.delphes_card_text[self.delphes_card.split('/')[-1]] = lines

        # also check for imported tcl files in the delphes card -- if there are any, add them too.
        pattern_with_capture = r'source\s+(\S*\.tcl)'
        filenames = re.findall(pattern_with_capture, lines)

        for filename in filenames:
            card_dir = '/'.join(self.delphes_card.split('/')[:-1])
            filename_full = '{}/{}'.format(card_dir,filename)
            with open(filename_full,'r') as f:
                lines = f.readlines()
            lines = ''.join(lines)
            self.delphes_card_text[filename] = lines
        return

    def GetDelphesCard(self):
        """
        Serialize delphes_card_text dictionary with JSON,
        for putting into HDF5 attributes. (dict not supported)
        """
        return json.dumps(self.delphes_card_text)