
import re, json
import numpy as np
from util.config import Configurator
from util.delphes.delphes import DelphesWrapper

from typing import TYPE_CHECKING
if TYPE_CHECKING: # Only imported during type checking -- limits the risk of circular imports, as the code is further developped
    from util.config import Configurator
    from util.meta import MetaDataHandler

class DetectorSimulator:
    """
    Abstract parent class for detector simulation.
    """
    def __init__(self):
        self.input_files = []
        self.output_files = []

        # for metadata
        self.metadata_handler = None
        self.citations = {}

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

    def GetCitations(self):
        return self.citations

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

class DelphesSimulator(DetectorSimulator):
    """
    Detector (pseudo)simulation using the Delphes fast detector
    simulation framework.
    """
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
        self._generate_citations()

    def _generate_citations(self):
        """
        Fills in citations (in BibTex format) for Delphes.
        """
        key = 'Delphes'
        self.citations[key] = [
            """
@article{deFavereau:2013fsa,
    author = "de Favereau, J. and Delaere, C. and Demin, P. and Giammanco, A. and Lema{\\^\\i}tre, V. and Mertens, A. and Selvaggi, M.",
    collaboration = "DELPHES 3",
    title = "{DELPHES 3, A modular framework for fast simulation of a generic collider experiment}",
    eprint = "1307.6346",
    archivePrefix = "arXiv",
    primaryClass = "hep-ex",
    doi = "10.1007/JHEP02(2014)057",
    journal = "JHEP",
    volume = "02",
    pages = "057",
    year = "2014"
}
            """,
            """
@article{Selvaggi:2014mya,
    author = "Selvaggi, Michele",
    editor = "Wang, Jianxiong",
    title = "{DELPHES 3: A modular framework for fast-simulation of generic collider experiments}",
    doi = "10.1088/1742-6596/523/1/012033",
    journal = "J. Phys. Conf. Ser.",
    volume = "523",
    pages = "012033",
    year = "2014"
}
            """,
            """
@article{Mertens:2015kba,
    author = "Mertens, Alexandre",
    editor = "Fiala, L. and Lokajicek, M. and Tumova, N.",
    title = "{New features in Delphes 3}",
    doi = "10.1088/1742-6596/608/1/012045",
    journal = "J. Phys. Conf. Ser.",
    volume = "608",
    number = "1",
    pages = "012045",
    year = "2015"
}
            """
        ]

    def _writeMetadata(self):
        if(self.metadata_handler is None):
            self._print('Unable to write metadata; no handler was provided.')
            return

        # Add information on Delphes card
        key = 'Metadata.Simulation.DelphesCard'
        self.FetchDelphesCard()
        self.metadata_handler.AddElement(key,self.GetDelphesCard())

        # Also add metadata on citations for algorithms.
        self.metadata_handler.AddCitations(self.GetCitations())

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

        self._writeMetadata()
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
        return self.delphes_card_text
        # return json.dumps(self.delphes_card_text)