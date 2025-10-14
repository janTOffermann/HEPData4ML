import glob, itertools
import numpy as np
import ROOT as rt
from util.reconstruction.setup import NtupleProducerSetup
from util.buffer.input import UprootTreeLoader
from util.buffer.output import OutputBuffer, RootOutputBuffer
from util.qol_utils.progress_bar import printProgressBarColor
from util.hepmc.hepmc import ExtractHepMCEvents, ExtractHepMCParticles, ParticleToProductionVertex, ParticleToEndVertex, ParticleToMomenta
from typing import Union, Optional, List, Any, TYPE_CHECKING
from util.misc.timing import profile_method, profile_block
import util.hdf5.hdf5 as h5util
import util.root.utils as rootutil


if(TYPE_CHECKING):
    import sys
    from util.config.config import Configurator
    from util.hepmc.setup import HepMCSetup
    setup = HepMCSetup(verbose=False)
    python_dir = setup.GetPythonDirectory()
    if(python_dir not in sys.path):
        sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm
    from util.metadata.meta import MetaDataHandler

class Processor:
    """
    Convert HepMC3 (+ optionally Delphes) file,
    into a ROOT n-tuple file.
    """
    def __init__(self, configurator:'Configurator'):
        self.configurator = configurator

        # Make sure that the NtupleProducer library is built and ready-to-go.
        self.setup = NtupleProducerSetup()
        self.setup.FullPreparation()
        self._init_ntuple_producer()

        self.print_prefix = 'Processor'

        # Set the output format.
        self.output_format = 'root'
        self.output_extension = 'root'

        self.SetProgressBarPrefix('Converting HepMC3 -> N-tuple:')
        self.suffix = 'Complete'
        self.bl = 50
        self.verbose = False

        self.outdir = ''

        self.SetStableTruthParticleName('StableTruthParticles') # a "special" name for the stable truth particles collection; this is always present
        self.SetTreeName('hepdata4ml_tree')

        self.SetPostProcessing()

        # Various integers for buffer size, number of particles read into memory from HepMC, number saved to file, etc.
        self.nevents = None

        # truth selector
        self.SetParticleSelection()

        # Metadata handling
        self.metadata_handler = None

        # Compression level for HDF5, 0 (least) to 9 (most)
        self.copts = 0

    def SetTreeName(self,val):
        self.tree_name = val
        if(self.processor is None):
            self._init_ntuple_producer()
        self.processor.SetTreeName(self.tree_name)

    def GetTreeName(self):
        return self.tree_name

    def SetStableTruthParticleName(self,val):
        self.stable_truth_particle_name = val
        if(self.processor is None):
            self._init_ntuple_producer()
        self.processor.SetStableTruthParticleName(self.stable_truth_particle_name)

    def _init_ntuple_producer(self):
        self.processor = rt.NtupleProducer.Converter()

        for entry in self.configurator.GetDelphesObjects():
            self.processor.AddDelphesObject(entry)

    def _set_output_format(self):
        self.output_format = self.configurator.GetReconstructionOutputFormat()

        self.output_format = 'root' # TODO: Only doing ROOT now (HDF5 conversion can come later)

        if(self.output_format.lower() == 'hdf5'):
            self.output_extension = 'h5'
        elif(self.output_format.lower() == 'root'):
            self.output_extension = 'root'
        else:
            self._print('Error: Output format {} not recognized.'.format(self.output_format))
            return

    def GetOutputExtension(self):
        return self.output_extension

    def SetH5Compression(self,val:int):
        if(val > 9):
            val = 9
        elif(val < 0):
            val = 0
        self.copts = val

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def SetParticleSelection(self):
        if(self.configurator is None):
            return
        else:
            self.truth_selection = self.configurator.GetParticleSelection()

        if(self.processor is None):
            self._init_ntuple_producer()

        if(self.truth_selection is not None):
            for key,val in self.truth_selection.items():
                self.processor.AddTruthParticleSelector(key,val.GetSelector())

    def SetProgressBarPrefix(self,text:str):
        self.prefix_level1 = text

    def SetPostProcessing(self,post_proc=None):
        if(post_proc is None): post_proc = self.configurator.GetPostProcessing()
        self.post_processing = post_proc

    def SetStatsFile(self,filename:str):
        self.stats_file = filename

    def SetVerbosity(self,flag:bool):
        self.verbose = flag

    def SetOutputDirectory(self,outdir:str):
        self.outdir = outdir

    def ProcessFull(self, hepmc_file:str, delphes_file:str, output_file:Optional[str]=None, verbosity:int=0):

        output_file = self.Process(hepmc_file,delphes_file, output_file,verbosity)
        self.PostProcess(hepmc_file,output_file)
        return output_file

    @profile_method('Processor.Process')
    def Process(self, hepmc_file:str, delphes_file:Optional[str], output_file:Optional[str]=None, verbosity:int=0):

        # Parse hepmc_files
        hepmc_file = '{}/{}'.format(self.outdir,hepmc_file)

        # Parse delphes_files
        if(delphes_file is not None):
            delphes_file = '{}/{}'.format(self.outdir,delphes_file)
        else:
            delphes_file = ""

        # Parser output_file
        if(output_file is None):
            output_file = 'events.hepdata4ml.{}'.format(self.output_extension)
        output_file = '{}/{}'.format(self.outdir,output_file)

        # TODO: Run the process here
        self.processor.Process(hepmc_file, delphes_file, output_file) # <- does the whole n-tuple conversion (to ROOT format)
        return

    @profile_method('Processor.PostProcess')
    def PostProcess(self,hepmc_file:str, ntuple_file:str, prepend_output_directory=False):
        if(self.post_processing is None):
            return
        for post_proc in self.post_processing:
            if(post_proc is None): continue
            post_proc.SetConfigurator(self.configurator)
            post_proc.SetMetadataHandler(self.metadata_handler)
            if(prepend_output_directory): # TODO: Clean this up? A little inconsistent
                ntuple_file = '{}/{}'.format(self.outdir,ntuple_file)
                hepmc_file = '{}/{}'.format(self.outdir,hepmc_file)
            post_proc(hepmc_file,ntuple_file,ntuple_file)
        return

    ################################################################
    # Some member functions that will behave like utilities,
    # i.e. we'll call them separate from the main Processor routine.
    ################################################################

    def ConcatenateNtuples(self,input_files:list,output_file, format=None):
            """
            This is a utility function for concatenating output ntuples.
            """

            if(format is None):
                format = self.output_format

            output_filepath = '/'.join((self.outdir,output_file))
            if(format.lower() == 'hdf5'):
                compression_opts = 1 # TODO fetch from configuration?
                h5util.ConcatenateH5(input_files,output_filepath,copts=compression_opts,delete_inputs=True,ignore_keys=['Event.Index'],verbose=False,silent_drop=True)
            elif(format.lower() == 'root'):
                rootutil.ConcatenateRootTreeFiles(input_files,output_filepath,self.tree_name,['Event.Index'])
            else:
                self._print('Error: ConcatenateNtuples() not implemented for format {}.'.format(self.output_format))
                return

    def AddEventIndices(self,input_file,key='Event.Index',offset=0):
        key = 'Event.Index' # TODO: Make this member variable?
        if(self.output_format.lower() == 'hdf5'):
            compression_opts = 1 # TODO fetch from configuration?
            h5util.AddEventIndices(input_file,cwd=self.outdir,copts=compression_opts,key=key,offset=offset)
        elif(self.output_format.lower() == 'root'):
            rootutil.AddEventIndices(input_file,self.tree_name,cwd=self.outdir,index_branch_name='Event.Index',offset=offset)
        else:
            self._print('Error: AddEventIndices() not implemented for format {}.'.format(self.output_format))
            return

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return
