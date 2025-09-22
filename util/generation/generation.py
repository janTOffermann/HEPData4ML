import os
import numpy as np
import subprocess as sub
from util.pythia.utils import PythiaWrapper
from util.hepmc.hepmc import Pythia8HepMC3Writer
from util.hepmc.Pythia8ToHepMC3 import PythiaToHepMC, PythiaToHepMCBatch
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from util.qol_utils.progress_bar import printProgressBarColor
from util.misc.timing import profile_method, profile_block
from typing import Optional,TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports, limits unnecessary imports
    from util.config.config import Configurator
    from util.metadata.meta import MetaDataHandler

class PythiaGenerator:
    """
    Generates events using Pythia8.
    """

    # TODO: Make a more general parent class, and use inheritance?

    def __init__(self, pt_min:float, pt_max:float, configurator:'Configurator', pythia_rng:Optional[int]=None, pythia_config_file:Optional[str]=None):
        self.configurator = configurator
        self.pt_min = pt_min
        self.pt_max = pt_max

        # Create our Pythia wrapper.
        self.pythia_rng = pythia_rng
        self.verbose = self.configurator.GetPythiaVerbosity()
        self.pythia = PythiaWrapper(verbose=self.verbose)
        self.ConfigPythia(config_file=pythia_config_file,verbose=self.verbose)

        # Set up HepMC, and create our HepMC converter
        self.hepmc_setup = HepMCSetup(self.configurator.GetHepMC3Directory(),verbose=False)
        # self.hepmc_setup.PrepHepMC()
        python_dir = self.hepmc_setup.GetPythonDirectory()
        # uncache_hepmc3()
        prepend_to_pythonpath(python_dir)

        # Also set the configurator's HepMC directory, so that in case it was "None" we don't end up
        # downloading HepMC3 multiple times.
        self.configurator.SetHepMC3Directory(self.hepmc_setup.GetDirectory())

        self.hepmc_converter = PythiaToHepMC(self.configurator.GetHepMC3Directory())
        self.hepmc_converter_batch = PythiaToHepMCBatch(self.configurator.GetHepMC3Directory())


        # Event filters. # TODO: May remove
        self.event_filter = None

        # Things for the progress bar.
        if(self.pt_min > 0. and self.pt_max > 0.):
            self.prefix = 'Generating events for pT^ in [{}, {}] GeV:'.format(self.pt_min,self.pt_max)
        elif(self.pt_min <= 0. and self.pt_max > 0.):
            self.prefix = 'Generating events for pT^ < {} GeV:'.format(self.pt_max)
        else:
            self.prefix = 'Generating events (default pT^ boundaries):'

        self.suffix = 'Complete'
        self.bl = 50

        self.SetOutputDirectory()

        self.writer = Pythia8HepMC3Writer()

        self.SetFilename('events.hepmc')
        self.filename_fullpath = None

        # self.diagnostic_plots = True
        # self.InitializeHistograms()

        # Containers for event-level information.
        self.xsecs = None

        self.progress_bar = True

        self.loop_number = 0 # used for keeping track of successful generation loops
        self.nevents = None # number of events requested, will be set in generation function

        self.hepev_buffer = []
        self.SetBufferSize(100) # TODO: Should this be configurable? Could be too much detail.
        self.buffername = None
        self.buffername_truth = None
        self.nevents_success = 0 # number of events successfully generated
        self.nevents_failed = 0 # number of events that failed (failure to pass basic cuts)

        self.header_status = False
        self.footer_status = False

        self.metadata_handler = None

        self.citations = {}
        self._generate_citations()

    def SetPtMin(self,pt_min):
        self.pt_min = pt_min

    def SetPtMax(self,pt_max):
        self.pt_max = pt_max

    def SetPt(self,pt_min,pt_max):
        self.SetPtMin(pt_min)
        self.SetPtMax(pt_max)

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def SetVerbose(self,val:bool):
        self.verbose = val
        self.pythia.SetVerbose(val)

    def SetEventFilter(self,filter):
        self.event_filter = filter

        if(filter is None):
            self.event_filter = self.configurator.GetEventFilter()

        if(self.event_filter is not None):
            self.event_filter.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import

    def SetPythiaConfigFile(self,file:Optional[str]=None):
        self.pythia_config_file = file
        if(file is None):
            self.pythia_config_file = self.configurator.GetPythiaConfigFile() # picked up from dictionary in config/config.py

    def SetDefaultFilenames(self):
        self.SetOutputDirectory(dir)
        self.SetFilename('events.hepmc')

    def SetProgressBar(self,flag:bool):
        self.progress_bar = flag

    def SetOutputDirectory(self,dir:Optional[str]=None):
        if(dir is None): dir = os.getcwd()
        self.outdir = dir

    def SetFilename(self,name:str):
        self.filename = name

        if(self.writer is not None):
            self.writer.SetFilename('{}/{}'.format(self.outdir,self.filename))
        return

    def GetFilename(self):
        return self.filename

    def GetHistFilename(self):
        return self.hist_filename

    def _generate_citations(self):
        """
        Fills in citations (in BibTex format) for Pythia8.
        """
        self.citations = {}
        key = 'Pythia8'
        self.citations[key] = """@article{Bierlich:2022pfr,
    author = "Bierlich, Christian and others",
    title = "{A comprehensive guide to the physics and usage of PYTHIA 8.3}",
    eprint = "2203.11601",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "LU-TP 22-16, MCNET-22-04, FERMILAB-PUB-22-227-SCD",
    doi = "10.21468/SciPostPhysCodeb.8",
    journal = "SciPost Phys. Codeb.",
    volume = "2022",
    pages = "8",
    year = "2022"
}"""

    def GetCitations(self):
        return self.citations

    def ConfigPythia(self, config_file:Optional[str]=None, verbose:bool=False):
        """
        Prepare and apply Pythia configuration. This turns our settings (from our config file)
        into a list of strings ready to be input to Pythia8.
        """
        self.pythia_config = self.configurator.GetPythiaConfig(self.pt_min,self.pt_max, verbose)
        self.SetPythiaConfigFile(file=config_file) # sets from the Python config file if "config_file = None"

        if(self.pythia_config_file is None): # special case
            return # not going to init the generator -- will need to handle this later on

        # Optionally set the Pythia RNG seed to something other than what's in the config.
        # TODO: Can we make this more tidy?
        if(self.pythia_rng is not None):
            self.pythia_config['Random:seed'] = self.pythia_rng

        # Now apply these configurations to our Generator's instance of PythiaWrapper.
        self.pythia.ClearConfigDict()
        self.pythia.AddToConfigDict(self.pythia_config)
        # self.pythia.ReadConfigDict()
        self.pythia.ReadStringsFromFile(self.pythia_config_file)
        self.pythia.InitializePythia()

    def SetBufferSize(self,size:int=100):
        self.buffer_size = size

    def GetCurrentBufferSize(self):
        return len(self.hepev_buffer)

    def ClearEventBuffer(self):
        self.hepev_buffer.clear()

    def AddToEventBuffer(self,hepev):
        if(isinstance(hepev,list)):
            self.hepev_buffer += hepev
            return
        self.hepev_buffer.append(hepev)

    def SetEventBuffer(self,list):
        self.hepev_buffer = list

    def WriteEventBufferToFile(self,header:bool=False,footer:bool=False):
        if(header): self.header_status = True
        if(footer): self.footer_status = True

        self.writer.Write(self.hepev_buffer)

        # PyHepMCOutput(self.hepev_buffer,self.buffername,self.filename_fullpath,header,footer)
        # HepMCOutputAscii(self.hepev_buffer,self.buffername,self.filename_fullpath,header,footer)

        self.ClearEventBuffer()

    def GenerateBatch(self, nevents, i_real:int=1, nevents_disp:Optional[int]=None):

        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that HepMC3's ASCII writing works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality, which
        # we can do with HepMC3's ROOT TTree format.
        self.filename_fullpath = '{}/{}'.format(self.outdir,self.filename)

        # For ASCII mode, create buffer file.
        if(self.writer.GetMode().lower() == 'ascii'):
            self.buffername = self.filename_fullpath.replace('.hepmc','_buffer.hepmc')

        # Determine how many events to actually generate on this call.
        # We base this on what was requested, but also on the current
        # buffer size.
        nevents_real = np.minimum(nevents,self.buffer_size - self.GetCurrentBufferSize())

        # Generate the events -- does the whole batch all at once!
        with profile_block('Generator.GenerateBatch - pythia'):
            self.pythia.GenerateBatch(nevents_real) # fills self.pythia.events

        # TODO: (Re)implement event filter logic.

        # Convert the Pythia8 events into HepMC3 events.
        hepmc_events = self.hepmc_converter_batch.fill_batch_events(self.pythia.events,self.pythia.GetPythia().infoPython(), i_real)

        # Fill the memory buffer with the event list.
        self.AddToEventBuffer(hepmc_events)

        # Write to the buffer, then flush it. Note that this is different than in GenerationLoop.
        self.WriteEventBufferToFile(header=True,footer=True)

        i_real += nevents_real # counter for number of successful events

        if(self.progress_bar): printProgressBarColor(i_real-1,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # Delete the buffer files, if relevant.
        if(self.buffername is not None):
            comm = ['rm', self.buffername]
            try: sub.check_call(comm,stderr=sub.DEVNULL)
            except: pass

        return i_real-1, 0 # note that i_real is using 1-indexing, which is what HepMC events use


    def GenerationLoop(self, nevents,i_real:int = 1, nevents_disp:Optional[int]=None):
        """
        This is the function where Pythia8 event generation happens, producing events.
        These are optionally filtered -- required to pass some condition(s) -- and then
        written to a HepMC3 file (either ASCII or ROOT format).

        This function operates as a Python loop. There is a batched version that uses
        Pythia8's batch generation and awkward arrays, which should be faster.
        """
        from pyHepMC3 import HepMC3 as hm # the HepMCSetup will have taken care of this -- so the package will be already cached
        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that HepMC3's ASCII writing works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality, which
        # we can do with HepMC3's ROOT TTree format.
        self.filename_fullpath = '{}/{}'.format(self.outdir,self.filename)

        # For ASCII mode, create buffer file.
        if(self.writer.GetMode().lower() == 'ascii'):
            self.buffername = self.filename_fullpath.replace('.hepmc','_buffer.hepmc')

        for i in range(nevents):

            if(self.pythia.IsInitialized()): # if false, pythia generator is not initialized -> will produce an empty HepMC event

                with profile_block('Generator.GenerationLoop - pythia'):
                    self.pythia.Generate() # generate an event!

                # ==========================================
                # Now we apply an (optional) "event filter", requiring that our event passes it.
                # If it does not, we will count the event as failed, and generate another.
                # ==========================================
                if(self.event_filter is not None):
                    passed_filter = self.event_filter(self.pythia)
                    if(not passed_filter):
                        n_fail += 1
                        continue

            # ==========================================
            # Now lets create the HepMC event.
            # ==========================================
            hepmc_event = hm.GenEvent()
            if(self.pythia.IsInitialized()):
                self.hepmc_converter.fill_next_event1(self.pythia.GetPythia(),hepmc_event,i_real)

            # Fill the memory buffer with this event.
            self.AddToEventBuffer(hepmc_event)

            # If buffer is full (i.e. has reached max size), write it to file & flush.
            if(self.GetCurrentBufferSize() == self.buffer_size):
                header = (self.loop_number == 0) and (not self.header_status)
                self.WriteEventBufferToFile(header=header,footer=False)

            if(self.progress_bar): printProgressBarColor(i_real,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)

            i_real += 1 # If success, increase i_real -- this is a counter for the number of successful events

        # Buffer gets written to the file in the loop above whenever it's full, but after exiting the loop
        # we should flush it once more -- even if empty, since this is still where we may write the footer.
        header = (self.loop_number == 0) and (not self.header_status) # for writing HepMC header
        footer = n_fail == 0
        self.WriteEventBufferToFile(header=header,footer=footer)

        # Delete the buffer files, if relevant.
        if(self.buffername is not None):
            comm = ['rm', self.buffername]
            try: sub.check_call(comm,stderr=sub.DEVNULL)
            except: pass

        return i_real-1, n_fail # note that i_real is using 1-indexing, which is what HepMC events use

    def GenerateSingle(self,event_number=1):
        """
        Just calls generation a single time.
        """
        from pyHepMC3 import HepMC3 as hm # the HepMCSetup will have taken care of this -- so the package will be already cached
        self.pythia.Generate()
        hepmc_event = hm.GenEvent()
        self.hepmc_converter.fill_next_event1(self.pythia.GetPythia(),hepmc_event,event_number)
        return hepmc_event

    # Generate a bunch of events in the given pT range,
    # and save them to a HepMC file.
    # We do perform event selection: Only certain particles are saved to the file to begin with.
    @profile_method('PythiaGenerator.Generate')
    def Generate(self,nevents:int):
        self.nevents = nevents # total number of events we request

        self.writer.InitializeWriter()

        if(self.progress_bar): printProgressBarColor(0,nevents, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # TODO: Eventually make this toggleable, or pick whichever method performs better.
        #       Right now, batch is only ~10% faster due to the HepMC3 conversion, but
        #       pythia.nextBatch() output has some issues with production vertices.
        #       See: https://gitlab.com/Pythia8/releases/-/issues/634
        use_batch = True

        if(not use_batch):
            # Loop in such a way as to guarantee that we get as many events as requested.
            # This logic is required as events could technically fail selections, e.g. not have the
            # requested truth particles (depends on requested truth particles & processes).
            self.nevents_success = 0
            n_fail = nevents
            self.loop_number = 0
            while(n_fail > 0):
                self.nevents_success, n_fail = self.GenerationLoop(
                    nevents-self.nevents_success,
                    i_real=self.nevents_success+1,
                    nevents_disp = nevents
                )
                self.loop_number += 1

        else:
            self.nevents_success = 0
            batch_size_default = self.buffer_size
            batch_size = batch_size_default
            self.loop_number = 0
            while(self.nevents_success < nevents):
                batch_size = np.minimum(batch_size_default, nevents - self.nevents_success)
                self.nevents_success, n_fail = self.GenerateBatch(batch_size,i_real = self.nevents_success+1,nevents_disp=nevents)
                self.loop_number += 1

        self.writer.Close()

        self.metadata_handler.AddCitations(self.GetCitations())
        return

    # This returns a list of all unique process codes encountered,
    # not a list of per-event process codes.
    def GetUniqueProcessCodes(self):
        return self.pythia.GetProcessCodes()

    # Get a dictionary containing cross sections (and their uncertainties)
    # for every process that was run, organized by process code.
    # Note that turning on a single Pythia process flag can in principle
    # turn on multiple processes, e.g. HardQCD will provide for many different
    # processes and the codes will distinguish between them.
    # Cross-sections are given in mb.
    # TODO: The estimates will possibly be inaccurate, as we throw out certain
    # events that Pythia8 has generated when they fail our particle selectors, but
    # those thrown out events are still included in computing the cross-section estimates.
    def GetSigmaDictionary(self):
        return self.pythia.GetSigmaDictionary()