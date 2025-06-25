import os
import numpy as np
import subprocess as sub
from util.pythia.utils import PythiaWrapper
from util.hepmc.hepmc import Pythia8HepMC3Writer
from util.hepmc.Pythia8ToHepMC3 import Pythia8ToHepMC3
from util.hepmc.setup import HepMCSetup, uncache_hepmc3, prepend_to_pythonpath

from util.qol_utils.progress_bar import printProgressBarColor
from typing import Optional,TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports, limits unnecessary imports
    from util.config import Configurator

class PythiaGenerator:
    """
    Generates events using Pythia8.
    """
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

        self.hepmc_converter = Pythia8ToHepMC3(self.configurator.GetHepMC3Directory())

        # Event filters. # TODO: May remove
        self.event_filter = None

        # Things for the progress bar.
        self.prefix = 'Generating events for pT bin [{},{}]:'.format(self.pt_min,self.pt_max)
        self.suffix = 'Complete'
        self.bl = 50

        self.SetOutputDirectory()

        self.writer = Pythia8HepMC3Writer()

        self.SetFilename('events.hepmc')
        self.SetEventFilterFlagFilename() # will give a default name
        self.filename_fullpath = None

        # self.diagnostic_plots = True
        # self.InitializeHistograms()

        # Containers for event-level information.
        self.weights = None
        self.process_codes = None
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

        # Variables related to handling of pileup.
        # TODO: Might want to eventually reorganize this, or move it into a separate class?
        # self.pileup_handler = self.configurator.GetPileupHandling()

        # self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

    def SetVerbose(self,val:bool):
        self.verbose = val
        self.pythia.SetVerbose(val)

    def SetEventFilter(self,filter):
        self.event_filter = filter
        if(self.event_filter is not None):
            self.event_filter.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import

    def SetEventFilterFlag(self,filter):
        self.event_filter_flag = filter
        if(self.event_filter_flag is not None):
            self.event_filter_flag.Initialize(self.configurator) # may be necessary for things like dynamic fastjet import
        else:
            self.event_filter_flag_filename = None

    def SetPythiaConfigFile(self,file:Optional[str]=None):
        self.pythia_config_file = file
        if(file is None):
            self.pythia_config_file = self.configurator.GetPythiaConfigFile() # picked up from dictionary in config/config.py

    def SetDefaultFilenames(self):
        self.SetOutputDirectory(dir)
        self.SetFilename('events.hepmc')
        self.SetEventFilterFlagFilename()

    def SetProgressBar(self,flag:bool):
        self.progress_bar = flag

    def SetOutputDirectory(self,dir:Optional[str]=None):
        if(dir is None): dir = os.getcwd()
        self.outdir = dir

    def SetFilename(self,name:str, rename_extra_files:bool=True):
        # if('.hepmc' not in name):
        #     name = '{}.hepmc'.format(name)
        self.filename = name
        if(rename_extra_files):
            self.SetEventFilterFlagFilename(None)

        if(self.writer is not None):
            self.writer.SetFilename('{}/{}'.format(self.outdir,self.filename))
        return

    def GetFilename(self):
        return self.filename

    def GetHistFilename(self):
        return self.hist_filename

    def SetEventFilterFlagFilename(self,name:Optional[str]=None):
        if(name is None): self.event_filter_flag_filename = self.filename.replace('.hepmc','_event_filter_flag.h5')
        else: self.event_filter_flag_filename = name
        return

    def GetEventFilterFlagFilename(self):
        return self.event_filter_flag_filename

    def ConfigPythia(self, config_file:Optional[str]=None, verbose:bool=False):
        """
        Prepare and apply Pythia configuration. This turns our settings (from our config file)
        into a list of strings ready to be input to Pythia8.
        """
        self.pythia_config = self.configurator.GetPythiaConfig(self.pt_min,self.pt_max, verbose)
        self.SetPythiaConfigFile(file=config_file)

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

    def FillEventBuffer(self,hepev_full):
        self.hepev_buffer.append(hepev_full)

    def WriteEventBufferToFile(self,header:bool=False,footer:bool=False):
        if(header): self.header_status = True
        if(footer): self.footer_status = True

        self.writer.Write(self.hepev_buffer)

        # PyHepMCOutput(self.hepev_buffer,self.buffername,self.filename_fullpath,header,footer)
        # HepMCOutputAscii(self.hepev_buffer,self.buffername,self.filename_fullpath,header,footer)

        self.ClearEventBuffer()

    def GenerationLoop(self, nevents,i_real:int = 1, nevents_disp:Optional[int]=None):
        """
        This is the function where Pythia8 matrix element generation + showering/hadronization happens.
        The results are filtered for the selected "truth" and "final state" particles, which are placed
        into HepMC3 events. These events are periodically written to a buffer file (which is then merged
        into the "master" HepMC3 file).
        """
        from pyHepMC3 import HepMC3 as hm # the HepMCSetup will have taken care of this -- so the package will be already cached
        n_fail = 0
        if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

        # The way that pyhepmc's WriterAscii works, writing an event will overwrite the whole file.
        # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
        # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
        # down things, so we ultimately want to find some way to do a write with "append" functionality.

        self.filename_fullpath = '{}/{}'.format(self.outdir,self.filename)

        # For ASCII mode, create buffer file.
        if(self.writer.GetMode().lower() == 'ascii'):
            self.buffername = self.filename_fullpath.replace('.hepmc','_buffer.hepmc')

        for i in range(nevents):

            self.pythia.Generate() # generate an event!

            # # Here, we optionally apply pileup (or more generally, merging in some other event(s)).
            # if(self.pileup_handler is not None):
            #     #TODO: This needs to be reworked.
            #     self.pileup_handler(self.pythia)

            # ==========================================
            # Now we apply an (optional) "event filter". This applies some condition to the set
            # of truth and final-state particles selected above, and requires that the event pass
            # this condition. If not we will count the event as failed, and generate another.
            # ==========================================
            if(self.event_filter is not None):
                passed_filter = self.event_filter(self.pythia)
                if(not passed_filter):
                    n_fail += 1
                    continue

            # # ==========================================
            # # Now we apply an (optional) "event filter" flag. This applies some condition to the set
            # # of truth and final-state particles selected above, and checks if the event passes
            # # this condition. The result is recorded in an HDF5 file that will be merged into the
            # # final dataset.
            # # ==========================================
            # # TODO: Rework or remove this. Don't want to create anything other than the HepMC at this stage.
            # if(self.event_filter_flag is not None):
            #     event_filter_flag_filename_full = '{}/{}'.format(self.outdir,self.event_filter_flag_filename)
            #     f = h5.File(event_filter_flag_filename_full,'a')
            #     key = self.event_filter_flag.GetName()
            #     filter_flag = self.event_filter_flag(self.pythia)
            #     if(i_real == 1):
            #         filter_flag = np.expand_dims(filter_flag,axis=0)
            #         f.create_dataset(key,data=filter_flag,compression='gzip',chunks=True,maxshape=(None,))
            #     else:
            #         f[key].resize(f[key].shape[0] + 1,axis=0)
            #         f[key][-1] = filter_flag
            #     f.close()

            # ==========================================
            # Now lets create the HepMC event.
            # ==========================================
            # hepmc_event = PythiaWrapperToPyHepMC(self.pythia,i_real)

            hepmc_event = hm.GenEvent()
            self.hepmc_converter.fill_next_event1(self.pythia.GetPythia(),hepmc_event,i_real)
            # hepmc_event = PythiaWrapperToHepMC(self.pythia,i_real)

            # Fill the memory buffer with this event.
            self.FillEventBuffer(hepmc_event)

            # If buffer is full (i.e. has reached max size), write it to file & flush.
            if(self.GetCurrentBufferSize() == self.buffer_size):
                header = (self.loop_number == 0) and (not self.header_status)
                self.WriteEventBufferToFile(header=header,footer=False)


            if(self.progress_bar): printProgressBarColor(i_real,nevents_disp, prefix=self.prefix, suffix=self.suffix, length=self.bl)

            # Record this event's weight and process code from Pythia.
            #TODO: Consider removing this, or adjusting how it is handled -- if the user provides
            #      externally-produced HepMC files, they won't have these things.
            self.weights[i_real-1] = self.pythia.GetEventWeight() # convert from 1-indexing to 0-indxing
            self.process_codes[i_real-1] = self.pythia.GetProcessCode() # convert from 1-indexing to 0-indxing

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

    def GenerateSingle(self):
        """
        Just calls generation a single time.
        To be used for testing.
        """
        self.pythia.Generate()
        return

    # Generate a bunch of events in the given pT range,
    # and save them to a HepMC file.
    # We do perform event selection: Only certain particles are saved to the file to begin with.
    def Generate(self,nevents:int):
        self.nevents = nevents # total number of events we request

        self.writer.InitializeWriter()

        if(self.progress_bar): printProgressBarColor(0,nevents, prefix=self.prefix, suffix=self.suffix, length=self.bl)

        # Loop in such a way as to guarantee that we get as many events as requested.
        # This logic is required as events could technically fail selections, e.g. not have the
        # requested truth particles (depends on requested truth particles & processes).
        self.nevents_success = 0
        n_fail = nevents
        self.loop_number = 0
        self.weights       = np.zeros(nevents)
        self.process_codes = np.zeros(nevents, dtype=int)
        while(n_fail > 0):
            self.nevents_success, n_fail = self.GenerationLoop(
                nevents-self.nevents_success,
                i_real=self.nevents_success+1,
                nevents_disp = nevents
            )
            self.loop_number += 1

        self.writer.Close()

        # # Fill out an array with each event's cross section (and the uncertainty on that cross-section).
        # self.xsecs = np.zeros((nevents,2))
        # xsec_dictionary = self.GetSigmaDictionary()
        # for i,pcode in enumerate(self.process_codes):
        #     self.xsecs[i,:] = xsec_dictionary[pcode]

        # # Create/append a stats file, and put in information on event weights & cross-sections.
        # #TODO: Will want to rework this.
        # f = h5.File('{}/{}'.format(self.outdir,self.stats_filename),'a')
        # compression = 'gzip'
        # copts = 9
        # f.create_dataset('mc_weight',data=self.weights,compression=compression,compression_opts=copts)
        # f.create_dataset('process_code',data=self.process_codes,compression=compression,compression_opts=copts)
        # f.create_dataset('cross_section',data=self.xsecs[:,0],compression=compression,compression_opts=copts)
        # f.create_dataset('cross_section_uncertainty',data=self.xsecs[:,1],compression=compression,compression_opts=copts)
        # f.close()

        return

    # Get the MC event weights (will typically just be 1 for each event).
    # This reads from a file, which is only written to when an event is written to disk.
    # Thus these weights will line up with the events we've saved, i.e. we don't have to worry
    # about events that Pythia8 successfully generated but which we threw out because they failed
    # one of our selectors.
    def GetEventWeights(self):
        return self.weights

    def GetProcessCodes(self):
        return self.process_codes

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