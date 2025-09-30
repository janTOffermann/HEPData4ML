#
# This file contains various functions for handling HepMC objects,
# as well as converting Pythia8 event listings -- as accessed by our
# custom `PythiaPythonWrapper` class -- into HepMC events.
# There is currently some redundancy, as we historically used the pyhepmc
# library, whereas we are now leveraging the official HepMC3 python bindings.
# For the time being, functions supporting the use of either package are available,
# in particular in case these are useful elsewhere.

import numpy as np
import subprocess as sub
import pathlib
from typing import Union, Optional, List, TYPE_CHECKING
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from util.hepmc.readers import ReaderAscii, ReaderRootTree
from util.hepmc.Pythia8ToHepMC3 import PythiaToHepMC
from util.misc.timing import profile_method

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.pythia.pythia import PythiaPythonWrapper
    from util.particle_selection.particle_selection import BaseSelector

    # make sure the correct Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=False)
    python_dir = setup.GetPythonDirectory()
    prepend_to_pythonpath(python_dir)
    from pyHepMC3 import HepMC3 as hm

class Pythia8HepMC3Writer:
    def __init__(self,hepmc_dir:Optional[str]=None, filename:Optional[str]=None):
        self.setup = HepMCSetup(hepmc_dir,verbose=False)
        # self.setup.PrepHepMC()
        python_dir = self.setup.GetPythonDirectory()
        # uncache_hepmc3()
        prepend_to_pythonpath(python_dir)

        # from pyHepMC3 import HepMC3 as hm

        self.mode = None
        self.initialized = False
        self.output = None

        self.SetFilename(filename)

    def SetFilename(self,filename:Optional[str]=None):
        self.filename = filename
        if(filename is not None):
            if('.root' in filename):
                self.SetMode('root')
            else:
                self.SetMode('ascii')

    def GetFilename(self) -> str:
        return self.filename

    def SetMode(self,mode=str):
        self.mode = mode

    def GetMode(self)->str:
        return self.mode

    def InitializeWriter(self):
        assert self.mode is not None

        if(self.mode == 'root'):
            self._init_root()
        else:
            self._init_ascii()

        self.initialized = True

    def _init_root(self):
        import pyHepMC3.rootIO.pyHepMC3rootIO as rio
        # from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import WriterRootTree # TODO: This broke recently -- why? Something off with pyHepMC dependency handling. -Jan
        self.output = rio.HepMC3.WriterRootTree(self.filename)

    def _init_ascii(self):
        from pyHepMC3 import HepMC3 as hm
        self.output = hm.WriterAscii(self.filename)

    def Close(self):
        self.output.close()

    def Write(self,hepev_list:Union[list,'hm.GenEvent']):
        if(type(hepev_list) is not list): hepev_list = [hepev_list]
        for hepev in hepev_list:
            self.output.write_event(hepev)
        return

# Some utility functions -- partly for interfacing with our Pythia8 wrapper.

def _extract_particle_momenta(pythia_wrapper:'PythiaPythonWrapper', particle_indices:Optional[np.ndarray]=None):
    momentum = pythia_wrapper.GetPxPyPzE(particle_indices)
    pdgid = pythia_wrapper.GetPdgId(particle_indices)
    status = pythia_wrapper.GetStatus(particle_indices, hepmc=True)
    return momentum, pdgid, status

def CopyHepMCBufferToFile(buffername:str,filename:str,header:bool=False,footer:bool=False):
    """
    This file copies HepMC events from file 'buffername' to file 'filename'.
    This is a workaround for the lack of an "append" writing mode for pyhepmc.

    TODO: Might want to rework this, it is very clunky. Is this necessary when using
          the official HepMC3 Python bindings? Also might eventually consider non-plaintext
          files (HepMC3/ROOT interface?).
    """
    # For speed, we do this using Unix commands (though it seems a bit hacky).

    # We need to determine whether or not to include the HepMC header/footer,
    # which will always be in the buffer file.
    upper_trim = 0 if header else 3
    lower_trim = 0 if footer else 3 # 2 if using Unix head (?)
    comm1 = 'tail -n +{} {}'.format(upper_trim, buffername).split(' ')

    proc = sub.Popen(comm1,stdout=sub.PIPE,text=True)
    if(lower_trim == 0): proc_out = proc.stdout.read().split('\n')
    else: proc_out = proc.stdout.read().split('\n')[:-lower_trim]

    with open(filename,'a') as f:
        f.writelines(line + '\n' for line in proc_out)

def CompressHepMC(files:list, delete:bool=True, cwd:Optional[str]=None):
    for file in files:
        compress_file = file.replace('.hepmc','.tar.gz')
        if(cwd is not None): compress_file = '{}/{}'.format(cwd,compress_file)
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-czf',compress_file.split('/')[-1],file.split('/')[-1]]
        # if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)
        if(delete):
            sub.check_call(['rm',file.split('/')[-1]],shell=False,cwd=cwd)
    return


#================================
# Functions below using real
# HepMC3 Python bindings.
#================================

def WriteHepMCEventsAscii(filename:str, hepev_list:Union[list,'hm.GenEvent']):
    if(type(hepev_list) is not list): hepev_list = [hepev_list]
    from pyHepMC3 import HepMC3 as hm
    output = hm.WriterAscii(filename)
    for hepev in hepev_list:
        output.write_event(hepev)
    output.close()
    return

def HepMCOutputAscii(hepev_list:Union[list,'hm.GenEvent'],buffername:str,filename:str,header:bool=False,footer:bool=False):
    WriteHepMCEventsAscii(buffername,hepev_list) # write the given event(s) to a buffer file
    CopyHepMCBufferToFile(buffername,filename,header,footer) # copy buffer file into the full file
    return

def PythiaPythonWrapperToHepMC(pythia_wrapper:'PythiaPythonWrapper', event_number:int) -> 'hm.GenEvent':
    """
    Convert a Pythia8 event index to a HepMC3 event, using the official HepMC3
    Python bindings. This takes our custom "PythiaPythonWrapper" object as an argument,
    though it really just interfaces with the underlying Pythia8 generator object.
    """
    from pyHepMC3 import HepMC3 as hm
    hepev = hm.GenEvent()
    converter = PythiaToHepMC()
    converter.fill_next_event1(pythia_wrapper.GetPythia(),hepev,event_number)
    return hepev

@profile_method('ExtractHepMCEvents')
def ExtractHepMCEvents(files:Union[str,List[str]],get_nevents:bool=False, silent:bool=False):
    if(isinstance(files,str)):
        files = [files]
    events = []
    nevents = 0
    for file in files:
        events_tmp = []
        nevents_tmp = 0
        if(file.split('.')[-1].lower() == 'root'):
            events_tmp, nevents_tmp = ExtractHepMCEventsROOT(file,True,silent)
        else:
            events_tmp, nevents_tmp = ExtractHepMCEventsAscii(file,True,silent)
        events += events_tmp
        nevents += nevents_tmp

    if(get_nevents): return events, nevents
    return events

def ExtractHepMCEventsAscii(files:Union[list,str],get_nevents:bool=False, silent:bool=False):
    from pyHepMC3 import HepMC3 as hm
    events = []
    nevents = 0
    if(isinstance(files,str)): files = [files]
    for file in files:
        # Check that the file exists.
        if(not pathlib.Path(file).exists()):
            if(not silent):
                print('Warning: Tried to access file {} but it does not exist!'.format(file))
            continue

        input = ReaderAscii(file)
        while(True):
            evt = hm.GenEvent()
            input.read_event(evt)
            if(input.failed()):
                break
            events.append(evt)
            if(get_nevents): nevents += 1
        input.close()

    if(get_nevents): return events, nevents
    return events

def ExtractHepMCEventsROOT(files:Union[list,str],get_nevents:bool=False, silent:bool=False):
    from pyHepMC3 import HepMC3 as hm
    events = []
    nevents = 0
    if(isinstance(files,str)): files = [files]
    for file in files:
        # Check that the file exists.
        if(not pathlib.Path(file).exists()):
            if(not silent):
                print('Warning: Tried to access file {} but it does not exist!'.format(file))
            continue

        input = ReaderRootTree(file)
        while(True):
            evt = hm.GenEvent()
            input.read_event(evt)
            if(input.failed()):
                break
            events.append(evt)
            if(get_nevents): nevents += 1
        input.close()

    if(get_nevents): return events, nevents
    return events

@profile_method('ExtractHepMCParticles')
def ExtractHepMCParticles(events: List['hm.GenEvent'], nparticles_max: Optional[int] = None, selection: Optional['BaseSelector'] = None):
    if selection is not None:
        particles = []
        for ev in events:
            ev_particles = ev.particles()
            indices = selection(ev)

            # Truncate indices first - avoid accessing unnecessary particles
            if nparticles_max is not None and len(indices) > nparticles_max:
                indices = indices[:nparticles_max]

            selected = [ev_particles[i] for i in indices]
            particles.append(selected)
    else:
        particles = [
            ev.particles()[:nparticles_max] if nparticles_max else ev.particles()
            for ev in events
        ]
    return particles

def ParticleToMomenta(particle:'hm.GenParticle'):
    """
    Produce 4-momenta in both Cartesian and cylindrical (pt,eta,phi,m) bases,
    stacked as one array.
    """
    momentum = particle.momentum()
    vectors = np.array(
        [
            [momentum.e(), momentum.px(), momentum.py(), momentum.pz()],
            [momentum.pt(), momentum.eta(), momentum.phi(), momentum.m()]
        ]
    )
    return vectors


def ParticleToProductionVertex(particle:'hm.GenParticle'):
    """
    A little more complex than handling production vertices -- particles may not have end vertices, if they are stable!
    """
    prod_vertex_position = particle.production_vertex().position()
    return np.array([prod_vertex_position.t(), prod_vertex_position.x(), prod_vertex_position.y(), prod_vertex_position.z()])

def ParticleToEndVertex(particle:'hm.GenParticle'):
    """
    A little more complex than handling production vertices -- particles may not have end vertices, if they are stable!
    """
    end_vertex = particle.end_vertex()
    if(end_vertex is None):
        return np.full(4,np.nan)

    end_vertex_position = end_vertex.position()
    return np.array([end_vertex_position.t(), end_vertex_position.x(), end_vertex_position.y(), end_vertex_position.z()])
