#
# This file contains various functions for handling HepMC objects,
# as well as converting Pythia8 event listings -- as accessed by our
# custom `PythiaWrapper` class -- into HepMC events.
# There is currently some redundancy, as we historically used the pyhepmc
# library, whereas we are now leveraging the official HepMC3 python bindings.
# For the time being, functions supporting the use of either package are available,
# in particular in case these are useful elsewhere.

import sys
import numpy as np
import ROOT as rt
import pyhepmc as hep # the "unofficial" Pythonic HepMC3 bindings -- useful for some things
# from pyHepMC3 import HepMC3 as hm # the official Pythonic HepMC3 bindings
import subprocess as sub
import pathlib
from typing import Union, Optional, List, TYPE_CHECKING
from util.hepmc.setup import HepMCSetup
from util.hepmc.Pythia8ToHepMC3 import Pythia8ToHepMC3

if TYPE_CHECKING: # Only imported during type checking -- avoids risk of circular imports
    from util.pythia.utils import PythiaWrapper
    from util.particle_selection.particle_selection import BaseSelector

    # make sure Python HepMC3 bindings are setup
    setup = HepMCSetup(verbose=True)
    setup.PrepHepMC() # will download/install if necessary
    sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm

class Pythia8HepMC3Writer:
    def __init__(self,filename:Optional[str]=None):

        self.setup = HepMCSetup(verbose=True)
        self.setup.PrepHepMC()
        if(self.setup.GetPythonDirectory() not in sys.path):
            sys.path = [self.setup.GetPythonDirectory()] + sys.path
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
        from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import WriterRootTree
        self.output = WriterRootTree(self.filename)

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

def _extract_particle_momenta(pythia_wrapper:'PythiaWrapper', particle_indices:Optional[np.ndarray]=None):
    momentum = pythia_wrapper.GetPxPyPzE(particle_indices)
    pdgid = pythia_wrapper.GetPdgId(particle_indices)
    status = pythia_wrapper.GetStatus(particle_indices, hepmc=True)
    return momentum, pdgid, status

def WritePyHepMCEventsAscii(filename:str, hepev_list:Union[list,hep.GenEvent]):
    if(type(hepev_list) is not list): hepev_list = [hepev_list]
    with hep.io.WriterAscii(filename) as f:
        for hepev in hepev_list:
            f.write_event(hepev)
    return

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

def PyHepMCOutputAscii(hepev_list:Union[list,hep.GenEvent],buffername:str,filename:str,header:bool=False,footer:bool=False):
    WritePyHepMCEventsAscii(buffername,hepev_list) # write the given event(s) to a buffer file
    CopyHepMCBufferToFile(buffername,filename,header,footer) # copy buffer file into the full file
    return

def CompressHepMC(files:list, delete:bool=True, cwd:Optional[str]=None):
    for file in files:
        compress_file = file.replace('.hepmc','.tar.bz2')
        if(cwd is not None): compress_file = '{}/{}'.format(cwd,compress_file)
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],file.split('/')[-1]]
        # if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)
        if(delete):
            sub.check_call(['rm',file.split('/')[-1]],shell=False,cwd=cwd)
    return

def ExtractPyHepMCEventsAscii(files:list,get_nevents:bool=False, silent:bool=False):
    events = []
    nevents = 0
    for file in files:
        # Check that the file exists.
        if(not pathlib.Path(file).exists()):
            if(not silent):
                print('Warning: Tried to access file {} but it does not exist!'.format(file))
            continue

        with hep.io.ReaderAscii(file) as f:
            for i,evt in enumerate(f):
                events.append(evt)
                if(get_nevents): nevents += 1

    if(get_nevents): return events, nevents
    return events

def PythiaWrapperToPyHepMC(pythia_wrapper:'PythiaWrapper', event_number:int) -> hep.GenEvent:
    """
    Together with the Pythia8 wrapper (another part of this package),
    this basically gives a (Pythonic) Pythia8->HepMC3 (ascii) interface.

    NOTE: This seems to break when ISR/FSR is turned on, so instead one
    should use the included Pythia8TOHepMC3 class, which is from the official
    HepMC3 repository.
    """
    hepev = hep.GenEvent()
    hepev.event_number = event_number
    momentum,pdgid,status = _extract_particle_momenta(pythia_wrapper) # fetch all the momenta
    mother_indices = pythia_wrapper.GetMothers()
    prod = pythia_wrapper.GetProd() # production vertices' positions (t,x,y,z)

    # Fetch the particles.
    particles = [hep.GenParticle(momentum=momentum[i], pid=pdgid[i], status=status[i]) for i in range(len(momentum))]

    # Now add the vertices, based on mother/daughter information.
    # Note that we do not want to create duplicate vertices, which will
    # happen if we naively loop over all particles' production vertices.
    vertices = []
    used_mothers = None # will store tuples of (vertex index, mother index)
    for i in range(len(momentum)):

        mother_idxs = mother_indices[i]
        nmothers = len(mother_idxs)
        if(nmothers == 0): continue

        new_vertex = True
        used_mother_idx = None
        if(used_mothers is not None):
            for j in range(nmothers):
                if(mother_idxs[j] in used_mothers[:,-1]):
                    new_vertex = False
                    used_mother_idx = mother_idxs[j]
                    break

        if(new_vertex):
            vertex = hep.GenVertex(np.roll(prod[i],-1)) # converting from (t,x,y,z) to (x,y,z,t) for constructor
            # add stuff
            vertex.add_particle_out(particles[i])
            for idx in mother_idxs:
                vertex.add_particle_in(particles[idx])
                mother_tuple = np.array([len(vertices),idx],dtype=int)
                if(used_mothers is None): used_mothers = np.array([mother_tuple])
                else: used_mothers = np.append(used_mothers,[mother_tuple],axis=0)
            vertices.append(vertex)
        else:
            # find the existing vertex
            vertex_idx = used_mothers[np.where(used_mothers[:,-1] == used_mother_idx)[0][0],0]
            vertices[vertex_idx].add_particle_out(particles[i])

        hepev.add_vertex(vertex)
    return hepev

def ExtractPyHepMCParticles(events:List[hep.GenEvent], nparticles_max:Optional[int]=None,selection:Optional['BaseSelector']=None):
    if(selection is not None):
        particles = [
            [ev.particles[x] for x in np.atleast_1d(selection(ev))] # TODO: use intertools.compress() -- tried before, but it didn't work as expected?
            for ev in events
        ]
        if(nparticles_max is not None):
            particles = [x[:int(np.minimum(len(x),nparticles_max))] for x in particles]

    else:
        particles = [
            ev.particles[:int(np.minimum(len(ev.particles),nparticles_max))]
            for ev in events
        ]
    return particles

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

def PythiaWrapperToHepMC(pythia_wrapper:'PythiaWrapper', event_number:int) -> 'hm.GenEvent':
    """
    Convert a Pythia8 event index to a HepMC3 event, using the official HepMC3
    Python bindings. This takes our custom "PythiaWrapper" object as an argument,
    though it really just interfaces with the underlying Pythia8 generator object.
    """
    from pyHepMC3 import HepMC3 as hm
    hepev = hm.GenEvent()
    converter = Pythia8ToHepMC3()
    converter.fill_next_event1(pythia_wrapper.GetPythia(),hepev,event_number)
    return hepev

def ExtractHepMCEvents(files:list,get_nevents:bool=False, silent:bool=False):
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

        input = hm.ReaderAscii(file)
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
    from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import ReaderRootTree
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

def ExtractHepMCParticles(events:List['hm.GenEvent'], nparticles_max:Optional[int]=None,selection:Optional['BaseSelector']=None):
    if(selection is not None):

        particles = [
            [ev.particles()[x] for x in np.atleast_1d(selection(ev))]# if x>0] # TODO: use intertools.compress() -- tried before, but it didn't work as expected?
            for ev in events
        ]
        if(nparticles_max is not None):
            particles = [x[:int(np.minimum(len(x),nparticles_max))] for x in particles]

    else:
        particles = [
            ev.particles()[:int(np.minimum(len(ev.particles()),nparticles_max))]
            for ev in events
        ]
    return particles

# ======================
# Some utility functions,
# Currently supporting both
# PyHepMC and HepMC
#=======================

def ParticleToVector(particle:Union[hep.GenParticle,'hm.GenParticle']):
    if(isinstance(particle,hep.GenParticle)):
        return rt.Math.PxPyPzEVector(particle.momentum.px, particle.momentum.py, particle.momentum.pz, particle.momentum.e)
    return rt.Math.PxPyPzEVector(particle.momentum().px(), particle.momentum().py(), particle.momentum().pz(), particle.momentum().e())

def GetParticleID(particle:Union[hep.GenParticle,'hm.GenParticle']):
    if(isinstance(particle,hep.GenParticle)):
        return particle.pid
    return particle.pid()

