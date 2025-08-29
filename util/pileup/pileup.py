import ROOT as rt
import uproot as ur
import numpy as np
import h5py as h5
import glob,sys,os,pathlib,itertools
import subprocess as sub
from util.qol_utils.progress_bar import printProgressBarColor
from util.qol_utils.pdg import DatabasePDG
from util.meta import MetaDataHandler
from util.config import Configurator
from util.math.rotations import RotateVector
from util.post_processing.jets import TruthJetFinder
from util.hepmc.setup import HepMCSetup, prepend_to_pythonpath
from util.hepmc.readers import ReaderAscii, ReaderRootTree # our wrappers for the HepMC3 reader classes
from util.hepmc.hepmc import ParticleToVector
from util.generation import PythiaGenerator
from typing import List, Union, Tuple, Optional, TYPE_CHECKING

if(TYPE_CHECKING):
    setup = HepMCSetup(verbose=False)
    setup.PrepHepMC() # will download/install if necessary
    python_dir = setup.GetPythonDirectory()
    if(python_dir not in sys.path):
        sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
    from pyHepMC3 import HepMC3 as hm

class PileupOverlay:
    """
    This class performs on-the-fly mixing, to mix in events from some pileup HepMC3 file(s)
    into "main events" (from some other HepMC file).
    It creates mixed events, which are output to a new file.

    In practice, you might want to use this to create pileup events for "pre-mixing",
    like what is done in Monte Carlo sample preparation for the CMS experiment.
    """

    def __init__(self, pileup_files:Optional[Union[str,list]]=None,rng_seed:int=1,mu_file:str=None):

        # Upon start, make sure that hepmc is set up
        self._init_hepmc()

        self.verbosity = 1
        self.print_prefix = 'PileupOverlay: '

        self.require_pileup_input = True
        self.files = None
        self.SetPileupFiles(pileup_files)

        self.selected_indices = None # transient storage for indices selected for a particular event
        self.allow_reuse = True
        self.pileup_indices = {} # accumulate the pileup event indices per event -> for output. Dictionary keys are filename

        self.SetRNGSeed(rng_seed)
        self.beam_spot_sigma = None # will store the beam spot size
        self.SetBeamSpotSigma() # sets some sensible defaults
        self.do_phi_rotations = True
        self.phi_rotations_transient = []
        self.phi_rotations = {} # accumulate the phi rotations per event -> for output. Dictionary keys are filename

        # Set up the distribution for # of interactions per crossing
        self.mu_file = mu_file
        self.available_mu_values = None
        self.mu_probabilities = None
        self._init_mu_distribution()
        self.mu = None # transient storage for the current value of mu
        self.mu_values = {} # will store a list of all used mu values

        self.pileup_events = None # transient storage for pileup events

        self.event_buffer_size = 10 # how many combined events to store in memory, before flushing to output file
        #TODO: Need to fix buffer behavior -- the HepMC writers overwrite the output_file each time they're initialized, as opposed to appending

        self.indir = None
        self.outdir = None

        # for particle charge lookup
        self.pdg_database = DatabasePDG()

        self.metadata_handler = None
        self.configurator = None

    def SetMetadataHandler(self,handler:'MetaDataHandler'):
        self.metadata_handler = handler

    def SetConfigurator(self,configurator:'Configurator'):
        self.configurator = configurator

    def GetRNGSeed(self):
        return self.rng_seed

    def SetRNGSeed(self,rng_seed:int):
        if(rng_seed < 0):
            rng_seed = 0
        self.rng_seed = rng_seed
        self.rng = np.random.default_rng(self.rng_seed)

    def SetVerbosity(self,val:int):
        self.verbosity = val

    def SetInputDirectory(self,val:str):
        self.indir = val

    def SetOutputDirectory(self,val:str):
        self.outdir = val

    def SetAllowReuse(self, val:bool):
        self.allow_reuse = val

    def SetUsePhiRotations(self, val:bool):
        self.do_phi_rotations = val

    def _init_hepmc(self):
        setup = HepMCSetup(verbose=False)
        python_dir = setup.GetPythonDirectory()
        prepend_to_pythonpath(python_dir)
        # now can do "from pyHepMC3 import HepMC3 as hm"

    def SetPileupFiles(self,files:Union[str,list]=None):
        if(files is None):
            return
        if(isinstance(files,str)):
            self.files = glob.glob(files)
            if(len(self.files) == 0): # maybe this was not an absolute path, but a relative one (w.r.t. cwd)
                self.files = glob.glob('{}/{}'.format(os.getcwd(),files))
        else:
            self.files = files
            for i,file in enumerate(self.files):
                if(not pathlib.Path(file).exists()):
                    self.files[i] = '{}/{}'.format(os.getcwd(),file)
        self.files = sorted([str(pathlib.Path(x).absolute()) for x in self.files if pathlib.Path(x).exists()])

        # determine the total number of events in these files
        self.n_dict = {}
        self.first_idx = {}
        self.n_total = 0
        self._index_files()

        self.event_indices = np.arange(self.n_total)
        self.event_mask = np.full(self.n_total,True,dtype=bool)
        return

    def _index_ascii(self,file):
        with open(file,'r') as f:
            n =  len([x for x in f.readlines() if 'E ' in x])
            self.first_idx[file] = self.n_total
            self.n_dict[file] = n
            self.n_total += n
        return

    def _index_root(self,file):
        f = rt.TFile(file,"READ")
        t = f.Get("hepmc3_tree")
        n = t.GetEntries()
        f.Close()
        self.first_idx[file] = self.n_total
        self.n_dict[file] = n
        self.n_total += n
        return

    def _index_files(self):
        for file in self.files:

            # indexing for ascii
            if(file.split('.')[-1].lower() != 'root'):
                self._index_ascii(file)

            # indexing for ROOT
            else:
                self._index_root(file)

    def SetBeamSpotSigma(self,dt:float=0.16,dx:float=0.01,dy:float=0.01,dz:float=35.):
        """
        Sets the size of the beamspot (standard deviation) in (t,x,y,z), in nanoseconds or millimeters (as appropriate).
        The default spatial values are based on ATLAS Run 2: https://twiki.cern.ch/twiki/bin/view/AtlasPublic/BeamSpotPublicResults#Run_2_25ns_pp_Collisions_s_13_Te .
        The default time value is based on the Delphes ATLAS detector card (with pileup).
        """
        self.beam_spot_sigma = (dt,dx,dy,dz)

    def _init_mu_distribution(self):
        hist_name = 'SimplePileup_mu'

        if(self.mu_file is None):
            # default
            nbins = 80
            self.mu_distribution = rt.TH1F(hist_name,'',nbins,0,nbins)
            # very approximate for Run 2,
            # see https://atlas.web.cern.ch/Atlas/GROUPS/DATAPREPARATION/PublicPlots/2018/DataSummary/figs/mu_2015_2018.png
            mu = 33.7
            sigma = 11.5
            for i in range(nbins):
                self.mu_distribution.SetBinContent(i+1,self._gaussian(i,mu,sigma))

        else:
            try:
                f = rt.TFile(self.mu_file,"READ")
                self.mu_distribution = f.Get(hist_name).Clone()
                self.mu_distribution.SetDirectory(0)
                f.Close()
            except:
                print("Error: SimplePileup failed to read histogram {} from file {}.".format(hist_name,self.mu_file))
                print("Falling back on default mu distribution.")
                self.mu_file = None
                self._init_mu_distribution()

        n = self.mu_distribution.GetNbinsX()
        mu_min = self.mu_distribution.GetXaxis().GetBinLowEdge(1)
        mu_max = self.mu_distribution.GetXaxis().GetBinLowEdge(n+1)
        self.available_mu_values = np.arange(mu_min,mu_max)
        self.mu_probabilities = np.array([self.mu_distribution.GetBinContent(i+1) for i in range(n)],dtype=float)
        self.mu_probabilities /= np.sum(self.mu_probabilities)

    def _sample_mu_distribution(self):
        self.mu = int(self.rng.choice(self.available_mu_values,p=self.mu_probabilities))

    def _pick_event_indices(self):
        self._sample_mu_distribution()
        self.selected_indices = self.rng.choice(self.event_indices[self.event_mask],size=self.mu)
        self.selected_indices.sort() # not strictly necessary, but nice since we do some bookkeeping
        if(not self.allow_reuse):
            self._update_mask()

    def _update_mask(self,indices=None):
        """
        Remove "indices" from pileup event indices that can be sampled from in the future.
        Can be used to achieve sampling without replacement.
        """
        if(indices is None): indices = self.selected_indices # use currently selected indices
        self.event_mask[self.selected_indices] = False

    def _fetch_event_single_file_ascii(self,filename:str,indices:Union[int,list,np.ndarray]):
        # NOTE: With ASCII, we require iterating through the whole file
        # and that will scale quite poorly as these files get large.

        from pyHepMC3 import HepMC3 as hm

        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)
        indices = sorted(indices)

        events = [hm.GenEvent()] * len(indices)
        counter = 0
        write_counter = 0

        reader = ReaderAscii(filename)

        while(not reader.failed()):

            if(write_counter == len(indices)):
                break

            # immediately skip to the 1st event of interest -- hopefully give some speedup.
            if(counter == 0):
                reader.skip(indices[0])
                counter += indices[0]

            evt = hm.GenEvent()
            reader.read_event(evt)
            if(reader.failed): break

            if(counter == indices[0]): # maybe faster than "counter in indices"? Requires the roll below
                indices = np.roll(indices,-1)
                events[write_counter] = evt
                write_counter += 1
            counter += 1
        reader.close()
        return events

    def _fetch_event_single_file_root(self,filename:str,indices:Union[int,list,np.ndarray]):
        # NOTE: Leverages functions currently only available in our custom fork of HepMC3,
        #       for accessing non-sequential events from a ROOT file.
        from pyHepMC3 import HepMC3 as hm

        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)
        indices = sorted(indices)
        indices = np.array([x - self.first_idx[filename] for x in indices],dtype=int) # adjust from global indexing, to indexing iwthin this file

        events = [hm.GenEvent()] * len(indices)
        reader = ReaderRootTree(filename)

        for i,idx in enumerate(indices):
            evt = hm.GenEvent()
            status = reader.read_event_at_index(evt,idx)
            if(reader.failed() or (not status)):
                assert False
            events[i] = evt

        reader.close()
        return events

    def _fetch_event_single_file(self,filename:str,indices:Union[int,list,np.ndarray]):
        # NOTE: Will not work nicely with ASCII, since that requires iterating
        # through the whole file and that will scale quite poorly as these files get large.
        # NOTE: For HepMC3/ROOT, random file access is in principle supported, but it is not
        #       clear to me how to actually accomplish this with the existing ReaderRootTree
        #       methods. Maybe we need to write something custom, and make sure it interfaces
        #       nicely in Python? - Jan

        from pyHepMC3 import HepMC3 as hm
        # from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import ReaderRootTree

        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)
        indices = sorted(indices)

        if(filename.split('.')[-1].lower() == 'root'):
            return self._fetch_event_single_file_root(filename,indices)
        else:
            return self._fetch_event_single_file_ascii(filename,indices)

    def _fetch_events(self,indices:Union[int,list,np.ndarray]):
        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)

        # determine which files we'll have to pull pileup events from
        file_indices = {}
        for idx in indices:
            file_idx_1 = np.where(idx >= np.array(list(self.first_idx.values()),dtype=int))[0]
            file_idx_2 = np.where(idx < np.array(list(self.first_idx.values())) + np.array(list(self.n_dict.values()),dtype=int))[0]
            file_idx = np.intersect1d(file_idx_1,file_idx_2)[0]

            if(file_idx not in file_indices.keys()):
                file_indices[file_idx] = []

            file_indices[file_idx].append(idx)

        events = []
        for key,val in file_indices.items():
            events += self._fetch_event_single_file(self.files[key],val)
        return events

    def _get_nevents_root(self,filename:str):
        f = ur.open(filename)
        t = f['hepmc3_tree']
        n = t.num_entries
        f.close()
        return n

    def __call__(self,input_file:str,output_file:Optional[str]=None):
        from pyHepMC3 import HepMC3 as hm

        input_file_extension = input_file.split('.')[-1]

        if(self.indir is not None):
            input_file = '{}/{}'.format(self.indir,input_file)

        if(output_file is None):
            output_file = input_file.replace('.{}'.format(input_file_extension),'_mixed.{}'.format(input_file_extension))
        elif(self.outdir is not None):
            output_file = '{}/{}'.format(self.outdir,output_file)

        # buffer file is only actually relevant for the ASCII mode
        # TODO: Have to check how ROOT writing works -- appends to existing file, or erases it?
        buffer_file = output_file.replace('.{}'.format(input_file_extension),'_buffer.{}'.format(input_file_extension))

        if(((self.files is None) or len(self.files) == 0) and self.require_pileup_input):
            print('Error: No pileup files for PileupOverlay.')
            return

        events = []

        # technically allows toggling mode from file to file... though one really should stick to ROOT.
        mode = 'root'
        if(input_file_extension.lower() != 'root'):
            mode = 'ascii'

        nevents = None
        if(mode=='root'):
            nevents = self._get_nevents_root(input_file)

        if(mode=='root'):
            reader = ReaderRootTree(input_file)
        else:
            reader = ReaderAscii(input_file)

        i = 0
        while(not reader.failed()):

            evt = hm.GenEvent()
            reader.read_event(evt)
            if(reader.failed()):
                # if(self.verbosity > 0 and nevents is not None):
                #     printProgressBarColor(nevents,nevents,'Adding pileup:',decimals=2)
                break

            self._pick_event_indices() # sets self.mu
            pileup_events = self._fetch_events(self.selected_indices)

            # overlay pileup on this event
            # for now, use automatic displacement -- will use self.beam_spot_sigma
            evt = self._combine_event_with_pileup(evt,pileup_events)

            # store the combined event in memory
            events.append(evt)

            # record the selected indices, and phi rotations
            self._record_pileup_info(output_file)

            # flush events from memory as needed
            if(len(events) >= self.event_buffer_size):

                # write to output file.
                # TODO: Is append functionality supported by hepmc writers?
                self._flush_to_file(events,output_file)

            i += 1

            if(self.verbosity > 0 and nevents is not None):
                printProgressBarColor(i,nevents,'Adding pileup:',decimals=2)

        # one last flush for any stragglers
        if(len(events) > 0):
            self._flush_to_file(events,output_file)

        try:
            # delete the buffer file
            os.unlink(buffer_file)
        except:
            pass

        return output_file

    def _record_pileup_info(self,output_file):
        key = output_file.split('/')[-1] # remove any leading directory, just use filename (this is generally our convention)
        if(key not in self.pileup_indices.keys()):
            self.pileup_indices[key] = []
            self.phi_rotations[key] = []
            self.mu_values[key] = []
        self.pileup_indices[key].append(self.selected_indices)
        self.mu_values[key].append(self.mu)

        # only record phi rotations if they were used -- no point in recording lots of zeros
        if(self.do_phi_rotations):
            self.phi_rotations[key].append(self.phi_rotations_transient)
        return

    def _migrate_pileup_info(self,old_key,new_key):
        """
        A bit of a messy function, a consequence of some file naming & I/O options.
        Simply renames a key in the info dictionaries.
        """
        for d in [self.mu_values,self.pileup_indices,self.phi_rotations]:
            if(old_key in d.keys()):
                d[new_key] = d[old_key]
                del d[old_key]

    def Process(self,inputs,outputs=None):
        replace_inputs = False
        if(outputs is None):
            replace_inputs = True
            final_outputs = []
            outputs = [] # temporary filenames
            for input in inputs:
                file_extension = input.split('.')[-1]
                basename = '.'.join(input.split('.')[:-1])
                outputs.append('{}.pileup.{}'.format(basename,file_extension))

        for (input,output) in zip(inputs,outputs):
            self.__call__(input,output)

            if(replace_inputs):
                input_fullpath = input
                if(self.indir is not None):
                    input_fullpath = '{}/{}'.format(self.indir,input)
                output_old = output
                output_new = input
                if(self.outdir is not None):
                    output_old = '{}/{}'.format(self.outdir,output)
                    output_new = '{}/{}'.format(self.outdir,input)

                os.unlink(input_fullpath)
                command = ['mv',output_old,output_new]
                sub.check_call(command)

                # we also need to rename some keys in a few dictionaries that were based on output_old
                # TODO: This is a bit messy, is there a better way?
                if(output in self.mu_values.keys()):
                    self._migrate_pileup_info(output,input)
                final_outputs.append(input) # TODO: a little messy in terms of code

        self._writeMetadata()
        return

    def _flush_to_file(self,events:List['hm.GenEvent'],output_file:str,buffername:str=None):
        from pyHepMC3 import HepMC3 as hm
        from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import WriterRootTree

        if(output_file.split('.')[-1].lower() == 'root'):
            writer = WriterRootTree(output_file, True) # uses our custom HepMC3 functionality for "append mode"
        else:
            assert False # for now, we dont't support ASCII output since there isn't (yet) append functionality

        for evt in events:
            writer.write_event(evt)

        writer.close()
        events.clear()
        return

    def _gaussian(self,x,mu,sig,A=None):
        if(A is None): A = 1. / (np.sqrt(2.0 * np.pi))
        return A * np.exp(-np.square((x - mu) / sig) / 2)

    def _sumpt2(self,evt:'hm.GenEvent'):
        """
        Compute the sum of pt2 of charged particles in the event.
        Used for scaling the x/y displacement of the primary vertex,
        see https://github.com/delphes/delphes/blob/d256775e652525b0c35929e72a8bf20252328696/modules/PileUpMerger.cc#L181.
        """

        sumpt2 = 0.
        for i,particle in enumerate(evt.particles()):
            charge = self.pdg_database.GetCharge(particle.pid()) # charge is in units of |e|/3
            if(np.abs(charge) > 1.0e-9): # abs might not be needed based on the above
                sumpt2 += np.square(particle.momentum().px()) + np.square(particle.momentum().py())
        return sumpt2

    def _combine_event_with_pileup(self,
        main_event: 'hm.GenEvent',
        pileup_events: Union['hm.GenEvent', List['hm.GenEvent']],
        pileup_displacements: Optional[Union[Tuple[float, float, float, float],
                                        List[Tuple[float, float, float, float]]]] = None,
        auto_displace: bool = True,
        beam_spot_sigma: Optional[Tuple[float,float,float,float]] = None
    ) -> 'hm.GenEvent':
        """
        Combine a main event with pileup events, applying vertex displacements.

        Parameters:
        -----------
        main_event : hm.GenEvent
            The primary hard-scatter event
        pileup_events : hm.GenEvent or List[hm.GenEvent]
            Single pileup event or list of pileup events to overlay
        pileup_displacements : tuple or list of tuples, optional
            (t, x, y, z) displacements for each pileup event in mm/ns units
            If None and auto_displace=True, random displacements are generated
        auto_displace : bool
            If True, automatically generate random displacements when none provided
        beam_spot_sigma : tuple, optional
            Standard deviation for random displacement generation (mm/ns), in (t,x,y,z).
            Defaults to self.beam_spot_sigma (preferred use).

        Returns:
        --------
        hm.GenEvent
            Combined event with displaced pileup vertices
        """

        from pyHepMC3 import HepMC3 as hm

        # Ensure pileup_events is a list
        if isinstance(pileup_events, hm.GenEvent):
            pileup_events = [pileup_events]

        # Generate or validate displacements
        if pileup_displacements is None:
            pileup_displacements = np.zeros((len(pileup_events),4))
            if auto_displace:
                if(beam_spot_sigma is None):
                    beam_spot_sigma = self.beam_spot_sigma
                beam_spot_sigma = list(beam_spot_sigma)
                beam_spot_sigma[0] *= 1.0e6 / rt.TMath.C() # convert from ns to mm/c
                pileup_displacements = np.random.normal(0,beam_spot_sigma,(len(pileup_events),4))

        if len(pileup_displacements) != len(pileup_events):
            raise ValueError("Number of displacements ({}) must match number of pileup events ({})".format(len(pileup_displacements),len(pileup_events)))

        # We displace the main event as well.
        # NOTE: Based on the approach taken in Delphes (https://github.com/delphes/delphes/blob/d256775e652525b0c35929e72a8bf20252328696/modules/PileUpMerger.cc#L181),
        # it looks like we should scale the x/y displacement by the sum of pT^2 of the event.
        # Note that there, the vertex position is basically some sort of pt2-weighted average of attached particles' "positions",
        # though it's not a perfect average since all particles contribute to the numerator (regardless of charge) but only
        # charged particles contribute to the denominator. I don't know how well-motivated this really is. - Jan
        combined_event = main_event

        main_event_displacement = np.random.normal(0,beam_spot_sigma,4)
        sumpt2 = self._sumpt2(main_event)
        main_event_displacement[1:3] /= sumpt2 # adjust the x- and y-displacements, to make them smaller based on sum of pt2 of charged particles. A little unclear on units/scale here...
        combined_event.shift_position_by(hm.FourVector(*np.roll(main_event_displacement,-1))) # using np.roll to get from (t,x,y,z) to (x,y,z,t)
        # self._displace_event(combined_event,*main_event_displacement)

        # Add each pileup event with displacement
        self.phi_rotations_transient.clear()
        for pileup_event, (dt, dx, dy, dz) in zip(pileup_events, pileup_displacements):
            self._add_displaced_event(combined_event, pileup_event, dt, dx, dy, dz)

        return combined_event

    def _add_displaced_event(self,
                            target_event: 'hm.GenEvent',
                            pileup_event: 'hm.GenEvent',
                            dt: float, dx: float, dy: float, dz: float):
        """Add a pileup event to the target event with vertex displacement."""
        from pyHepMC3 import HepMC3 as hm

        # TODO: Consider using hm.GenEvent.shift_position_by() here? Maybe not possible if doing a phi rotation.

        # Add a random rotation in phi to the input pileup event.
        # (rotate, then translate)
        phi_rotation_angle = self.rng.uniform(0.,2 * np.pi)
        if(self.do_phi_rotations):
            self.phi_rotations_transient.append(phi_rotation_angle)

        # Create mapping from old vertex objects to new vertex objects using list index
        vertex_map = {}
        vertices_list = list(pileup_event.vertices())

        # Copy vertices with displacement
        for i, vertex in enumerate(vertices_list):
            # Apply displacement to vertex position

            old_position = np.array([getattr(vertex.position(),method)() for method in ['t','x','y','z']])
            if(self.do_phi_rotations):
                new_position = RotateVector(old_position,phi_rotation_angle,0.,0.)
            else:
                new_position = old_position
            new_position += np.array([dt,dx,dy,dz])
            new_position = hm.FourVector(*np.roll(new_position,-1)) # using np.roll to get from (t,x,y,z) to (x,y,z,t) for the constructor

            new_vertex = hm.GenVertex(new_position)
            new_vertex.set_status(vertex.status())
            # target_event.add_vertex(new_vertex)
            vertex_map[i] = new_vertex

        # Copy particles and establish relationships
        for particle in pileup_event.particles():

            old_momentum = np.array([getattr(particle.momentum(),method)() for method in ['e','px','py','pz']])
            if(self.do_phi_rotations):
                new_momentum = hm.FourVector(*np.roll(RotateVector(old_momentum,phi_rotation_angle,0.,0.),-1))
            else:
                new_momentum = hm.FourVector(*np.roll(old_momentum,-1)) # using np.roll to get from (e,px,py,pz) to (px,py,pz,e) for the constructor

            new_particle = hm.GenParticle(
                new_momentum,
                particle.pid(),
                particle.status()
            )
            new_particle.set_generated_mass(particle.generated_mass())

            # Set production vertex if it exists
            if particle.production_vertex():
                try:
                    vertex_idx = vertices_list.index(particle.production_vertex())
                    prod_vertex = vertex_map[vertex_idx]
                    prod_vertex.add_particle_out(new_particle)
                except ValueError:
                    pass # expected to be triggered, by beam particles

            # Set end vertex if it exists
            if particle.end_vertex():
                try:
                    vertex_idx = vertices_list.index(particle.end_vertex())
                    end_vertex = vertex_map[vertex_idx]
                    end_vertex.add_particle_in(new_particle)
                except ValueError:
                    print("Warning: Could not find end vertex for particle {}".format(particle.pid))

        for i,vtx in vertex_map.items():
            target_event.add_vertex(vtx)

    def _writeMetadata(self):
        pileup_metadata = {}
        if(self.metadata_handler is None):
            self._print('Unable to write metadata; no handler was provided.')
            return

        # Here, the metadata from the pileup files themselves will be fetched.
        if(self.files is not None):
            for pileup_file in self.files:
                pileup_file_extension = pileup_file.split('.')[-1]
                if(pileup_file_extension.lower() != 'root'):
                    self._print('Warning: Cannot read in metadata from pileup file {} .'.format(pileup_file))
                    continue
                key = pileup_file.split('/')[-1]
                pileup_metadata[key] = self.metadata_handler.ReadMetaDataFromROOTFile(pileup_file)
            self.metadata_handler.AddElement('Metadata.Pileup.InputMetadata',pileup_metadata)

        # Also add info on the mu distribution that was used
        mu_bins = np.array([self.mu_distribution.GetBinLowEdge(x+1) for x in range(self.mu_distribution.GetNbinsX())]) # left edges)
        mu_weights = np.array([self.mu_distribution.GetBinContent(x+1) for x in range(self.mu_distribution.GetNbinsX())])
        self.metadata_handler.AddElement('Metadata.Pileup.MuDistribution.BinEdgesLeft',mu_bins)
        self.metadata_handler.AddElement('Metadata.Pileup.MuDistribution.BinContents',mu_weights)

        return

    def GetMuValues(self):
        return self.mu_values

    def GetPileupInfo(self):
        data = {}
        data['Pileup.Mu'] = self.mu_values
        data['Pileup.Index'] = self.pileup_indices
        if(self.do_phi_rotations):
            data['Pileup.PhiRotation'] = self.phi_rotations
        return data

    def AddPileupInfoToH5(self,h5_file,file_key, cwd=None,copts=9):
        if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

        f = h5.File(h5_file,'r+')
        keys = list(f.keys())
        nevents = f[keys[0]].shape[0]
        data = self.GetPileupInfo()

        for key,value_dict in data.items():
            # for lists of (variable-length) arrays, need to embed these in some fixed length array

            value = value_dict[file_key]

            is_nested_list = False
            if(isinstance(value[0],list)):
                is_nested_list = True
            elif(isinstance(value[0],np.ndarray)):
                is_nested_list = True

            if(is_nested_list):
                # create an array
                max_length = np.max([len(x) for x in value])
                value_array = np.array([np.pad(row, (0, max_length-len(row))) for row in value])
                f.create_dataset(key,data=value_array,compression='gzip',compression_opts=copts)
            else:
                f.create_dataset(key,data=value,compression='gzip',compression_opts=copts)

        f.close()

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

class PileupOverlaySingle(PileupOverlay):

    """
    An extension of the PileupOverlay class, that just picks
    a single event to overlay. In practice, this can be used
    for "pre-mixing", where this class is sampling from some
    pre-formed pileup events. This may be faster than on-the-fly
    mixing, as it allows one to pre-compute the full pileup.
    """
    def _sample_mu_distribution(self):
        self.mu = 1 # just forces sampling of a single event from the file

class PileupOverlayPtFilter(PileupOverlay):
    """
    A modified version of the PileupOverlay class. This will
    allow for recycling of events if they have leading jet pt
    below some user-defined cut; events with leading jet pt
    above this cut can only be used once -- but once all
    such events have been used once, we then allow recycling
    (so that we don't just run out of high-pT events).

    Warning: Keep in mind that when running in parallel
    (e.g. HTCondor), each job is unaware of what pileup
    events another job may have already used; this only
    limits recycling of high-pt jets *within* a single job.
    """
    # In practice, you'll want to use enough events that
    # there isn't *much* recycling, but this method should
    # still ensure that the effective "recycling rate" of
    # the higher-pt events is comparitively low; at the cost
    # of running more slowly (since this involves an extra layer
    # of jet clustering).

    def __init__(self, pileup_files:Optional[Union[str,list]]=None,rng_seed:int=1,mu_file:str=None, pt_cut=35., precompute=False):

        super(PileupOverlayPtFilter,self).__init__(pileup_files,rng_seed,mu_file)

        self.pt_cut = pt_cut # GeV
        self.jet_finder = TruthJetFinder('anti_kt',0.4)
        self.event_usage_mask = np.full(self.n_total,False,dtype=bool) # keep track of whether or not an event has been used
        self.leading_pt = np.full(self.n_total,-999.)
        self.precompute = precompute

        # Create an internal histogram, to keep track of the leading jet kinematics -- can be interesting for diagnosis of issues.
        self._initialize_histogram()

    def Initialize(self):
        self.jet_finder.SetConfigurator(self.configurator)
        self.jet_finder.Initialize()

        if(self.precompute):
            self._precompute_leading_pt()

    def _initialize_histogram(self):
        binning = (100,0.,100.)
        self.hist_pt = rt.TH1D('h_PileupOverlayPtFilter_leading_pt',';p_{T} [GeV];Count (unweighted)',*binning)

    def _pick_event_indices(self):
        self._sample_mu_distribution()

        if(self.mu > np.sum(self.event_mask)): # running out of events; need to allow recycling
            self._reset_mask()

        if(np.sum(self.event_usage_mask) == len(self.event_usage_mask)):
            self._reset_mask() # TODO: Maybe needs some more work, or a different approach? -Jan

        self.selected_indices = self.rng.choice(self.event_indices[self.event_mask],size=self.mu)
        self.event_usage_mask[self.selected_indices] = True

        self._update_mask()

    def _reset_mask(self):
        self.event_mask = np.full(self.n_total,True,dtype=bool)

    def _fill_leading_pt(self,indices=None):

        if(indices is None):
            indices = self.selected_indices

        pileup_events = self._fetch_events(indices)

        for event,index in zip(pileup_events,indices):
            if(self.leading_pt[index] > 0.): continue

            event_particles = event.particles()
            status = np.array([x.status() for x in event_particles])
            stable_particles = list(itertools.compress(event_particles, status == 1))
            stable_particle_vecs = [ParticleToVector(x) for x in stable_particles]

            momenta = np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_vecs]
                    for method in ['E','Px','Py','Pz']
                ]).T

            self.jet_finder.Process(momenta)
            self.leading_pt[index] = np.max( [momentum[0] for momentum in self.jet_finder.jet_vectors_cyl.values() ] )
            self.hist_pt.Fill(self.leading_pt[index])
            if(self.leading_pt[index] > self.pt_cut):
                self._print('pt > {} GeV @ index = {}'.format(self.pt_cut,index))
        return

    def _update_mask(self):
        # make sure we know the leading jet pt of these events
        self._fill_leading_pt() # only fills if it's not already computed
        pt_mask = self.leading_pt[self.selected_indices] > self.pt_cut
        masked_indices = list(itertools.compress(self.selected_indices,pt_mask))
        self.event_mask[masked_indices] = False

    def _precompute_leading_pt(self):
        # We'll batch things.
        batch_size = 100 # configurable?
        index_ranges = [(i, min(i + batch_size, self.n_total)) for i in range(0, self.n_total, batch_size)]

        self._print('Pre-computing leading jet pT for all input pileup events.')

        for i,index_range in enumerate(index_ranges):
            indices = range(*index_range)
            self._fill_leading_pt(indices)
            printProgressBarColor(i+1,len(index_ranges),'Pre-computing jet pT:',decimals=2)

    def Process(self, inputs, outputs=None):
        super().Process(inputs, outputs)

        # Save the histogram.
        histogram_filename = '{}/pileup_histograms.root'.format(self.outdir)
        self.hist_pt.SaveAs(histogram_filename)

        return

class PileupOverlayFromGenerator(PileupOverlay):
    """
    This is a modification of the PileupOverlay class,
    which generates pileup events on-the-fly using a Pythia8
    instance taken from the configuration of the generation step.
    This has the potential advantage of using whatever RNG seed
    was given to the generator -- so that for parallelized jobs
    (HTCondor), if the generation RNG seeds are different then
    the jobs will each generate their unique pileup events
    on-the-fly.
    """

    def __init__(self, rng_seed:int=1,mu_file:Optional[str]=None, pythia_config_file:Optional[str]=None):

        super(PileupOverlayFromGenerator,self).__init__(None,rng_seed,mu_file)
        self.require_pileup_input = False
        self.pythia_config_file = pythia_config_file
        self.SetUsePhiRotations(False) # no reason to use these as we're generating new events on-the-fly

    def SetGenerator(self,generator:'PythiaGenerator'):
        self.generator = generator

    def Initialize(self):
        """
        Here, we will reconfigure the Pythia8 generator.
        This should be OK since in practice, we're not going
        back to the generation step (so it's OK that we are
        modifying the generator).
        """

        # We need to (re)configure the generator, to make sure that it is making SoftQCD.
        # TODO: This is a bit hacky, maybe could consider reworking it?
        if(self.pythia_config_file is None):
            this_dir = os.path.dirname(os.path.abspath(__file__))
            self.pythia_config_file = str(pathlib.Path('{}/../pythia/pythia_templates/SoftQCD.txt'.format(this_dir)).absolute())
        self.generator.SetPt(-1.,-1.) # remove any \hat{p_{T}} cuts
        self.generator.ConfigPythia(self.pythia_config_file)

    def _pick_event_indices(self):
        self._sample_mu_distribution()
        self.selected_indices = np.full(self.mu,0) # actual entries will be unused, just need the right size

    def _fetch_events(self,indices:Union[int,list,np.ndarray]):
        """
        Redefined w.r.t. PileupOverlay. Here, the provided indices
        argument is almost vestigial; it just tells us how many
        pileup events to generate.
        """
        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)
        n_pileup_events = len(indices)

        return self._generate(n_pileup_events)

    def _generate(self,n_events:int):
        from pyHepMC3 import HepMC3 as hm

        events = []
        for i in range(n_events):
            hepmc_event = self.generator.GenerateSingle() # generate an event!
            events.append(hepmc_event)
        return events
