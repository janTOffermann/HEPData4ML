import ROOT as rt
import numpy as np
import glob,sys,os

sys.path.append("/home/jofferma/projects/HEPData4ML")

from util.hepmc.setup import HepMCSetup
from util.hepmc.readers import ReaderAscii, ReaderRootTree
from typing import List, Union, Tuple, Optional, TYPE_CHECKING

if(TYPE_CHECKING):
    setup = HepMCSetup(verbose=True)
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
    """

    def __init__(self, pileup_files:Optional[Union[str,list]]=None,rng_seed:int=1,mu_file:str=None):

        # Upon start, make sure that hepmc is set up
        self._init_hepmc()

        self.files = None
        self.SetPileupFiles(pileup_files)

        self.files = sorted(glob.glob(pileup_files))

        # determine the total number of events in these files
        self.n_dict = {}
        self.first_idx = {}
        self.n_total = 0
        self._index_files()

        self.event_indices = np.arange(self.n_total)
        self.event_mask = np.full(self.n_total,True,dtype=bool)
        self.selected_indices = None # transient storage for indices selected for a particular event

        self.rng_seed = rng_seed
        self.rng = np.random.default_rng(self.rng_seed)
        self.beam_spot_sigma = None # will store the beam spot size
        self.SetBeamSpotSigma() # sets some sensible defaults

        # Set up the distribution for # of interactions per crossing
        self.mu_file = mu_file
        self.mu_values = None
        self.mu_probs = None
        self._init_mu_distribution()
        self.mu = None # transient storage for the current value of mu

        self.pileup_events = None # transient storage for pileup events

        self.event_buffer_size = 10 # how many combined events to store in memory, before flushing to output file
        #TODO: Need to fix buffer behavior -- the HepMC writers overwrite the output_file each time they're initialized, as opposed to appending

    def _init_hepmc(self):
        setup = HepMCSetup(verbose=True)
        setup.PrepHepMC() # will download/install if necessary
        python_dir = setup.GetPythonDirectory()
        if(python_dir not in sys.path):
            sys.path = [setup.GetPythonDirectory()] + sys.path # prepend, to make sure we pick this one up first
        # now can do "from pyHepMC3 import HepMC3 as hm"

    def SetPileupFiles(self,files:Union[str,list]=None):
        if(files is None):
            return
        if(isinstance(files,str)):
            self.files = glob.glob(files)
        else:
            self.files = files
        self.files = sorted(self.files)
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
        # print('Running InitializeMuDistibution')
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
        self.mu_values = np.arange(mu_min,mu_max)
        self.mu_probs = np.array([self.mu_distribution.GetBinContent(i+1) for i in range(n)],dtype=float)
        self.mu_probs /= np.sum(self.mu_probs)

    def _sample_mu_distribution(self):
        self.mu = int(self.rng.choice(self.mu_values,p=self.mu_probs))

    def _pick_event_indices(self, update_mask=False):
        self._sample_mu_distribution()
        self.selected_indices = self.rng.choice(self.event_indices[self.event_mask],size=self.mu)
        if(update_mask): self._update_mask()

        print('self.selected_indices = ', self.selected_indices)

    def _update_mask(self,indices=None):
        """
        Remove "indices" from pileup event indices that can be sampled from in the future.
        Can be used to achieve sampling without replacement.
        """
        if(indices is None): indices = self.selected_indices # use currently selected indices
        self.event_mask[self.selected_indices] = False

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

        events = [hm.GenEvent()] * len(indices)
        counter = 0
        write_counter = 0

        if(filename.split('.')[-1].lower() == 'root'):
            reader = ReaderRootTree(filename)
        else:
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

    def _fetch_events(self,indices:Union[int,list,np.ndarray]):
        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)

        # determine which files we'll have to pull pileup events from
        file_indices = {}
        for idx in indices:
            file_idx = np.where(idx >= np.array(list(self.first_idx.values()),dtype=int))[0][0]

            if(file_idx not in file_indices.keys()):
                file_indices[file_idx] = []

            file_indices[file_idx].append(idx)

        events = []
        for key,val in file_indices.items():
            events += self._fetch_event_single_file(self.files[key],val)
        return events

    def __call__(self,input_file:str,output_file:Optional[str]=None):
        from pyHepMC3 import HepMC3 as hm
        # from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import ReaderRootTree

        input_file_extension = input_file.split('.')[-1]

        if(output_file is None):
            output_file = input_file.replace('.{}'.format(input_file_extension),'_mixed.{}'.format(input_file_extension))

        # buffer file is only actually relevant for the ASCII mode
        # TODO: Have to check how ROOT writing works -- appends to existing file, or erases it?
        buffer_file = output_file.replace('.{}'.format(input_file_extension),'_buffer.{}'.format(input_file_extension))

        if((self.files is None) or len(self.files) == 0):
            print('Error: No pileup files for PileupOverlay.')
            return

        events = []

        if(input_file_extension.lower() == 'root'):
            reader = ReaderRootTree(input_file)
        else:
            reader = ReaderAscii(input_file)

        i = 0
        while(not reader.failed()):

            evt = hm.GenEvent()
            reader.read_event(evt)

            # print('i = {}, with event {}'.format(i,evt))

            self._pick_event_indices()
            pileup_events = self._fetch_events(self.selected_indices)

            # overlay pileup on this event
            # for now, use automatic displacement -- will use self.beam_spot_sigma
            self._combine_event_with_pileup(evt,pileup_events)

            # store the combined event in memory
            events.append(evt)

            # flush events from memory as needed
            if(len(events) >= self.event_buffer_size):

                # write to output file.
                # TODO: Is append functionality supported by hepmc writers?
                self._flush_to_file(events,output_file)

            i += 1

        # one last flush for any stragglers
        if(len(events) > 0):
            self._flush_to_file(events,output_file)

        try:
            # delete the buffer file
            os.unlink(buffer_file)
        except:
            pass

        return output_file

    def _flush_to_file(self,events:List['hm.GenEvent'],output_file:str,buffername:str=None):
        from pyHepMC3 import HepMC3 as hm
        from pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 import WriterRootTree

        if(output_file.split('.')[-1].lower() == 'root'):
            writer = WriterRootTree(output_file)
        else:
            writer = hm.WriterAscii(output_file) # TODO: double-check if there is append functionality, esp. for ascii (may be missing)

        for evt in events:
            writer.write_event(evt)

        writer.close()
        events.clear()
        return


    def _gaussian(self,x,mu,sig,A=None):
        if(A is None): A = 1. / (np.sqrt(2.0 * np.pi))
        return A * np.exp(-np.square((x - mu) / sig) / 2)

    def _combine_event_with_pileup(self,
        main_event: 'hm.GenEvent',
        pileup_events: Union['hm.GenEvent', List['hm.GenEvent']],
        pileup_displacements: Optional[Union[Tuple[float, float, float, float],
                                        List[Tuple[float, float, float, float]]]] = None,
        auto_displace: bool = True,
        beam_spot_sigma: Optional[Tuple[float,float,float,float]] = None,
        in_place: bool = True
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
        in_place: bool
            If True, modifies main_event in-place. If False, copies main_event first.
            For the foreseen use cases, this should be set to True (it is faster).

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

        # Copy all particles and vertices from main event to combined_event
        # NOTE: Doing things in-place appears to be 30% faster in basic tests. So this is preferred based on usage.
        if(in_place): #
            combined_event = main_event
        else:
            combined_event = self._copy_event_content(main_event)

        # Add each pileup event with displacement
        for pileup_event, (dt, dx, dy, dz) in zip(pileup_events, pileup_displacements):
            self._add_displaced_event(combined_event, pileup_event, dt, dx, dy, dz)

        return combined_event

    def _copy_event_content(self,source_event: 'hm.GenEvent', target_event: 'hm.GenEvent' = None):
        """Copy all particles and vertices from source to target event."""
        from pyHepMC3 import HepMC3 as hm

        if(target_event is None):
            target_event = hm.GenEvent(source_event.momentum_unit(), source_event.length_unit())

        # Copy main event properties
        target_event.set_event_number(source_event.event_number())
        target_event.set_cross_section(source_event.cross_section())
        # target_event.weights = source_event.weights.copy() # TODO: haven't figured out best way to copy this

        # Create mapping from old vertex objects to new vertex objects using list index
        vertex_map = {}
        vertices_list = list(source_event.vertices())

        # Copy vertices first
        for i, vertex in enumerate(vertices_list):
            new_vertex = hm.GenVertex(vertex.position())
            new_vertex.status = vertex.status()
            target_event.add_vertex(new_vertex)
            vertex_map[i] = new_vertex

        # Copy particles and establish relationships
        for particle in source_event.particles():
            new_particle = hm.GenParticle(
                particle.momentum(),
                particle.pid(),
                particle.status()
            )
            new_particle.set_generated_mass(particle.generated_mass())

            # Set production vertex if it exists
            if particle.production_vertex():
                # Find the index of the production vertex in the original list
                try:
                    vertex_idx = vertices_list.index(particle.production_vertex())
                    prod_vertex = vertex_map[vertex_idx]
                    prod_vertex.add_particle_out(new_particle)
                except ValueError:
                    pass # expected to be triggered by beam particles
                    # print(f"Warning: Could not find production vertex for particle {particle.pid}")

            # Set end vertex if it exists
            if particle.end_vertex():
                try:
                    vertex_idx = vertices_list.index(particle.end_vertex())
                    end_vertex = vertex_map[vertex_idx]
                    end_vertex.add_particle_in(new_particle)
                except ValueError:
                    print(f"Warning: Could not find end vertex for particle {particle.pid}")
        return target_event

    def _add_displaced_event(self,
                            target_event: 'hm.GenEvent',
                            pileup_event: 'hm.GenEvent',
                            dt: float, dx: float, dy: float, dz: float):
        """Add a pileup event to the target event with vertex displacement."""
        from pyHepMC3 import HepMC3 as hm

        # Create mapping from old vertex objects to new vertex objects using list index
        vertex_map = {}
        vertices_list = list(pileup_event.vertices())

        # Copy vertices with displacement
        for i, vertex in enumerate(vertices_list):
            # Apply displacement to vertex position
            old_pos = vertex.position
            new_position = hm.FourVector(
                old_pos.x + dx,
                old_pos.y + dy,
                old_pos.z + dz,
                old_pos.t + dt
            )

            new_vertex = hm.GenVertex(new_position)
            new_vertex.set_status(vertex.status())
            target_event.add_vertex(new_vertex)
            vertex_map[i] = new_vertex

        # Copy particles and establish relationships
        for particle in pileup_event.particles():
            new_particle = hm.GenParticle(
                particle.momentum(),
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
                    # print(f"Warning: Could not find production vertex for particle {particle.pid}")

            # Set end vertex if it exists
            if particle.end_vertex():
                try:
                    vertex_idx = vertices_list.index(particle.end_vertex())
                    end_vertex = vertex_map[vertex_idx]
                    end_vertex.add_particle_in(new_particle)
                except ValueError:
                    print(f"Warning: Could not find end vertex for particle {particle.pid}")

def main(args):
    pileup = PileupOverlay(pileup_files="/home/jofferma/projects/HEPData4ML/test/pileup/pileup.root")
    _ = pileup(input_file="/home/jofferma/projects/HEPData4ML/test/pileup/events.root",output_file="/home/jofferma/projects/HEPData4ML/test/pileup/output.root")


if(__name__=='__main__'):
    main(sys.argv)






