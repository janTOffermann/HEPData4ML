import ROOT as rt
import pyhepmc as hep
import numpy as np
import glob,sys,os
from util.hepmc.hepmc import PyHepMCOutputAscii
from typing import List, Union, Tuple, Optional

class PileupOverlay:
    """
    This class performs on-the-fly mixing, to mix in events from some pileup HepMC3 file(s)
    into "main events" (from some other HepMC file).
    It creates mixed events, which are output to a new file.
    """

    def __init__(self, pileup_files:Optional[Union[str,list]]=None,rng_seed:int=1,mu_file:str=None):

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

    def SetPileupFiles(self,files:Union[str,list]=None):
        if(files is None):
            return
        if(isinstance(files,str)):
            self.files = glob.glob(files)
        else:
            self.files = files
        self.files = sorted(self.files)
        return

    def _index_files(self):
        for file in self.files:
            with open(file,'r') as f:
                n =  len([x for x in f.readlines() if 'E ' in x])
                self.first_idx[file] = self.n_total
                self.n_dict[file] = n
                self.n_total += n

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

        print(self.selected_indices)

    def _update_mask(self,indices=None):
        """
        Remove "indices" from pileup event indices that can be sampled from in the future.
        Can be used to achieve sampling without replacement.
        """
        if(indices is None): indices = self.selected_indices # use currently selected indices
        self.event_mask[self.selected_indices] = False

    def _fetch_event_single_file(self,filename:str,indices:Union[int,list,np.ndarray]):
        # TODO: Would be nice to find alternative, but for now we have to iterate through
        #       the given HepMC3 file since pyhepmc treats files as iterable but does
        #       not have any way to just fetch the i'th event. This will scale poorly
        #       with large HepMC3 files.

        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)
        indices = sorted(indices)

        events = [hep.GenEvent()] * len(indices)
        counter = 0

        with hep.io.ReaderAscii(filename) as f:
            for i,evt in enumerate(f):
                if(i < indices[0]): continue
                if(i in indices):
                    events[counter] = evt
                    counter += 1
                if(i > indices[-1]):
                    break
        return events

    def _fetch_events(self,indices:Union[int,list,np.ndarray]):
        if(isinstance(indices,int)):
            indices = np.atleast_1d(int)

        # determine which files we'll have to pull pileup events from
        file_indices = {}
        for idx in indices:
            file_idx = np.where(idx > np.array(list(self.first_idx.values()),dtype=int))[0][0]

            if(file_idx not in file_indices.keys()):
                file_indices[file_idx] = []

            file_indices[file_idx].append(idx)

        events = []
        for key,val in file_indices.items():
            events += self._fetch_event_single_file(self.files[key],val)
        return events

    def __call__(self,input_file:str,output_file:Optional[str]=None):

        if(output_file is None):
            output_file = input_file.replace('.hepmc','_mixed.hepmc')
        buffer_file = output_file.replace('.hepmc','_buffer.hepmc')

        if((self.files is None) or len(self.files) == 0):
            print('Error: No pileup files for PileupOverlay.')
            return

        events = []
        with hep.io.ReaderAscii(input_file) as f:
            for i,evt in enumerate(f):

                # get the pileup
                self._pick_event_indices()
                pileup_events = self._fetch_events(self.selected_indices)

                # overlay pileup on this event
                # for now, use automatic displacement -- will use self.beam_spot_sigma
                self._combine_event_with_pileup(evt,pileup_events)

                # store the combined event in memory
                events.append(evt)

                # flush events from memory as needed
                if(len(events) >= self.event_buffer_size):
                    PyHepMCOutputAscii(events,buffername=buffer_file,filename=output_file)
                    events.clear()

        # one last flush for any stragglers
        if(len(events) > 0):
            PyHepMCOutputAscii(events,buffername=buffer_file,filename=output_file)

        # delete the buffer file
        os.unlink(buffer_file)
        return output_file

    def _gaussian(self,x,mu,sig,A=None):
        if(A is None): A = 1. / (np.sqrt(2.0 * np.pi))
        return A * np.exp(-np.square((x - mu) / sig) / 2)

    def _combine_event_with_pileup(self,
        main_event: hep.GenEvent,
        pileup_events: Union[hep.GenEvent, List[hep.GenEvent]],
        pileup_displacements: Optional[Union[Tuple[float, float, float, float],
                                        List[Tuple[float, float, float, float]]]] = None,
        auto_displace: bool = True,
        beam_spot_sigma: Optional[Tuple[float,float,float,float]] = None,
        in_place: bool = True
    ) -> hep.GenEvent:
        """
        Combine a main event with pileup events, applying vertex displacements.

        Parameters:
        -----------
        main_event : hep.GenEvent
            The primary hard-scatter event
        pileup_events : hep.GenEvent or List[hep.GenEvent]
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
        hep.GenEvent
            Combined event with displaced pileup vertices
        """

        # Ensure pileup_events is a list
        if isinstance(pileup_events, hep.GenEvent):
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

    def _copy_event_content(self,source_event: hep.GenEvent, target_event: hep.GenEvent = None):
        """Copy all particles and vertices from source to target event."""

        if(target_event is None):
            target_event = hep.GenEvent(source_event.momentum_unit, source_event.length_unit)

        # Copy main event properties
        target_event.event_number = source_event.event_number
        target_event.cross_section = source_event.cross_section
        target_event.weights = source_event.weights.copy()

        # Create mapping from old vertex objects to new vertex objects using list index
        vertex_map = {}
        vertices_list = list(source_event.vertices)

        # Copy vertices first
        for i, vertex in enumerate(vertices_list):
            new_vertex = hep.GenVertex(vertex.position)
            new_vertex.status = vertex.status
            target_event.add_vertex(new_vertex)
            vertex_map[i] = new_vertex

        # Copy particles and establish relationships
        for particle in source_event.particles:
            new_particle = hep.GenParticle(
                particle.momentum,
                particle.pid,
                particle.status
            )
            new_particle.generated_mass = particle.generated_mass

            # Set production vertex if it exists
            if particle.production_vertex:
                # Find the index of the production vertex in the original list
                try:
                    vertex_idx = vertices_list.index(particle.production_vertex)
                    prod_vertex = vertex_map[vertex_idx]
                    prod_vertex.add_particle_out(new_particle)
                except ValueError:
                    pass # expected to be triggered by beam particles
                    # print(f"Warning: Could not find production vertex for particle {particle.pid}")

            # Set end vertex if it exists
            if particle.end_vertex:
                try:
                    vertex_idx = vertices_list.index(particle.end_vertex)
                    end_vertex = vertex_map[vertex_idx]
                    end_vertex.add_particle_in(new_particle)
                except ValueError:
                    print(f"Warning: Could not find end vertex for particle {particle.pid}")
        return target_event

    def _add_displaced_event(self,
                            target_event: hep.GenEvent,
                            pileup_event: hep.GenEvent,
                            dt: float, dx: float, dy: float, dz: float):
        """Add a pileup event to the target event with vertex displacement."""

        # Create mapping from old vertex objects to new vertex objects using list index
        vertex_map = {}
        vertices_list = list(pileup_event.vertices)

        # Copy vertices with displacement
        for i, vertex in enumerate(vertices_list):
            # Apply displacement to vertex position
            old_pos = vertex.position
            new_position = hep.FourVector(
                old_pos.x + dx,
                old_pos.y + dy,
                old_pos.z + dz,
                old_pos.t + dt
            )

            new_vertex = hep.GenVertex(new_position)
            new_vertex.status = vertex.status
            target_event.add_vertex(new_vertex)
            vertex_map[i] = new_vertex

        # Copy particles and establish relationships
        for particle in pileup_event.particles:
            new_particle = hep.GenParticle(
                particle.momentum,
                particle.pid,
                particle.status
            )
            new_particle.generated_mass = particle.generated_mass

            # Set production vertex if it exists
            if particle.production_vertex:
                try:
                    vertex_idx = vertices_list.index(particle.production_vertex)
                    prod_vertex = vertex_map[vertex_idx]
                    prod_vertex.add_particle_out(new_particle)
                except ValueError:
                    pass # expected to be triggered, by beam particles
                    # print(f"Warning: Could not find production vertex for particle {particle.pid}")

            # Set end vertex if it exists
            if particle.end_vertex:
                try:
                    vertex_idx = vertices_list.index(particle.end_vertex)
                    end_vertex = vertex_map[vertex_idx]
                    end_vertex.add_particle_in(new_particle)
                except ValueError:
                    print(f"Warning: Could not find end vertex for particle {particle.pid}")

def main(args):
    pileup = PileupOverlay(pileup_files="/home/jofferma/projects/HEPData4ML/tmp/pileup_events.hepmc3")

    _ = pileup()


if(__name__=='__main__'):
    main(sys.argv)






