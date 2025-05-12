import sys, glob, uuid, itertools
import numpy as np
import h5py as h5
import ROOT as rt
import uproot as ur
# import pyhepmc as hep
import subprocess as sub
from util.calcs import Calculator
from util.qol_utils.qol_util import printProgressBarColor, RN
from util.qol_utils.pdg import pdg_plotcodes, pdg_names, FillPdgHist
# from util.particle_selection.algos import IsNeutrino
from util.fastjet import FastJetSetup
from util.hepmc import ExtractHepMCEvents

def embed_array(array, target_shape,padding_value=0):
    """
    A generic function for embedding an input array into some target shape.
    The array will be truncated or padded as needed.

    """
    # Create the output array
    result = np.full(target_shape, padding_value, dtype=array.dtype)

    # Calculate the effective shape (minimum dimensions)
    effective_shape = tuple(min(s, t) for s, t in zip(array.shape, target_shape))

    slices = tuple(slice(0, dim) for dim in effective_shape)
    result[slices] = array[slices]
    return result

class Processor:
    """
    Convert HepMC or Delphes/ROOT file, representing an event, into an HDF5 file where each entry/event corresponds
    with a single jet. (Can also skip jet clustering entirely).
    """
    def __init__(self, configurator):
        self.configurator = configurator
        self.delphes = False # will be set to True if SetDelphesFiles() is called

        self.SetProgressBarPrefix('\tClustering jets & preparing data:')
        self.suffix = 'Complete'
        self.bl = 50
        self.verbose = False

        self.SetupFastJet()
        sys.path.append(self.fastjet_dir)
        import fastjet as fj # This is where fastjet is really imported & cached by Python, but there are other import statements peppered throughout since this has limited scope.
        #TODO: Is there a nicer way to handle this import, which requires a filepath passed to this class's constructor?

        self.jet_config, self.jetdef = self.InitFastJet()
        self.jet_selection = self.jet_config['jet_selection']
        if(self.jet_selection is not None):
            self.jet_selection.Initialize(self.configurator)
        self.n_constituents = self.configurator.GetNPars()['jet_n_par']

        self.outdir = ''

        self.hist_filename = 'hists.root'
        self.diagnostic_plots = True
        self.InitializeHists()

        self.separate_truth_particles = False
        self.n_separate_truth_particles = -1

        self.cluster_sequence = None
        self.jets = None
        self.jets_filtered = None

        self.SetRecordFinalStateIndices(False)
        self.SetPostProcessing()

        self.SetRecordFullFinalState(False) # by default we don't want this, it may significantly increase filesize!

        self.calculator = Calculator(use_vectorcalcs=self.configurator.GetUseVectorCalcs())

        self.data = {}

        # Various integers for buffer size, number of particles read into memory from HepMC, number saved to file, etc.
        self.nevents_per_chunk = int(1e3) # Can be configured. Affects memory footprint.
        self.nparticles_max = int(2e3) # TODO: This is some hardcoded max number of particles to be read in from HepMC. Should be plenty.
        self.nparticles_stable = self.configurator.GetNPars()['n_stable']
        self.nparticles_truth_selected = self.configurator.GetNPars()['n_truth']
        self.n_delphes = self.configurator.GetNPars()['n_delphes']

        # self.n_constituents = self.configurator.GetNPars()['jet_n_par']

    def SetNentriesPerChunk(self,val):
        self.nevents_per_chunk = int(val)

    def SetDelphesFiles(self,val):
        if(len(val) > 0):
            self.delphes = True
            self.delphes_files = val

    def ResetDataContainers(self):

        self.event_record_vector_components = {
            'px':None,
            'py':None,
            'pz':None,
            'e':None,
            'pt':None,
            'eta':None,
            'phi':None,
            'm':None,
            'pdg':None,
            'status':None
        }

        self.stable_particle_components = {
            'px':None,
            'py':None,
            'pz':None,
            'e':None,
            'pt':None,
            'eta':None,
            'phi':None,
            'm':None,
            'pdg':None
        }

        self.truth_particle_components = {
            'px':None,
            'py':None,
            'pz':None,
            'e':None,
            'pt':None,
            'eta':None,
            'phi':None,
            'm':None,
            'pdg':None
        }

    def SetRecordFullFinalState(self,flag):
        self.record_full_final_state = flag

    def SetProgressBarPrefix(self,text):
        self.prefix_level1 = text

    def SetRecordFinalStateIndices(self,flag):
        self.record_final_state_indices = flag

    def SetPostProcessing(self,post_proc=None):
        if(post_proc is None): post_proc = self.configurator.GetPostProcessing()
        self.post_processing = post_proc
        if(self.post_processing is None):
            return

        # Determine whether or not to generate files containing indices
        # of particles that were both in the truth and final-state selections.
        # Useful for things like particle daughter tracing.
        # TODO: There might be a neater way to handle this.
        self.SetRecordFinalStateIndices(False)
        for p in post_proc:
            if(p is None): continue
            if(p.RequiresIndexing()):
                self.SetRecordFinalStateIndices(True)
                break

    def SetStatsFile(self,filename):
        self.stats_file = filename

    def SetSeparateTruthParticles(self,flag):
        self.separate_truth_particles = flag

    def SetNSeparateTruthParticles(self,n):
        self.n_separate_truth_particles = n

    def SetVerbosity(self,flag):
        self.verbose = flag

    def SetOutputDirectory(self,outdir):
        self.outdir = outdir

    def SetHistFilename(self,name):
        self.hist_filename = name

    def SetDiagnosticPlots(self,flag):
        self.diagnostic_plots = flag

    def InitializeHists(self):
        self.hists = []
        self.pdg_hist = rt.TH1D(RN(),'Selected final-state particles;Particle;Count',47,0.,47.)
        self.hists.append(self.pdg_hist)
        self.jet_truth_dr_hists = {} # Entries get initialized and filled in FillDataBuffer

        self.jet_hists = { # Jet kinematic histograms. Initialized/filled in FillDataBuffer.
            'pt'  : rt.TH1D(RN(), 'Jet p_{T};p_{T} [GeV];Count', 500,0.,1000.),
            'eta' : rt.TH1D(RN(), 'Jet #eta;#eta;Count',         200,-4.,4.),
            'phi' : rt.TH1D(RN(), 'Jet #phi;#phi;Count',         200,-np.pi,np.pi),
            'm'   : rt.TH1D(RN(), 'Jet m;m [GeV];Count',         250,0.,500.),
            'e'   : rt.TH1D(RN(), 'Jet E;E [GeV];Count',         500,0.,1000.),
            'et'  : rt.TH1D(RN(), 'Jet E_{T};E [GeV];Count',     500,0.,1000.)
        }

        self.pmu_sum_hists = {} # Ditto, for the sum of jet constituents.
        for key,hist in self.jet_hists.items():
            h = rt.TH1D(hist)
            h.SetName(RN())
            h.SetTitle(h.GetTitle().replace('Jet','#sum constituents'))
            self.pmu_sum_hists[key] = h

        self.jet_diff_hists = {
            'dm': rt.TH1D(RN(),  ';(m_{consts} - m_{j}) / m_{j};Count'      , 400,-1.,1.),
            'de': rt.TH1D(RN(),  ';(E_{consts} - E_{j}) / E_{j};Count'      , 400,-1.,1.),
            'dpt': rt.TH1D(RN(), ';(p_{consts} - p_{T,j}) / p_{T,j};Count', 400,-1.,1.),
            'det': rt.TH1D(RN(), ';(E_{consts} - E_{T,j}) / E_{T,j};Count', 400,-1.,1.),
        } # Kinematic differences between the jets, and the sum of their constituents. (Sanity check!)

        self.jet_diff_hists_2d = {
            'dm': rt.TH2D(RN(),  ';m_{j};(m_{consts} - m_{j}) / m_{j};Count',         250,0.,500. ,400,-1.,1.),
            'de': rt.TH2D(RN(),  ';E_{j};(E_{consts} - E_{j}) / E_{j};Count',         500,0.,1000.,400,-1.,1.),
            'dpt': rt.TH2D(RN(), ';p_{T,j};(p_{consts} - p_{T,j}) / p_{T,j};Count', 500,0.,1000.,400,-1.,1.),
            'det': rt.TH2D(RN(), ';E_{T,j};(E_{T,consts} - E_{T,j}) / E_{T,j};Count', 500,0.,1000.,400,-1.,1.),
        } # Kinematic differences between the jets, and the sum of their constituents. (Sanity check!)

    def OutputHistograms(self):
        for key,h in self.jet_truth_dr_hists.items():
            self.hists.append(h)

        for key,h in self.pmu_sum_hists.items():
            self.hists.append(h)

        for key,h in self.jet_hists.items():
            self.hists.append(h)

        for key,h in self.jet_diff_hists.items():
            self.hists.append(h)

        for key,h in self.jet_diff_hists_2d.items():
            self.hists.append(h)

        rt.gStyle.SetOptStat(0)
        offset = 10 # TODO: This is a bit hacky
        hist_filename = '{}/{}'.format(self.outdir,self.hist_filename)
        f = rt.TFile(hist_filename,'UPDATE')

        # We separately handle the pdg_hist (first in list), since it requires
        # modifying the axis labels after drawing.
        for i,h in enumerate(self.hists):
            canvas_name = 'c_{}'.format(offset + i)
            hist_name = 'h_{}'.format(offset + i)
            c = rt.TCanvas(canvas_name,'',800,600)

            if(i == 0): # pdg_hist

                pdg_names_inv = {}
                for key,val in pdg_plotcodes.items():
                    code = val
                    name = pdg_names[key]
                    pdg_names_inv[code] = name

                n = h.GetNbinsX()
                xaxis = h.GetXaxis()
                for j in range(n):
                    if(j == 0): name = ''
                    else:
                        try: name = pdg_names_inv[j]
                        except: name = ''
                    xaxis.SetBinLabel(j+1,name)
                xaxis.SetLabelSize(2.0e-2)
                rt.gPad.SetGrid(1,0)
            else: rt.gPad.SetGrid()

            if(type(h) in [rt.TH1F, rt.TH1D, rt.TH1I, rt.TH1S]): h.Draw('HIST')
            else:
                h.Draw('COLZ')
                rt.gPad.SetRightMargin(0.2)
            rt.gPad.SetLogy()
            c.Draw()
            f.cd() # unnecessary
            # c.Write(canvas_name)
            h.Write(hist_name)
        f.Close()
        return

    def Process(self, hepmc_files, h5_file=None, verbosity=0):
        if(type(hepmc_files) == list): hepmc_files = ['{}/{}'.format(self.outdir,x) for x in hepmc_files]
        else: hepmc_files = '{}/{}'.format(self.outdir,hepmc_files)

        if(h5_file is None):
            if(type(hepmc_files) == list): h5_file = hepmc_files[0]
            else: h5_file = hepmc_files
            if(self.delphes): h5_file =  h5_file.replace('*','').replace('.root','.h5')
            else: h5_file =  h5_file.replace('*','').replace('.root','.h5') # TODO: Probably a bug
        else:
            h5_file = '{}/{}'.format(self.outdir,h5_file)

        if(type(hepmc_files) == str): hepmc_files = glob.glob(hepmc_files,recursive=True)

        if(self.delphes):
            # NOTE: It's important that the Delphes files and truth HepMC files line up!
            #       The way we usually do this, it will be guaranteed.
            #       However, if you borrow functions from here you will have to keep this in mind,
            #       as things will go wrong if final_state_fiels and delphes_files aren't sorted
            #       the same way (or are of different lengths).
            delphes_arr,var_map = self.PrepDelphesArrays()
            nentries = len(delphes_arr)

        ## Extract final-state truth particle info from the HepMC files.
        hepmc_events, nentries = ExtractHepMCEvents(hepmc_files,get_nevents=True)

        # self.PrepDataBuffer(self.nentries_per_chunk,self.separate_truth_particles,self.n_separate_truth_particles,self.record_final_state_indices,self.record_full_final_state)

        # Some indexing preparation for writing in chunks.
        start_idxs,stop_idxs,ranges = self.PrepIndexRanges(nentries,self.nevents_per_chunk)
        nchunks = len(ranges)

        dsets = None

        # if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        for i in range(nchunks):
            # Clear the buffer (for safety).
            for key in self.data.keys(): self.data[key][:] = 0

            # TODO: Start here. Things like "final_state_events" will become "events" etc.

            # Extract all the particles from the HepMC events, into memory.
            particles = self.ExtractHepMCParticles(hepmc_events[start_idxs[i]:stop_idxs[i]],self.nparticles_max)

            # 0) Write some keys that are the same across all events in this chunk.
            self.WriteToDataBuffer(None,'is_signal',self.configurator.GetSignalFlag())

            for j,event_particles in enumerate(particles): # Loop over events in this chunk

                # 1) Save the stable truth-level particles from the event.
                status = np.array([x.status for x in event_particles])
                stable_particles = list(itertools.compress(event_particles, status == 1))
                stable_particle_vecs = [rt.Math.PxPyPzEVector(x.momentum.px, x.momentum.py, x.momentum.pz, x.momentum.e) for x in stable_particles]

                # print('[{}], {} stable particles.'.format(j,len(stable_particle_vecs)))

                self.WriteToDataBuffer(j,'Pmu_Nobj',len(stable_particle_vecs))

                self.WriteToDataBuffer(j, 'Pmu', np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_vecs]
                    for method in ['E','Px','Py','Pz']
                ]).T,
                                       dimensions={1:self.nparticles_stable}
                )

                self.WriteToDataBuffer(j, 'Pmu_cyl', np.vstack([
                    [getattr(vec, method)() for vec in stable_particle_vecs]
                    for method in ['Pt','Eta','Phi','M']
                ]).T,
                                       dimensions={1:self.nparticles_stable}
                )

                self.WriteToDataBuffer(j,'Pmu_Pdg',[x.pid for x in stable_particles],
                                       dimensions={1:self.nparticles_stable}
                )

                # 2) Extract the filtered truth record from the events.
                truth_selected_status = np.full(len(event_particles),True) # TODO: New truth selection goes here.
                truth_selected_particles = list(itertools.compress(event_particles, truth_selected_status == 1))
                truth_selected_particle_vecs = [rt.Math.PxPyPzEVector(x.momentum.px, x.momentum.py, x.momentum.pz, x.momentum.e) for x in truth_selected_particles]

                self.WriteToDataBuffer(j,'truth_Pmu_Nobj',len(truth_selected_particle_vecs))

                self.WriteToDataBuffer(j, 'truth_Pmu', np.vstack([
                    [getattr(vec, method)() for vec in truth_selected_particle_vecs]
                    for method in ['E','Px','Py','Pz']
                ]).T,
                                       dimensions={1:self.nparticles_truth_selected}
                )

                self.WriteToDataBuffer(j, 'truth_Pmu_cyl', np.vstack([
                    [getattr(vec, method)() for vec in truth_selected_particle_vecs]
                    for method in ['Pt','Eta','Phi','M']
                ]).T,
                                        dimensions={1:self.nparticles_truth_selected}
                )

                self.WriteToDataBuffer(j,'truth_Pmu_Pdg',[x.pid for x in truth_selected_particles],
                                        dimensions={1:self.nparticles_truth_selected}
                )

            # 3) Optional: Extract the full event record, store it separately. This is both particles and vertices.

            # 4) If Delphes was run, we will also extract the relevant information.
            #    Note that PrepDelphesArrays() has been called earlier, if delphes=True.
            #    That's where the Delphes ROOT files are already read and prepared for access
            #    via uproot.
            #
            #    Note that we choose to store these objects as four-momenta, rather than explicitly
            #    storing their components. This is possibly a relevant detail as some objects have
            #    attributes such as pt, while others have Et. These are the same *if* we assume
            #    the objects themselves to be massless.
            if(self.delphes):
                for j in range(len(particles)): # TODO: reusing len(particles) (= number of events in chunk), OK but looks kind of hacky
                    for k,delphes_type in enumerate(var_map.keys()): # loop over different kinds of Delphes collections

                        delphes_pt  = delphes_arr[var_map[delphes_type]['pt' ]][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                        delphes_eta = delphes_arr[var_map[delphes_type]['eta']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                        delphes_phi = delphes_arr[var_map[delphes_type]['phi']][start_idxs[i]:stop_idxs[i]][j].to_numpy()
                        delphes_m   = np.zeros(delphes_pt.shape)

                        self.WriteToDataBuffer(j, '{}.Pmu'.format(delphes_type), np.vstack([
                            delphes_pt,
                            delphes_eta,
                            delphes_phi,
                            delphes_m
                        ]).T,
                                            dimensions={1:self.n_delphes[k]}
                        )

            # We have now filled a chunk, time to write it.
            # If this is the first instance of the loop, we will initialize the HDF5 file.
            # NOTE: We assume that after this first loop, we've generated all the necessary keys.
            #       Probably a safe assumption for now.

            if(i == 0):
                dsets = self.PrepH5File(h5_file,nentries,self.data)

            with h5.File(h5_file, 'a') as f:
                for key in dsets.keys():
                    dset = f[key]
                    dset[start_idxs[i]:stop_idxs[i]] = self.data[key][:ranges[i]]
            # if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        # if(self.diagnostic_plots): self.OutputHistograms()
        return h5_file

    # --- Utility functions for the class ---
    def SetupFastJet(self):
        verbose = self.configurator.GetPrintFastjet()
        self.fastjet_setup = FastJetSetup(self.configurator.GetFastjetDirectory(),full_setup=True,verbose=verbose)
        self.configurator.SetPrintFastjet(False)
        self.fastjet_dir = self.fastjet_setup.GetPythonDirectory()
        return

    def InitFastJet(self):
        import fastjet as fj # hacky, but will work because SetupFastJet() was run in __init__()
        # Print the FastJet banner -- it's unavoidable (other packages don't do this!).
        fj.ClusterSequence.print_banner()
        jet_config = self.configurator.GetJetConfig()
        jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
        return jet_config,jetdef

    def ClusterJets(self,vecs, jet_config=None):
        if(jet_config is None): jet_config = self.jet_config
        import fastjet as fj # hacky, but will work because SetupFastJet() was run in __init__()

        cartesian = not self.delphes
        if(cartesian): # vecs has format (E,px,py,pz) -- FastJet uses (px,py,pz,E) so we must modify it. Using np.roll.
            pj = [fj.PseudoJet(*x) for x in np.roll(vecs,-1,axis=-1)]

        else: # vecs has format(pt,eta,phi,m)
            pj = [fj.PseudoJet(0.,0.,0.,0.) for i in range(len(vecs))]
            for i,x in enumerate(vecs):
                pj[i].reset_momentum_PtYPhiM(*x) # TODO: massless -> okay to use eta for Y

        # Attach indices to the pseudojet objects, so that we can trace them through jet clustering.
        # Indices will correspond to the order they were input (with zero-indexing).
        for i,pseudojet in enumerate(pj):
            pseudojet.set_user_index(i)

        # selector = fj.SelectorPtMin(jet_config['jet_min_pt']) & fj.SelectorAbsEtaMax(jet_config['jet_max_eta'])
        # Note: Switched from the old method, this is more verbose but seems to do the same thing anyway.
        self.cluster_sequence = fj.ClusterSequence(pj, self.jetdef) # member of class, otherwise goes out-of-scope when ref'd later
        self.jets = self.cluster_sequence.inclusive_jets()

        # Now we will apply our (optional) pt and eta cuts to the jets.
        # TODO: Should be achievable with Fastjet selector classes too.
        jet_eta = np.array([jet.eta() for jet in self.jets])
        jet_pt  = np.array([jet.pt()  for jet in self.jets])

        selected_eta = np.where(np.abs(jet_eta) <= jet_config['jet_max_eta'])[0]
        selected_pt  = np.where(jet_pt          >= jet_config['jet_min_pt'] )[0]
        selected = np.sort(np.intersect1d(selected_eta, selected_pt)) # TODO: Is the sort needed?
        self.jets_filtered = [self.jets[x] for x in selected]

    def FetchJetConstituents(self,jet,n_constituents):
        pt,px,py,pz,e = np.hsplit(np.array([[x.pt(), x.px(),x.py(),x.pz(),x.e()] for x in jet.constituents()]),5)
        pt = pt.flatten() # use this for sorting
        px = px.flatten()
        py = py.flatten()
        pz = pz.flatten()
        e  =  e.flatten()

        # The indices of the jet constituents, corresponding with the order in which they
        # were passed to jet clustering.
        indices = np.array([x.user_index() for x in jet.constituents()],dtype=int).flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt)
        l = int(np.minimum(n_constituents,len(pt)))

        px = px[sorting][:l]
        py = py[sorting][:l]
        pz = pz[sorting][:l]
        e  =  e[sorting][:l]
        indices = indices[sorting][:l]

        vecs = np.array([e,px,py,pz],dtype=np.dtype('f8')).T
        return vecs, indices

    def SelectFinalStateParticles(self, truth_particles, j, debug_print=False):
        # Note: "vecs" will be Cartesian if not using Delphes. If using Delphes, it will be in (pt,eta,phi,m).

        if(self.delphes):
            vecs = np.vstack((
                self.stable_particle_components['pt'],
                self.stable_particle_components['eta'],
                self.stable_particle_components['phi'],
                self.stable_particle_components['m'],
            )).T
        else:
            vecs = np.vstack((
                self.stable_particle_components['e'],
                self.stable_particle_components['px'],
                self.stable_particle_components['py'],
                self.stable_particle_components['pz'],
            )).T

        vecs_copy = np.copy(vecs) # TODO: Is this necessary?

        if(self.jet_selection is None):
            # If we're using Delphes, vecs is in cylindrical and we need to convert it to Cartesian since we're not
            # calling the ClusterJets() routine (which internally handles the conversion via FastJet).
            if(self.delphes): vecs = self.calculator.PtEtaPhiMToEPxPyPz(vecs[:,0], vecs[:,1],vecs[:,2],vecs[:,3]) #TODO: Is this conditional okay?
            self.FillDataBuffer(j,vecs,None,truth_particles)
        else:
            # Any necessary cylindrical->Cartesian conversion (i.e. Delphes case) is handled within ClusterJets().
            self.ClusterJets(vecs)
            jets = self.jets_filtered

            if(len(jets) == 0): # If there are no jets, throw out this event.
                self.data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return

            try: truth = truth_particles[j]
            except: truth = None # needed in case there was no truth selection

            # Now we apply our jet selection algorithm.
            selected_jet_idx = self.jet_selection(truth=truth, jets=jets)

            if(selected_jet_idx < 0): # If we didn't select any jet, throw out this event.
                self.data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return
            jet = jets[selected_jet_idx]

            # Get the constituents of our selected jet, as well as their indices w.r.t. the final state.
            selected_vecs, selected_indices = self.FetchJetConstituents(jet,self.n_constituents)

            self.FillDataBuffer(j,selected_vecs,jet,truth_particles,final_state_indices=selected_indices)

        # If PDG code information was provided, let's determine the codes of the selected particles.
        # TODO: This will involve looping -> will be intensive. Is there a better way?
        # See PseudoJet::set_user_index : http://www.fastjet.fr/repo/fastjet-doc-3.4.0.pdf
        pdg = self.stable_particle_components['pdg']
        if(pdg is not None and self.diagnostic_plots):
            pdg_matching_tolerance = 1.0e-4
            selected_pdgs = np.zeros(pdg.shape,dtype=bool)
            candidate_particles = vecs_copy
            for particle in vecs:
                for k,candidate_particle in enumerate(candidate_particles):
                    if(selected_pdgs[k] == 1): continue
                    match_value = np.dot(candidate_particle-particle,candidate_particle-particle)
                    if match_value < pdg_matching_tolerance:
                        selected_pdgs[k] = 1
                        break
            selected_pdgs = pdg[selected_pdgs]
            FillPdgHist(self.pdg_hist,selected_pdgs)
        return

    def AddKeyToDataBuffer(self,key,value,dtype=None,dimensions=None):
        if(key in self.data.keys()):
            return
        value_array = np.asarray(value)

        if(dtype is None):
            dtype = value_array.dtype

        if value_array.shape == ():  # Scalar value
            # Create shape (N,)
            self.data[key] = np.zeros(self.nevents_per_chunk, dtype=dtype)
        else:
            # Create the buffer shape.
            # The user can optionally specify dimensions via a dictionary,
            # otherwise they are inferred.
            buffer_shape = list((self.nevents_per_chunk,) + value_array.shape)
            if(dimensions is not None):
                for idx,val in dimensions.items():
                    try:
                        buffer_shape[idx] = val
                    except:
                        pass # TODO: Add warning
            self.data[key] = np.zeros(buffer_shape, dtype=dtype)
        print('Added key {} to buffer.'.format(key))
        return

    def WriteToDataBuffer(self,event_index,key,value,dtype=None,dimensions=None):
        if(key not in self.data.keys()):
            self.AddKeyToDataBuffer(key,value,dtype,dimensions)

        value_array = np.asarray(value) # TODO: not sure if needed?
        # print('for key = {}, value_array has shape {}'.format(key,value_array.shape))
        # print('self.nentries_per_chunk = ',self.nevents_per_chunk)
        # print('\tself.data[{}] has shape {}'.format(key,self.data[key].shape))
        if(event_index is not None):

            self.data[key][event_index] = embed_array(value_array,self.data[key][event_index].shape)
        else:
            self.data[key][:] = embed_array(value_array,self.data[key].shape)
        return

    def FillDataBuffer(self,j,vecs,jet,truth_particles,final_state_indices=None):
        npars = self.configurator.GetNPars()
        truth_selection = self.configurator.GetTruthSelection()
        n_truth = npars['n_truth']
        l = vecs.shape[0]
        self.data['Nobj'][j] = l
        self.data['Pmu'][j,:l,:] = vecs
        try: n_truth_local = len(truth_particles[j])
        except: n_truth_local = 0
        self.data['truth_Nobj'][j] = n_truth_local

        # Filling the truth_Pmu array.
        # For background generation, we may want to set "truth_selection" to None,
        # but still fill some empty numpy arrays so that the data format matches
        # some signal sample (if the arrays have different shapes, we may not be
        # able to simply concatenate signal and background HDF5 files).
        if(truth_selection is None):
            self.data['truth_Pmu'][j,:,:] = np.zeros((n_truth,4))
            self.data['truth_Pdg'][j,:] = np.zeros(n_truth)
        else:
            # Some truth-level selections may give a variable number of particles,
            # the buffer array is of a fixed length so there will be zero-padding.
            self.data['truth_Pmu'][j,:n_truth_local,:] = [
                [par.momentum.e, par.momentum.px, par.momentum.py, par.momentum.pz]
                for par in truth_particles[j]
            ]
            self.data['truth_Pdg'][j,:n_truth_local] = [par.pid for par in truth_particles[j]]

        if(jet is None):
            self.data['jet_Pmu'][j,:] = 0.
            self.data['jet_Pmu_cyl'][j,:] = 0.
        else:
            self.data['jet_Pmu'][j,:]     = [jet.e(), jet.px(), jet.py(), jet.pz()]
            self.data['jet_Pmu_cyl'][j,:] = [jet.pt(), jet.eta(), self.calculator.AdjustPhi(jet.phi()), jet.m()]

            if(self.diagnostic_plots):
                # Histogram the dR between the jet and each of the truth-level particles.
                for i, pid in enumerate(self.data['truth_Pdg'][j,:]):
                    if(pid not in self.jet_truth_dr_hists.keys()):
                        self.jet_truth_dr_hists[pid] = rt.TH1D('jet_truth_dr_{}'.format(pid),'',200,0.,2.)
                        title = '#Delta R #left( jet, {}_{} #right)'.format(pdg_names[pid],'{truth}')
                        self.jet_truth_dr_hists[pid].SetTitle(title)
                        self.jet_truth_dr_hists[pid].GetXaxis().SetTitle('#Delta R')
                        self.jet_truth_dr_hists[pid].GetYaxis().SetTitle('Count')
                    h = self.jet_truth_dr_hists[pid]
                    dr = np.sqrt(self.calculator.DeltaR2(truth_particles[j][i].momentum.eta(),truth_particles[j][i].momentum.phi(), self.data['jet_Pmu_cyl'][j,1],self.data['jet_Pmu_cyl'][j,2]))
                    h.Fill(dr)

        if(self.separate_truth_particles):
            n_separate = self.n_separate_truth_particles
            if(n_separate <= 0): n_separate = n_truth
            for k in range(n_separate):
                if(truth_selection is None): # TODO: Is this redundant? Should be initialized to zeros.
                    self.data['truth_Pmu_{}'.format(k)][j,:] = np.zeros(4)
                else:
                    if(k >= n_truth_local): break
                    # print('k = {}, len(truth_particles[{}]) = {}'.format(k,j,len(truth_particles[j])))
                    self.data['truth_Pmu_{}'.format(k)][j,:] = [truth_particles[j][k].momentum.e,truth_particles[j][k].momentum.px,truth_particles[j][k].momentum.py,truth_particles[j][k].momentum.pz]


        if(self.record_final_state_indices and final_state_indices is not None):
            self.data['final_state_idx'][j,:] = -1
            self.data['final_state_idx'][j,:l] = final_state_indices

        if(self.record_full_final_state):
            # We want to get all final-state vectors going into jet-clustering.

            if(self.delphes):
                vecs = np.vstack((
                    self.stable_particle_components['pt'],
                    self.stable_particle_components['eta'],
                    self.stable_particle_components['phi'],
                    self.stable_particle_components['m'],
                )).T
                vecs = self.calculator.PtEtaPhiMToEPxPyPz(vecs[:,0], vecs[:,1],vecs[:,2],vecs[:,3])

            else:
                vecs = np.vstack((
                    self.stable_particle_components['e'],
                    self.stable_particle_components['px'],
                    self.stable_particle_components['py'],
                    self.stable_particle_components['pz'],
                )).T
            key = 'Pmu_full_final_state'
            n = int(np.minimum(vecs.shape[0],self.data[key].shape[1]))
            self.data['Pmu_full_final_state'][j,:n,:] = vecs[:n]

        # Let's make some more diagnostic plots, to look at the jet kinematics, and the kinematics of the sum of jet constituents.
        if(self.diagnostic_plots):

            # Kinematic histograms for sum of Pmu (jet constituents).
            s = np.sum(vecs,axis=0) # E, px, py, pz
            s_cyl = self.calculator.EPxPyPzToPtEtaPhiM_single(*s)
            s_pt,s_eta,s_phi,s_m = s_cyl
            s_e = s[0]
            y = self.calculator.EPzToRap(s[0],s[3])
            s_et = s[0]/np.cosh(y)
            self.pmu_sum_hists['pt'].Fill(s_pt)
            self.pmu_sum_hists['eta'].Fill(s_eta)
            self.pmu_sum_hists['phi'].Fill(s_phi)
            self.pmu_sum_hists['m'].Fill(s_m)
            self.pmu_sum_hists['e'].Fill(s_e)
            self.pmu_sum_hists['et'].Fill(s_et)

            if(jet is not None):

                # Kinematic histograms for jets.
                jet_pt = self.data['jet_Pmu_cyl'][j,0]
                jet_eta = self.data['jet_Pmu_cyl'][j,1]
                jet_phi = self.data['jet_Pmu_cyl'][j,2]
                jet_m = self.data['jet_Pmu_cyl'][j,3]
                jet_e = jet.e()
                jet_et = jet.Et()

                self.jet_hists['pt'].Fill(jet_pt)
                self.jet_hists['eta'].Fill(jet_eta)
                self.jet_hists['phi'].Fill(jet_phi)
                self.jet_hists['m'].Fill(jet_m)
                self.jet_hists['e'].Fill(jet_e)
                self.jet_hists['et'].Fill(jet_et)

                # Kinematic difference histograms.
                delta_e = (s_e - jet_e) / jet_e
                delta_m = (s_m - jet_m) / jet_m
                delta_pt = (s_pt - jet_pt) / jet_pt
                delta_et = (s_et - jet_et) / jet_et

                self.jet_diff_hists['de'].Fill(delta_e)
                self.jet_diff_hists['dm'].Fill(delta_m)
                self.jet_diff_hists['dpt'].Fill(delta_pt)
                self.jet_diff_hists['det'].Fill(delta_et)

                self.jet_diff_hists_2d['de'].Fill(jet_e,delta_e)
                self.jet_diff_hists_2d['dm'].Fill(jet_m,delta_m)
                self.jet_diff_hists_2d['dpt'].Fill(jet_pt,delta_pt)
                self.jet_diff_hists_2d['det'].Fill(jet_et,delta_et)
        return

    def PostProcess(self,hepmc_files, h5_files=None, indices_files=None):
        if(self.post_processing is None):
            return
        nfiles = len(hepmc_files)
        # if(truth_files is not None): assert(nfiles == len(truth_files)) # TODO: may have to rework some of this logic
        if(h5_files is not None): assert(nfiles == len(h5_files))
        if(indices_files is not None): assert(nfiles == len(indices_files))

        for post_proc in self.post_processing:
            if(post_proc is None): continue
            post_proc.SetConfigurator(self.configurator)
            for i in range(nfiles):
                # truth_file = None
                h5_file = None
                idx_file = None
                fs_file = '{}/{}'.format(self.outdir,hepmc_files[i])
                # if(truth_files is not None): truth_file = '{}/{}'.format(self.outdir,truth_files[i])
                if(h5_files is not None): h5_file = '{}/{}'.format(self.outdir,h5_files[i])
                if(indices_files is not None): idx_file = '{}/{}'.format(self.outdir,indices_files[i])
                post_proc(fs_file,h5_file,idx_file)
        return

    def MergeEventFilterFlag(self,h5_file,event_filter_flag_files,copts=9):
        if(type(event_filter_flag_files) not in [list,tuple]):
            event_filter_flag_files = [event_filter_flag_files]
        f = h5.File('{}/{}'.format(self.outdir,h5_file),'a')
        keys = list(f.keys())

        try:
            g = [h5.File('{}/{}'.format(self.outdir,x),'r') for x in event_filter_flag_files]
        except:
            pass
            f.close()
            return

        gkeys = list(g[0].keys())

        keys_to_merge = []
        for key in gkeys:
            if(key in keys):
                print('\tWarning: Found key {} among the event_filter_flag keys, but this matches an existing key in the dataset. Skipping.'.format(key))
                continue
            keys_to_merge.append(key)

        for key in keys_to_merge:
            data = np.concatenate([x[key][:] for x in g],axis=0)
            f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)

        f.close()
        for x in g:
            x.close()
        return

    def ExtractHepMCParticles(self,events,n_par):
        particles = [
            ev.particles[:int(np.minimum(len(ev.particles),n_par))]
            for ev in events
        ]
        return particles

    def PrepDelphesArrays(self,):
        # TODO: This uses uproot.lazy, which is part of uproot 4 but has been
        #       deprecated and removed in uproot 5, in favor of uproot.dask.
        #       We'll probably need to support both versions soon.

        types = self.configurator.GetDelphesObjects()
        components = ['PT','Eta','Phi','ET'] # not all types have all components, this is OK
        delphes_keys = ['{x}.{y}'.format(x=x,y=y) for x in types for y in components]
        delphes_tree = 'Delphes'
        delphes_trees = ['{}:{}'.format('{}/{}'.format(self.outdir,x),delphes_tree) for x in self.delphes_files]
        delphes_arr = ur.lazy(delphes_trees, full_paths=False, filter_branch=lambda b: b.name in delphes_keys)
        delphes_keys = delphes_arr.fields # only keep branches that actually exist

        # Some convenient "mapping" between branch names and variable names we'll use.
        var_map = {key:{} for key in types}
        for branch in delphes_keys:
            key,var = branch.split('.')
            if('Eta' in var): var_map[key]['eta'] = branch
            elif('Phi' in var): var_map[key]['phi'] = branch
            else: var_map[key]['pt'] = branch
        return delphes_arr,var_map

    # def PrepDataBuffer(self,nentries_per_chunk,separate_truth_particles=False, n_separate=-1,final_state_indices=False, full_final_state=False):
    #     nentries_per_chunk = int(nentries_per_chunk)
    #     npars = self.configurator.GetNPars()
    #     n_constituents = npars['jet_n_par']
    #     n_truth = npars['n_truth']

    #     # Create our Numpy buffers to hold data. This dictates the structure of the HDF5 file we're making,
    #     # each key here corresponds to an HDF5 dataset. The shapes of the datasets will be what is shown here,
    #     # except with nentries_per_chunk -> nentries.
    #     self.data = {
    #         'Nobj'        : np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of jet constituents
    #         'Pmu'         : np.zeros((nentries_per_chunk,n_constituents,4),dtype=np.dtype('f8')), # 4-momenta of jet constituents (E,px,py,pz)
    #         'truth_Nobj'  : np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of truth-level particles (this is somewhat redundant -- it will typically be constant)
    #         'truth_Pdg'   : np.zeros((nentries_per_chunk,n_truth),dtype=np.dtype('i4')), # PDG codes to ID truth particles
    #         'truth_Pmu'   : np.zeros((nentries_per_chunk,n_truth,4),dtype=np.dtype('f8')), # truth-level particle 4-momenta
    #         'is_signal'   : np.zeros(nentries_per_chunk,dtype=np.dtype('i1')), # signal flag (0 = background, 1 = signal)
    #         'jet_Pmu'     : np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')), # jet 4-momentum, in Cartesian coordinates (E, px, py, pz)
    #         'jet_Pmu_cyl' : np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')) # jet 4-momentum, in cylindrical coordinates (pt,eta,phi,m)
    #     }

    #     # In addition to the above, we will also store the truth-level 4-momenta with each in its own HDF5 dataset.
    #     # In practice, this might be a more convenient storage format for things like neural network training.
    #     # Note, however, that the number of these keys will depend on n_truth (the number of truth particles saved per event).
    #     # TODO: Would be nice to accomplish this with references if possible, see: https://docs.h5py.org/en/stable/refs.html#refs .
    #     #       Not clear if that is actually feasible.
    #     if(separate_truth_particles):
    #         if(n_separate <= 0): n_separate = n_truth
    #         for i in range(n_separate):
    #             key = 'truth_Pmu_{}'.format(i)
    #             self.data[key] = np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8'))

    #     # We can optionally record the index of each Pmu with respect to all the four-momenta that were passed to jet clustering.
    #     # This may be useful if we are trying to somehow trace certain information through our data generation pipeline,
    #     # e.g. if we're trying to keep track of whether daughter particles of a particular decay are in a given jet.
    #     # We'll keep this optional since, if we're not doing something like that, this extra info may be useless or even confusing.
    #     # This uses zero-indexing, so we'll fill things with -1 to avoid any confusion.
    #     #TODO: The filling with the -1's is done when filling the buffer with data. Otherwise the -1's get changed to 0's, I can't remember how/why.
    #     if(final_state_indices):
    #         key = 'final_state_idx'
    #         self.data[key] = np.zeros((nentries_per_chunk,n_constituents),dtype=np.dtype('i4'))

    #     if(full_final_state): # This is really only for debugging purposes! We truncate to a fixed length.
    #         key = 'Pmu_full_final_state'
    #         n_max = 400
    #         self.data[key] = np.zeros((nentries_per_chunk,n_max,4),dtype=np.dtype('f8'))
    #     return

    def PrepH5File(self,filename,nentries,data_buffer):
        dsets = {}
        with h5.File(filename, 'w') as f:
            for key, val in data_buffer.items():
                shape = list(val.shape)
                shape[0] = nentries
                shape = tuple(shape)
                dsets[key] = f.create_dataset(key, shape, val.dtype,compression='gzip')
        return dsets

    def PrepIndexRanges(self,nentries,nentries_per_chunk):
        nchunks = int(np.ceil(nentries / nentries_per_chunk))
        start_idxs = np.zeros(nchunks,dtype = np.dtype('i8'))
        for i in range(1,start_idxs.shape[0]): start_idxs[i] = start_idxs[i-1] + nentries_per_chunk
        stop_idxs = start_idxs + nentries_per_chunk
        stop_idxs[-1] = nentries
        ranges = stop_idxs - start_idxs
        return start_idxs,stop_idxs,ranges

# Generic function for concatenating HDF5 files.
def ConcatenateH5(input_files,output_file,cwd=None,delete_inputs=False, compression='gzip', copts=9,ignore_keys=None,verbose=False,silent_drop=False):
    if(cwd is not None):
        input_files = ['{}/{}'.format(cwd,x) for x in input_files]
    infiles = [h5.File(f,'r') for f in input_files]
    keys = list(infiles[0].keys())
    if(cwd is not None):
        output_file = '{}/{}'.format(cwd,output_file)
    outfile = h5.File(output_file,'w')

    if(verbose):
        print('\n\t' + 21*"#")
        print("\t### ConcatenateH5 ###")
        for i,infile in enumerate(input_files):
            print('\t  Input {}:\t{}'.format(i,infile))
        print('\tOutput:\t{}'.format(output_file))

    if(ignore_keys is not None):
        keys = [k for k in keys if k not in ignore_keys]
        if(verbose and not silent_drop):
            print('\tExcluding the following keys from concatenation, these will be dropped:')
            for k in ignore_keys:
                print('\t\t{}'.format(k))

    for key in keys:
        data = np.concatenate([f[key] for f in infiles],axis=0)
        outfile.create_dataset(key,data=data,compression=compression,compression_opts=copts)

    outfile.close()

    for f in infiles:
        f.close()

    if(verbose):
        if(delete_inputs):
            print('\tDeleting input files.')
        print('\t' + 21*"#" + '\n')

    if(delete_inputs):
        for f in input_files:
            sub.check_call(['rm',f])

def MergeH5(target_file, input_file, cwd=None, delete_stats_file=False, compression='gzip',copts=9):
    if(cwd is not None):
        target_file = '{}/{}'.format(cwd,target_file)
        input_file = '{}/{}'.format(cwd,input_file)

    f = h5.File(target_file,'a')
    g = h5.File(input_file,'r')
    keys = list(g.keys())
    f_keys = list(f.keys())

    for key in keys:
        if key in f_keys:
            print('Warning: key {} found in {} before merging in {}.'.format(key,target_file,input_file))
            continue
        f.create_dataset(key,data=g[key][:],compression=compression,compression_opts=copts)

    f.close()
    g.close()
    if(delete_stats_file):
        sub.check_call(['rm',input_file])
    return

def RemoveFailedFromHDF5(h5_file, cwd=None, copts=9):
    """
    Remove any failed events from the HDF5 file -- these are identified
    as those with a negative value in the "is_signal" dataset.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

    fname_tmp = h5_file.replace('.h5','_{}.h5'.format(str(uuid.uuid4())))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    keep_indices = f['is_signal'][:] >= 0
    if(np.sum(keep_indices) == N): return

    g = h5.File(fname_tmp,'w')

    for i,key in enumerate(keys):
        g.create_dataset(key, data = f[key][:][keep_indices],compression='gzip',compression_opts=copts)

    f.close()
    g.close()

    sub.check_call(['rm',h5_file])
    sub.check_call(['mv',fname_tmp, h5_file])
    return

# Add a column with event indices.
def AddEventIndices(h5_file,cwd=None,copts=9,key='event_idx',offset=0):
    """
    Add a dataset with event indices to an HDF5 file. This is a sequential index,
    which may be useful if some entries will later be dropped and one wants to
    keep track of which events these were.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['is_signal'].shape[0]
    event_indices = np.arange(offset,nevents+offset,dtype=int)
    f.create_dataset(key,data=event_indices, compression='gzip', compression_opts=copts)
    f.close()

def AddConstantValue(h5_file,cwd=None,copts=9,value=0,key='constant_value',dtype=None):
    """
    Add a dataset with some constant value to the HDF5 file. In practice this can
    be used to attach metadata to events -- like the Pythia8 RNG seed that was used.
    This may become useful when multiple datasets are combined.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['is_signal'].shape[0]

    if(dtype is None): dtype = np.dtype('f8')
    if(dtype != str):
        data = np.full((nevents),value)
    else:
        data = nevents * [value]
    f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
    f.close()

def AddMetaData(h5_file,cwd=None,value='',key='metadata'):
    """
    Generic function for adding a value to the HDF5 file
    metadata container.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    f.attrs[key] = value
    f.close()

def AddMetaDataWithReference(h5_file,cwd=None,value='',key='metadata',copts=9):
    """
    Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
    Also creates a column in the dataset that will point to this metadata's index.
    Somewhat redundant for file generation but this type of logic will be useful when concatenating files
    with different entries in the metadata fields.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['is_signal'].shape[0]
    metadata = f.attrs
    if(key not in metadata.keys()): f.attrs[key] = [value]
    else: f.attrs[key] = list(f.attrs[key]) + [value] # I think the list <-> array stuff should be OK here
    idx = len(f.attrs[key]) - 1
    f.create_dataset(key,data=np.full(nevents,idx,dtype=np.dtype('i4')),compression='gzip',compression_opts=copts)
    f.close()

def SplitH5(h5_file, split_ratio = (7,2,1), train_name=None, val_name=None, test_name=None, cwd=None, copts=9,verbose=False, seed=0):
    """
    Split an HDF5 file into training, testing and validation samples.
    The split is random (with the RNG seed provided).
    Note that this copies over the full metadata, though with the splitting
    it is technically possible that one of the files has some metadata
    entries with no corresponding events.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    file_dir = '/'.join(h5_file.split('/')[:-1])

    if(train_name is None):
        train_name = 'train.h5'
    if(val_name is None):
        val_name   =   'valid.h5'
    if(test_name is None):
        test_name  =  'test.h5'

    if(cwd is not None):
        train_name = '{}/{}'.format(cwd,train_name)
        val_name = '{}/{}'.format(cwd,val_name)
        test_name = '{}/{}'.format(cwd,test_name)

    if(verbose):
        print("\n\t" + 15*"#")
        print('\t### SplitH5 ###')
        print('\tInput: {}'.format(h5_file))
        print('\tOutputs:')
        for (name,frac) in zip((train_name,val_name,test_name),split_ratio):
            print('\t\t {} \t (fraction of events = {:.2e})'.format(name,frac))
        print('\tSplitting using seed: {}'.format(seed))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    split_ratio = np.array(split_ratio)
    tot = np.sum(np.array(split_ratio))
    split_ratio = split_ratio / tot
    split_ratio[-1] = 1. - np.sum(split_ratio[:-1])

    names = [train_name, val_name, test_name]
    n = np.array(N * split_ratio,dtype=int)
    diff = np.sum(n) - N
    n[-1] -= diff

    rng = np.random.default_rng(seed)
    index_list = np.array(range(N),dtype=int)
    rng.shuffle(index_list,axis=0)
    indices = []
    start = 0
    stop = 0
    for i in range(len(names)):
        stop += n[i]
        idxs = index_list[start:stop]
        start = stop
        indices.append(idxs)

    for i in range(len(names)):
        g = h5.File(names[i],'w')
        for j,key in enumerate(keys):
                # Putting the data from f into memory (numpy array)
                # since h5py doesn't like non-ordered indices.
                g.create_dataset(key, data=f[key][:][indices[i]],compression='gzip', compression_opts=copts)
        for key in list(f.attrs.keys()):
            g.attrs[key] = f.attrs[key]
        g.close()
    f.close()

    if(verbose):
        print("\t" + 15*"#" + "\n")
    return