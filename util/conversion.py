import sys, glob, uuid
import numpy as np
import h5py as h5
import subprocess as sub
from util.config import GetNPars, GetJetConfig, GetInvisiblesFlag, GetSignalFlag
from util.calcs import PtEtaPhiMToPxPyPzE, PtEtaPhiMToEPxPyPz, AdjustPhi
from util.fastjet import BuildFastjet
from util.qol_utils.qol_util import printProgressBarColor
from util.conv_utils.utils import ExtractHepMCParticles, InitFastJet, PrepDataBuffer, ExtractHepMCEvents,PrepDelphesArrays, PrepH5File, PrepIndexRanges
from util.particle_selection import IsNeutrino
# --- FASTJET IMPORT ---
# TODO: Can this be done more nicely?
fastjet_dir = BuildFastjet(j=8)
fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
import fastjet as fj
# ----------------------

# Convert HepMC / Delphes-ROOT files into our final output -- possibly performing jet clustering along the way.
class Processor:
    def __init__(self, use_delphes=False):
        self.delphes = use_delphes

        self.prefix_level1 = '\tClustering jets & preparing data:'
        self.suffix = 'Complete'
        self.bl = 50

        self.jet_config, self.jetdef = InitFastJet()

    def Process(self, final_state_files, truth_files=None, h5_file=None, nentries_per_chunk=1e4, verbosity=0, separate_truth_particles=True):

        if(h5_file is None):
            if(type(final_state_files) == list): h5_file = final_state_files[0]
            else: h5_file = final_state_files
            if(self.delphes): h5_file =  h5_file.replace('*','').replace('.root','.h5')
            else: h5_file =  h5_file.replace('*','').replace('.root','.h5')
        npars = GetNPars()
        n_truth = npars['n_truth']

            # Note: It's important that the Delphes files and truth HepMC files line up!
        #       The way we usually do this, it will be guaranteed.
        if(type(final_state_files) == str): final_state_files = glob.glob(final_state_files,recursive=True)
        if(truth_files is None):
            if(self.delphes): truth_files = [x.replace('.root','_truth.hepmc') for x in final_state_files]
            else: truth_files = [x.replace('.hepmc','_truth.hepmc') for x in final_state_files]

        if(self.delphes):
            delphes_arr,var_map = PrepDelphesArrays(final_state_files)
            nentries = len(delphes_arr)
        else:
            ## Extract final-state truth particle info from the HepMC files.
            final_state_events, nentries = ExtractHepMCEvents(final_state_files,get_nevents=True)

        # Also extract truth particle info from the HepMC files.
        truth_events = ExtractHepMCEvents(truth_files)

        nentries_per_chunk = int(nentries_per_chunk)
        data = PrepDataBuffer(nentries_per_chunk,separate_truth_particles)

        # Prepare the HDF5 file.
        dsets = PrepH5File(h5_file,nentries,data)

        # Some indexing preparation for writing in chunks.
        start_idxs,stop_idxs,ranges = PrepIndexRanges(nentries,nentries_per_chunk)
        nchunks = len(ranges)

        if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)

        for i in range(nchunks):
            # Clear the buffer (for safety).
            for key in data.keys(): data[key][:] = 0

            if(self.delphes):
                # Get the momenta of the final-state particles, first in cylindrical coordinates.
                # We will be setting m=0 for all these particles.
                momenta = {}
                types = var_map.keys() # TODO: Can this be placed elsewhere? Renamed?
                for delphes_type in types:
                    momenta[delphes_type] = {
                        x:delphes_arr[var_map[delphes_type][x]][start_idxs[i]:stop_idxs[i]]
                        for x in ['pt','eta','phi']
                    }
            else:
                ## Get final-state particles.
                npar_fs = 1e4 # TODO: This is some hardcoded max number of final-state particles. Should be plenty.
                final_state_particles = ExtractHepMCParticles(final_state_events[start_idxs[i]:stop_idxs[i]],npar_fs)

            # Get the truth particles.
            truth_particles = ExtractHepMCParticles(truth_events[start_idxs[i]:stop_idxs[i]],n_truth)

            prefix_nzero = int(np.ceil(np.log10(nchunks))) + 1
            prefix_level2 = '\tClustering jets & preparing data for chunk {}/{}:'.format(str(i+1).zfill(prefix_nzero),nchunks)
            if(verbosity==2): printProgressBarColor(0,ranges[i], prefix=prefix_level2, suffix=self.suffix, length=self.bl)

            # For each event we must combine tracks and neutral hadrons, perform jet clustering on them,
            # select a single jet (based on user criteria), and select that jet's leading constituents.
            for j in range(ranges[i]):
                data['is_signal'][j] = GetSignalFlag()

                if(self.delphes):
                    pt  = np.concatenate([momenta[x]['pt'][j].to_numpy() for x in types])
                    eta = np.concatenate([momenta[x]['eta'][j].to_numpy() for x in types])
                    phi = np.concatenate([momenta[x]['phi'][j].to_numpy() for x in types])
                    m   = np.zeros(pt.shape)

                    vecs = PtEtaPhiMToPxPyPzE(pt,eta,phi,m) # convert 4-vectors to (px, py, pz, e) for jet clustering.
                    px = vecs[:,0]
                    py = vecs[:,1]
                    pz = vecs[:,2]
                    e  = vecs[:,3]

                else:
                    l = len(final_state_particles[j])

                    # Optionally exclude invisibles (neutrinos). If using Delphes, this should be handled by Delphes internally.
                    if(GetInvisiblesFlag):
                        visibles = np.arange(l)
                    else:
                        invisibles = np.array([IsNeutrino(final_state_particles[j][k]) for k in range(l)])
                        visibles = np.where(invisibles == 0)[0]
                        del invisibles

                    px = np.array([final_state_particles[j][k].momentum.px for k in visibles])
                    py = np.array([final_state_particles[j][k].momentum.py for k in visibles])
                    pz = np.array([final_state_particles[j][k].momentum.pz for k in visibles])
                    e  = np.array([final_state_particles[j][k].momentum.e  for k in visibles])

                self.SelectFinalStateParticles(px,py,pz,e, self.jetdef, truth_particles, data, j, separate_truth_particles)

                if(verbosity==2): printProgressBarColor(j+1,ranges[i], prefix=prefix_level2, suffix=self.suffix, length=self.bl)

            # We have now filled a chunk, time to write it.
            with h5.File(h5_file, 'a') as f:
                for key in dsets.keys():
                    dset = f[key]
                    dset[start_idxs[i]:stop_idxs[i]] = data[key][:ranges[i]]
            if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=self.prefix_level1, suffix=self.suffix, length=self.bl)
        return h5_file


    # --- Utility functions for the class ---
    def InitFastJet(self):
        fj.ClusterSequence.print_banner() # Get the Fastjet banner out of the way.
        jet_config = GetJetConfig()
        jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
        return jet_config,jetdef

    def ClusterJets(self,vecs, jetdef, jet_config):
        pj = [fj.PseudoJet(*x) for x in vecs]
        jets = jetdef(pj)

        # Apply optional minimum jet pT cut.
        jet_pt = np.array([jet.pt() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[jet_pt >= jet_config['jet_min_pt']]
        jets = [jets[i] for i in jet_indices]

        # Apply optional maximum |eta| cut.
        jet_eta = np.array([jet.eta() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[np.abs(jet_eta) <= jet_config['jet_max_eta']]
        jets = [jets[i] for i in jet_indices]
        return jets

    def FetchJetConstituents(self,jet,n_constituents,zero_mass=False):
        if(zero_mass):
            pt,eta,phi,m = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), 0.] for x in jet.constituents()]),4)
        else:
            pt,eta,phi,m = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m()] for x in jet.constituents()]),4)

        pt = pt.flatten()
        eta = eta.flatten()
        phi = phi.flatten()
        m = m.flatten()

        # Sort by decreasing pt, and only keep leading constituents.
        sorting = np.argsort(-pt)
        l = int(np.minimum(n_constituents,len(pt)))
        vec = PtEtaPhiMToEPxPyPz(pt=pt[sorting][:l], # Vector will be of format (E, px, py, pz)
                                    eta=eta[sorting][:l],
                                    phi=phi[sorting][:l],
                                    m=m[sorting][:1]
                                )
        return vec

    def SelectFinalStateParticles(self,px,py,pz,e, jetdef, truth_particles, data, j, separate_truth_particles):
        jet_config = GetJetConfig()
        n_constituents = GetNPars()['jet_n_par']
        jet_sel = jet_config['jet_selection']

        if(jet_sel is None):
            vecs = np.vstack((e,px,py,pz)).T # note energy is given as 0th component
            self.FillDataBuffer(data,j,vecs,None,truth_particles,separate_truth_particles)
        else:
            vecs = np.vstack((px,py,pz,e)).T # note energy is given as last component, needed for fastjet
            jets = self.ClusterJets(vecs,jetdef,jet_config)

            njets = len(jets)
            if(njets == 0):
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return

            # selected_jet_idx = jet_config['jet_selection'](truth=truth_particles[j], jets=jets, use_hepmc=True)
            selected_jet_idx = jet_sel(truth=truth_particles[j], jets=jets)
            if(selected_jet_idx < 0):
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                return
            jet = jets[selected_jet_idx]
            # Get the constituents of our selected jet. Assuming they are massless -> don't need to fetch the mass.
            vecs = self.FetchJetConstituents(jet,n_constituents)
            self.FillDataBuffer(data,j,vecs,jet,truth_particles,separate_truth_particles)
        return

    def FillDataBuffer(self,data_buffer,j,vecs,jet,truth_particles,separate_truth_particles=False):
        npars = GetNPars()
        n_truth = npars['n_truth']
        l = vecs.shape[0]
        data_buffer['Nobj'][j] = l
        data_buffer['Pmu'][j,:l,:] = vecs
        data_buffer['truth_Nobj'][j] = len(truth_particles[j])
        data_buffer['truth_Pdg'][j,:] = [par.pid for par in truth_particles[j]]
        data_buffer['truth_Pmu'][j,:,:] = [
            [par.momentum.e, par.momentum.px, par.momentum.py, par.momentum.pz]
            for par in truth_particles[j]
        ]

        if(jet is None):
            data_buffer['jet_Pmu'][j,:] = 0.
            data_buffer['jet_Pmu_cyl'][j,:] = 0.
        else:
            data_buffer['jet_Pmu'][j,:]     = [jet.e(), jet.px(), jet.py(), jet.pz()]
            data_buffer['jet_Pmu_cyl'][j,:] = [jet.pt(), jet.eta(), AdjustPhi(jet.phi()), jet.m()]

        if(separate_truth_particles):
            for k in range(n_truth):
                data_buffer['truth_Pmu_{}'.format(k)][j,:] = [truth_particles[j][k].momentum.e,truth_particles[j][k].momentum.px,truth_particles[j][k].momentum.py,truth_particles[j][k].momentum.pz]
        return

# Remove any "failed events". They will be marked with a negative signal flag.
def RemoveFailedFromHDF5(h5_file):
    fname_tmp = h5_file.replace('.h5','_{}.h5'.format(str(uuid.uuid4())))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    keep_indices = f['is_signal'][:] >= 0
    if(np.sum(keep_indices) == N): return

    g = h5.File(fname_tmp,'w')

    for i,key in enumerate(keys):
        g.create_dataset(key, data = f[key][:][keep_indices],compression='gzip')

    f.close()
    g.close()

    sub.check_call(['rm',h5_file])
    sub.check_call(['mv',fname_tmp, h5_file])
    return

# Split an HDF5 file into training, testing and validation samples.
def SplitHDF5(h5_file, split_ratio = (7,2,1), train_name=None, val_name=None, test_name=None):

    file_dir = '/'.join(h5_file.split('/')[:-1])

    if(train_name is None):
        train_name = 'train.h5'
        if(file_dir != ''): train_name = '{}/{}'.format(file_dir,train_name)
    if(val_name is None):
        val_name   =   'valid.h5'
        if(file_dir != ''): val_name = '{}/{}'.format(file_dir,val_name)
    if(test_name is None):
        test_name  =  'test.h5'
        if(file_dir != ''): test_name = '{}/{}'.format(file_dir,test_name)

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

    rng = np.random.default_rng()
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
                g.create_dataset(key, data=f[key][:][indices[i]],compression='gzip')
        g.close()
    f.close()
    return