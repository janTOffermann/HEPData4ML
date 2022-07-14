import sys, os, glob, uuid
import numpy as np
import pyhepmc_ng as hep
import uproot as ur
import h5py as h5
import subprocess as sub
from util.fastjet import BuildFastjet, ParticleInfo
from util.config import GetNPars, GetJetConfig
from util.calcs import PtEtaPhiMToPxPyPzE, PtEtaPhiMToEPxPyPz, AdjustPhi
from util.qol_util import printProgressBarColor

# Convert a set of Delphes ROOT files (modified with truth info!) files to a single HDF5 file.
# This also performs jet clustering of the Delphes output (and associated selections).
def DelphesWithTruthToHDF5(delphes_files, truth_files=None, h5_file=None, nentries_per_chunk=1e4, verbosity=0, separate_truth_particles=True):

    # ----- FASTJET SETUP -----
    fastjet_dir = BuildFastjet(j=8)
    fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
    if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
    import fastjet as fj

    fj.ClusterSequence.print_banner() # Get the Fastjet banner out of the way. TODO: Can we get rid of the banner? It is annoying.
    jet_config = GetJetConfig()
    jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
    # ----- END FASTJET SETUP -----

    if(h5_file is None):
        if(type(delphes_files) == list): h5_file = delphes_files[0]
        else: h5_file = delphes_files
        h5_file =  h5_file.replace('*','').replace('.root','.h5')
    npars = GetNPars()
    n_constituents = npars['jet_n_par']
    n_truth = npars['n_truth']

    # Note: It's important that the Delphes files and truth HepMC files line up!
    #       The way we usually do this, it will be guaranteed.
    if(type(delphes_files) == str): delphes_files = glob.glob(delphes_files,recursive=True)
    if(truth_files is None): truth_files = [x.replace('.root','_truth.hepmc') for x in delphes_files]

    # We are interested in extracting information on particles in the Delphes files.
    types = ['EFlowTrack','EFlowNeutralHadron']
    components = ['PT','Eta','Phi','ET'] # not all types have all components, this is OK
    delphes_keys = ['{x}.{y}'.format(x=x,y=y) for x in types for y in components]
    delphes_tree = 'Delphes'
    delphes_trees = ['{}:{}'.format(x,delphes_tree) for x in delphes_files]
    delphes_arr = ur.lazy(delphes_trees, full_paths=False, filter_branch=lambda b: b.name in delphes_keys)
    delphes_keys = delphes_arr.fields # only keep branches that actually exist
    nentries = len(delphes_arr)

    # Also extract truth particle info from the HepMC files.
    truth_events = []
    for truth_file in truth_files:
        with hep.ReaderAscii(truth_file) as f:
            for evt in f:
                truth_events.append(evt)

    # Some convenient "mapping" between branch names and variable names we'll use.
    var_map = {key:{} for key in types}
    for branch in delphes_keys:
        key,var = branch.split('.')
        if('Eta' in var): var_map[key]['eta'] = branch
        elif('Phi' in var): var_map[key]['phi'] = branch
        else: var_map[key]['pt'] = branch

    # Now create our Numpy buffers to hold data. This dictates the structure of the HDF5 file we're making,
    # each key here corresponds to an HDF5 dataset. The shapes of the datasets will be what is shown here,
    # except with nentries_per_chunk -> nentries.
    nentries_per_chunk = int(nentries_per_chunk)
    data = {
        'Nobj':np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of jet constituents
        'Pmu':np.zeros((nentries_per_chunk,n_constituents,4),dtype=np.dtype('f8')), # 4-momenta of jet constituents (E,px,py,pz)
        'truth_Nobj':np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of truth-level particles (this is somewhat redundant -- it will typically be constant)
        'truth_Pdg':np.zeros((nentries_per_chunk,n_truth),dtype=np.dtype('i4')), # PDG codes to ID truth particles
        'truth_Pmu':np.zeros((nentries_per_chunk,n_truth,4),dtype=np.dtype('f8')), # truth-level particle 4-momenta
        'is_signal':np.zeros(nentries_per_chunk,dtype=np.dtype('i1')), # signal flag (0 = background, 1 = signal)
        'jet_Pmu':np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')), # jet 4-momentum, in Cartesian coordinates (E, px, py, pz)
        'jet_Pmu_cyl':np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')) # jet 4-momentum, in cylindrical coordinates (pt,eta,phi,m)
    }

    # In addition to the above, we will also store the truth-level 4-momenta with each in its own HDF5 dataset.
    # In practice, this might be a more convenient storage format for things like neural network training.
    # Note, however, that the number of these keys will depend on n_truth (the number of truth particles saved per event).
    # TODO: Would be nice to accomplish this with references if possible, see: https://docs.h5py.org/en/stable/refs.html#refs .
    #       Not clear if that is actually feasible.
    if(separate_truth_particles):
        for i in range(n_truth):
            key = 'truth_Pmu_{}'.format(i)
            data[key] = np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8'))

    # Prepare the HDF5 file.
    dsets = {}
    with h5.File(h5_file, 'w') as f:
        for key, val in data.items():
            shape = list(val.shape)
            shape[0] = nentries
            shape = tuple(shape)
            dsets[key] = f.create_dataset(key, shape, val.dtype,compression='gzip')

    # Some indexing preparation for writing in chunks.
    nchunks = int(np.ceil(nentries / nentries_per_chunk))
    start_idxs = np.zeros(nchunks,dtype = np.dtype('i8'))
    for i in range(1,start_idxs.shape[0]): start_idxs[i] = start_idxs[i-1] + nentries_per_chunk
    stop_idxs = start_idxs + nentries_per_chunk
    stop_idxs[-1] = nentries
    ranges = stop_idxs - start_idxs

    # Some setup for (optional printouts).
    prefix_level1 = '\tClustering jets & preparing data:'
    suffix = 'Complete'
    bl = 50

    if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=prefix_level1, suffix=suffix, length=bl)

    for i in range(nchunks):
        # Clear the buffer (for safety).
        for key in data.keys(): data[key][:] = 0

        # Get the momenta of the final-state particles, first in cylindrical coordinates.
        # We will be setting m=0 for all these particles.
        momenta = {}

        for delphes_type in types:
            momenta[delphes_type] = {
                x:delphes_arr[var_map[delphes_type][x]][start_idxs[i]:stop_idxs[i]]
                for x in ['pt','eta','phi']
            }

        # Get the truth particles.
        truth_particles = [
            ev.particles[:int(np.minimum(len(ev.particles),n_truth))]
            for ev in truth_events[start_idxs[i]:stop_idxs[i]]
        ]

        prefix_nzero = int(np.ceil(np.log10(nchunks))) + 1
        prefix_level2 = '\tClustering jets & preparing data for chunk {}/{}:'.format(str(i+1).zfill(prefix_nzero),nchunks)
        if(verbosity==2): printProgressBarColor(0,ranges[i], prefix=prefix_level2, suffix=suffix, length=bl)

        # For each event we must combine tracks and neutral hadrons, perform jet clustering on them,
        # select a single jet (based on user criteria), and select that jet's leading constituents.
        for j in range(ranges[i]):
            data['is_signal'][j] = 1 # TODO: redundant for now
            pt  = np.concatenate([momenta[x]['pt'][j].to_numpy() for x in types])
            eta = np.concatenate([momenta[x]['eta'][j].to_numpy() for x in types])
            phi = np.concatenate([momenta[x]['phi'][j].to_numpy() for x in types])
            m   = np.zeros(pt.shape)

            vecs = PtEtaPhiMToPxPyPzE(pt,eta,phi,m) # convert 4-vectors to (px, py, pz, e) for jet clustering.


            ## We already have the px, py, pz, e, and we do not need to adjust it, just extract the information directly from the hempc file.


            # ----- Begin jet clustering -----
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

            # Check if there are any jets left. We may have executed a "jet check" before Delphes,
            # but after passing things through Delphes we may no longer have jets that meet our cuts.

            njets = len(jets)
            if(njets == 0):
                #TODO: Determine what's best to do here.
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                continue

            # Now select a jet.
            selected_jet_idx = jet_config['jet_selection'](truth = truth_particles[j], jets = jets, use_hepmc=True)
            if(selected_jet_idx < 0):
                # TODO: Determine what's best to do here.
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                continue
            jet = jets[selected_jet_idx]

            # Get the constituents of our selected jet. Assuming they are massless -> don't need to fetch the mass.
            pt,eta,phi = np.hsplit(np.array([[x.pt(), x.eta(), x.phi()] for x in jet.constituents()]),3)
            pt = pt.flatten()
            eta = eta.flatten()
            phi = phi.flatten()
            # ----- End jet clustering -----

            # Sort by decreasing pt, and only keep leading constituents.
            sorting = np.argsort(-pt)
            l = int(np.minimum(n_constituents,len(pt)))
            vec = PtEtaPhiMToEPxPyPz(pt=pt[sorting][:l], # Vector will be of format (E, px, py, pz)
                                     eta=eta[sorting][:l],
                                     phi=phi[sorting][:l],
                                     m=np.zeros(l)
                                    )

            data['Nobj'][j] = l
            data['Pmu'][j,:l,:] = vec[:l,:] # indexing for vec is redundant
            # Now take care of truth particles. We probably expect the same number of truth particles
            # per event, but to be safe we will always truncate the list to the requested number at most.
            data['truth_Nobj'][j] = len(truth_particles[j])
            data['truth_Pdg'][j,:] = [par.pid for par in truth_particles[j]]
            data['truth_Pmu'][j,:,:] = [
                [par.momentum.e, par.momentum.px, par.momentum.py, par.momentum.pz]
                for par in truth_particles[j]
            ]
            # Also fill our separate buffers for each truth particle, if requested.
            if(separate_truth_particles):
                for k in range(n_truth):
                    data['truth_Pmu_{}'.format(k)][j,:] = [truth_particles[j][k].momentum.e,truth_particles[j][k].momentum.px,truth_particles[j][k].momentum.py,truth_particles[j][k].momentum.pz]

            # Also store some jet-level information for debugging purposes.
            # We will also store the jet momentum in cylindrical coordinates (pt,eta,phi,m).
            # We make sure that the jet's phi is in [-pi ,pi] for consistency.
            data['jet_Pmu'][j,:]     = [jet.e(), jet.px(), jet.py(), jet.pz()]
            data['jet_Pmu_cyl'][j,:] = [jet.pt(), jet.eta(), AdjustPhi(jet.phi()), jet.m()]

            if(verbosity==2): printProgressBarColor(j+1,ranges[i], prefix=prefix_level2, suffix=suffix, length=bl)

        # We have now filled a chunk, time to write it.
        with h5.File(h5_file, 'a') as f:
            for key in dsets.keys():
                dset = f[key]
                dset[start_idxs[i]:stop_idxs[i]] = data[key][:ranges[i]]
        if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=prefix_level1, suffix=suffix, length=bl)
    return h5_file

# Convert a set of Delphes ROOT files (modified with truth info!) files to a single HDF5 file.
# This also performs jet clustering of the Delphes output (and associated selections).
def HepMCWithTruthToHDF5(final_state_files, truth_files=None, h5_file=None, nentries_per_chunk=1e4, verbosity=0, separate_truth_particles=True):

    # ----- FASTJET SETUP -----
    fastjet_dir = BuildFastjet(j=8)
    fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
    if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
    import fastjet as fj

    fj.ClusterSequence.print_banner() # Get the Fastjet banner out of the way. TODO: Can we get rid of the banner? It is annoying.
    jet_config = GetJetConfig()
    jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
    # ----- END FASTJET SETUP -----

    if(h5_file is None):
        if(type(final_state_files) == list): h5_file = final_state_files[0]
        else: h5_file = final_state_files
        h5_file =  h5_file.replace('*','').replace('.hepmc','.h5')
    npars = GetNPars()
    n_constituents = npars['jet_n_par']
    n_truth = npars['n_truth']

    # Note: It's important that the Delphes files and truth HepMC files line up!
    #       The way we usually do this, it will be guaranteed.
    if(type(final_state_files) == str): final_state_files = glob.glob(final_state_files,recursive=True)
    if(truth_files is None): truth_files = [x.replace('.hepmc','_truth.hepmc') for x in final_state_files]

    # nentries = len(delphes_arr)

    # Also extract truth particle info from the HepMC files.
    truth_events = []
    for truth_file in truth_files:
        with hep.ReaderAscii(truth_file) as f:
            for evt in f:
                truth_events.append(evt)

    ## Extract final-state truth particle info from the HepMC files.
    nentries = 0
    final_state_events = []
    for fs_file in final_state_files:
        with hep.ReaderAscii(fs_file) as f:
            for evt in f:
                final_state_events.append(evt)
                nentries += 1

    # Now create our Numpy buffers to hold data. This dictates the structure of the HDF5 file we're making,
    # each key here corresponds to an HDF5 dataset. The shapes of the datasets will be what is shown here,
    # except with nentries_per_chunk -> nentries.
    nentries_per_chunk = int(nentries_per_chunk)
    data = {
        'Nobj':np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of jet constituents
        'Pmu':np.zeros((nentries_per_chunk,n_constituents,4),dtype=np.dtype('f8')), # 4-momenta of jet constituents (E,px,py,pz)
        'truth_Nobj':np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of truth-level particles (this is somewhat redundant -- it will typically be constant)
        'truth_Pdg':np.zeros((nentries_per_chunk,n_truth),dtype=np.dtype('i4')), # PDG codes to ID truth particles
        'truth_Pmu':np.zeros((nentries_per_chunk,n_truth,4),dtype=np.dtype('f8')), # truth-level particle 4-momenta
        'is_signal':np.zeros(nentries_per_chunk,dtype=np.dtype('i1')), # signal flag (0 = background, 1 = signal)
        'jet_Pmu':np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')), # jet 4-momentum, in Cartesian coordinates (E, px, py, pz)
        'jet_Pmu_cyl':np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')) # jet 4-momentum, in cylindrical coordinates (pt,eta,phi,m)
    }

    # In addition to the above, we will also store the truth-level 4-momenta with each in its own HDF5 dataset.
    # In practice, this might be a more convenient storage format for things like neural network training.
    # Note, however, that the number of these keys will depend on n_truth (the number of truth particles saved per event).
    if(separate_truth_particles):
        for i in range(n_truth):
            key = 'truth_Pmu_{}'.format(i)
            data[key] = np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8'))

    # Prepare the HDF5 file.
    dsets = {}
    with h5.File(h5_file, 'w') as f:
        for key, val in data.items():
            shape = list(val.shape)
            shape[0] = nentries
            shape = tuple(shape)
            dsets[key] = f.create_dataset(key, shape, val.dtype,compression='gzip')

    # Some indexing preparation for writing in chunks.
    nchunks = int(np.ceil(nentries / nentries_per_chunk))
    start_idxs = np.zeros(nchunks,dtype = np.dtype('i8'))
    for i in range(1,start_idxs.shape[0]): start_idxs[i] = start_idxs[i-1] + nentries_per_chunk
    stop_idxs = start_idxs + nentries_per_chunk
    stop_idxs[-1] = nentries
    ranges = stop_idxs - start_idxs

    # Some setup for (optional printouts).
    prefix_level1 = '\tClustering jets & preparing data:'
    suffix = 'Complete'
    bl = 50

    if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=prefix_level1, suffix=suffix, length=bl)

    for i in range(nchunks):
        # Clear the buffer (for safety).
        for key in data.keys(): data[key][:] = 0

        # Get the truth particles.
        truth_particles = [
            ev.particles[:int(np.minimum(len(ev.particles),n_truth))]
            for ev in truth_events[start_idxs[i]:stop_idxs[i]]
        ]

        ## Get final-state particles.
        final_state_particles = [
            ev.particles[:int(np.minimum(len(ev.particles), 1e4))]
            for ev in final_state_events[start_idxs[i]:stop_idxs[i]]
        ]

        prefix_nzero = int(np.ceil(np.log10(nchunks))) + 1
        prefix_level2 = '\tClustering jets & preparing data for chunk {}/{}:'.format(str(i+1).zfill(prefix_nzero),nchunks)
        if(verbosity==2): printProgressBarColor(0,ranges[i], prefix=prefix_level2, suffix=suffix, length=bl)

        # For each event we must combine tracks and neutral hadrons, perform jet clustering on them,
        # select a single jet (based on user criteria), and select that jet's leading constituents.
        for j in range(ranges[i]):
            data['is_signal'][j] = 1 # TODO: redundant for now

            vecs = np.concatenate([[[final_state_particles[j][k].momentum.px,
                                       final_state_particles[j][k].momentum.py,final_state_particles[j][k].momentum.pz,
                                       final_state_particles[j][k].momentum.e]] for k in range(len(final_state_particles[j]))])

            # print(vecs)
            ## We already have the px, py, pz, e, and we do not need to adjust it, just extract the information directly from the hempc file.

            # ----- Begin jet clustering -----
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

            # Check if there are any jets left. We may have executed a "jet check" before Delphes,
            # but after passing things through Delphes we may no longer have jets that meet our cuts.

            njets = len(jets)
            if(njets == 0):
                #TODO: Determine what's best to do here.
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                continue

            # Now select a jet.
            selected_jet_idx = jet_config['jet_selection'](truth = truth_particles[j], jets = jets, use_hepmc=True)
            if(selected_jet_idx < 0):
                # TODO: Determine what's best to do here.
                data['is_signal'][j] = -1 # Using -1 to mark as "no event". (To be discarded.)
                continue
            jet = jets[selected_jet_idx]

            # Get the constituents of our selected jet. Assuming they are massless -> don't need to fetch the mass.
            pt,eta,phi,m = np.hsplit(np.array([[x.pt(), x.eta(), x.phi(), x.m()] for x in jet.constituents()]),4)
            pt = pt.flatten()
            eta = eta.flatten()
            phi = phi.flatten()
            m = m.flatten()
            # ----- End jet clustering -----

            # Sort by decreasing pt, and only keep leading constituents.
            sorting = np.argsort(-pt)
            l = int(np.minimum(n_constituents,len(pt)))
            vec = PtEtaPhiMToEPxPyPz(pt=pt[sorting][:l], # Vector will be of format (E, px, py, pz)
                                     eta=eta[sorting][:l],
                                     phi=phi[sorting][:l],
                                     m=m[sorting][:1]
                                    )

            data['Nobj'][j] = l
            data['Pmu'][j,:l,:] = vec[:l,:] # indexing for vec is redundant
            # Now take care of truth particles. We probably expect the same number of truth particles
            # per event, but to be safe we will always truncate the list to the requested number at most.
            data['truth_Nobj'][j] = len(truth_particles[j])
            data['truth_Pdg'][j,:] = [par.pid for par in truth_particles[j]]
            data['truth_Pmu'][j,:,:] = [
                [par.momentum.e, par.momentum.px, par.momentum.py, par.momentum.pz]
                for par in truth_particles[j]
            ]
            # Also fill our separate buffers for each truth particle, if requested.
            if(separate_truth_particles):
                for k in range(n_truth):
                    data['truth_Pmu_{}'.format(k)][j,:] = [truth_particles[j][k].momentum.e,truth_particles[j][k].momentum.px,truth_particles[j][k].momentum.py,truth_particles[j][k].momentum.pz]

            # Also store some jet-level information for debugging purposes.
            # We will also store the jet momentum in cylindrical coordinates (pt,eta,phi,m).
            # We make sure that the jet's phi is in [-pi ,pi] for consistency.
            data['jet_Pmu'][j,:]     = [jet.e(), jet.px(), jet.py(), jet.pz()]
            data['jet_Pmu_cyl'][j,:] = [jet.pt(), jet.eta(), AdjustPhi(jet.phi()), jet.m()]

            if(verbosity==2): printProgressBarColor(j+1,ranges[i], prefix=prefix_level2, suffix=suffix, length=bl)

        # We have now filled a chunk, time to write it.
        with h5.File(h5_file, 'a') as f:
            for key in dsets.keys():
                dset = f[key]
                dset[start_idxs[i]:stop_idxs[i]] = data[key][:ranges[i]]
        if(verbosity == 1): printProgressBarColor(i+1,nchunks, prefix=prefix_level1, suffix=suffix, length=bl)
    return h5_file


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

