import sys, os, glob, uuid
from webbrowser import get
import numpy as np
import pyhepmc_ng as hep
import uproot as ur
import h5py as h5
import subprocess as sub
from util.fastjet import BuildFastjet, ParticleInfo
from util.config import GetNPars, GetJetConfig
from util.calcs import PtEtaPhiMToPxPyPzE, PtEtaPhiMToEPxPyPz, AdjustPhi
from util.qol_util import printProgressBarColor
from util.conv_utils.utils import ClusterJets, ExtractHepMCParticles, FetchJetConstituents, FillDataBuffer, InitFastJet, PrepDataBuffer, ExtractHepMCEvents,PrepDelphesArrays, PrepH5File, PrepIndexRanges, SelectFinalStateParticles

# Convert a set of Delphes ROOT files (modified with truth info!) files to a single HDF5 file.
# This also performs jet clustering of the Delphes output (and associated selections).
def DelphesWithTruthToHDF5(delphes_files, truth_files=None, h5_file=None, nentries_per_chunk=1e4, verbosity=0, separate_truth_particles=True):

    jet_config,jetdef = InitFastJet()

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

    delphes_arr,var_map = PrepDelphesArrays(delphes_files)
    nentries = len(delphes_arr)

    # Also extract truth particle info from the HepMC files.
    truth_events = ExtractHepMCEvents(truth_files)

    nentries_per_chunk = int(nentries_per_chunk)
    data = PrepDataBuffer(nentries_per_chunk,separate_truth_particles)

    # Prepare the HDF5 file.
    dsets = PrepH5File(h5_file,nentries,data)

    # Some indexing preparation for writing in chunks.
    start_idxs,stop_idxs,ranges = PrepIndexRanges(nentries,nentries_per_chunk)
    nchunks = len(ranges)

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
        types = var_map.keys() # TODO: Can this be placed elsewhere? Renamed?
        for delphes_type in types:
            momenta[delphes_type] = {
                x:delphes_arr[var_map[delphes_type][x]][start_idxs[i]:stop_idxs[i]]
                for x in ['pt','eta','phi']
            }

        # Get the truth particles.
        truth_particles = ExtractHepMCParticles(truth_events[start_idxs[i]:stop_idxs[i]],n_truth)

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

            # TODO: Use the below
            # SelectFinalStateParticles(px,py,pz,e, jetdef, truth_particles, data, j, separate_truth_particles)

            jets = ClusterJets(vecs,jetdef,jet_config)

            # Check if there are any jets left. We may have executed a "jet check" before Delphes,
            # but after passing things through Delphes we may no longer have jets that meet our cuts.
            njets = len(jets)
            if(njets == 0):
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
            vecs = FetchJetConstituents(jet,n_constituents,zero_mass=True)
            FillDataBuffer(data,j,vecs,jet,truth_particles,separate_truth_particles)

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

    jet_config,jetdef = InitFastJet()

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

    # Extract truth particle info from the HepMC files.
    truth_events = ExtractHepMCEvents(truth_files,get_nevents=False)

    ## Extract final-state truth particle info from the HepMC files.
    final_state_events, nentries = ExtractHepMCEvents(final_state_files,get_nevents=True)

    # Now create our Numpy buffers to hold data. This dictates the structure of the HDF5 file we're making,
    # each key here corresponds to an HDF5 dataset. The shapes of the datasets will be what is shown here,
    # except with nentries_per_chunk -> nentries.
    nentries_per_chunk = int(nentries_per_chunk)
    data = PrepDataBuffer(nentries_per_chunk,separate_truth_particles)

    # Prepare the HDF5 file.
    dsets = PrepH5File(h5_file,nentries,data)

    # Some indexing preparation for writing in chunks.
    start_idxs,stop_idxs,ranges = PrepIndexRanges(nentries,nentries_per_chunk)
    nchunks = len(ranges)

    # Some setup for (optional printouts).
    prefix_level1 = '\tClustering jets & preparing data:'
    suffix = 'Complete'
    bl = 50

    if(verbosity == 1): printProgressBarColor(0,nchunks, prefix=prefix_level1, suffix=suffix, length=bl)

    for i in range(nchunks):

        # Clear the buffer (for safety).
        for key in data.keys(): data[key][:] = 0

        # Get the truth particles.
        truth_particles = ExtractHepMCParticles(truth_events[start_idxs[i]:stop_idxs[i]],n_truth)

        ## Get final-state particles.
        npar_fs = 1e4
        final_state_particles = ExtractHepMCParticles(final_state_events[start_idxs[i]:stop_idxs[i]],npar_fs)

        prefix_nzero = int(np.ceil(np.log10(nchunks))) + 1
        prefix_level2 = '\tClustering jets & preparing data for chunk {}/{}:'.format(str(i+1).zfill(prefix_nzero),nchunks)
        if(verbosity==2): printProgressBarColor(0,ranges[i], prefix=prefix_level2, suffix=suffix, length=bl)

        # For each event we must select some final-state particles, typically via jet clustering.
        for j in range(ranges[i]):

            l = len(final_state_particles[j])
            px = np.array([final_state_particles[j][k].momentum.px for k in range(l)])
            py = np.array([final_state_particles[j][k].momentum.py for k in range(l)])
            pz = np.array([final_state_particles[j][k].momentum.pz for k in range(l)])
            e  = np.array([final_state_particles[j][k].momentum.e  for k in range(l)])
            SelectFinalStateParticles(px,py,pz,e, jetdef, truth_particles, data, j, separate_truth_particles)

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

