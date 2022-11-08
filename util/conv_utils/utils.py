import sys, glob
import h5py as h5
import numpy as np
import uproot as ur
import pyhepmc as hep
from util.fastjet import BuildFastjet
from util.config import GetNPars, GetJetConfig

# --- FASTJET IMPORT ---
# TODO: Can this be done more nicely?
fastjet_dir = BuildFastjet(j=8)
fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
import fastjet as fj
# ----------------------

def InitFastJet():
    fj.ClusterSequence.print_banner() # Get the Fastjet banner out of the way. TODO: Can we get rid of the banner? It is annoying.
    jet_config = GetJetConfig()
    jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])
    return jet_config,jetdef

def ExtractHepMCEvents(files,get_nevents=False):
    events = []
    nevents = 0
    for file in files:
        with hep.io.ReaderAscii(file) as f:
            for evt in f:
                events.append(evt)
                if(get_nevents): nevents += 1
    if(get_nevents): return events, nevents
    return events

def ExtractHepMCParticles(events,n_par):
    particles = [
        ev.particles[:int(np.minimum(len(ev.particles),n_par))]
        for ev in events
    ]
    return particles

def PrepDelphesArrays(delphes_files):
    # types = ['EFlowTrack','EFlowNeutralHadron','EFlowPhoton']
    types = ['Tower'] # TODO: Is there a better way to prepare Delphes jet constituents? We want PFlow/EFlow?
    components = ['PT','Eta','Phi','ET'] # not all types have all components, this is OK
    delphes_keys = ['{x}.{y}'.format(x=x,y=y) for x in types for y in components]
    delphes_tree = 'Delphes'
    delphes_trees = ['{}:{}'.format(x,delphes_tree) for x in delphes_files]
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

def PrepDataBuffer(nentries_per_chunk,separate_truth_particles=False, n_separate=-1):
    nentries_per_chunk = int(nentries_per_chunk)
    npars = GetNPars()
    n_constituents = npars['jet_n_par']
    n_truth = npars['n_truth']

    # Create our Numpy buffers to hold data. This dictates the structure of the HDF5 file we're making,
    # each key here corresponds to an HDF5 dataset. The shapes of the datasets will be what is shown here,
    # except with nentries_per_chunk -> nentries.
    data = {
        'Nobj'        : np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of jet constituents
        'Pmu'         : np.zeros((nentries_per_chunk,n_constituents,4),dtype=np.dtype('f8')), # 4-momenta of jet constituents (E,px,py,pz)
        'truth_Nobj'  : np.zeros(nentries_per_chunk,dtype=np.dtype('i2')), # number of truth-level particles (this is somewhat redundant -- it will typically be constant)
        'truth_Pdg'   : np.zeros((nentries_per_chunk,n_truth),dtype=np.dtype('i4')), # PDG codes to ID truth particles
        'truth_Pmu'   : np.zeros((nentries_per_chunk,n_truth,4),dtype=np.dtype('f8')), # truth-level particle 4-momenta
        'is_signal'   : np.zeros(nentries_per_chunk,dtype=np.dtype('i1')), # signal flag (0 = background, 1 = signal)
        'jet_Pmu'     : np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')), # jet 4-momentum, in Cartesian coordinates (E, px, py, pz)
        'jet_Pmu_cyl' : np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8')) # jet 4-momentum, in cylindrical coordinates (pt,eta,phi,m)
    }

    # In addition to the above, we will also store the truth-level 4-momenta with each in its own HDF5 dataset.
    # In practice, this might be a more convenient storage format for things like neural network training.
    # Note, however, that the number of these keys will depend on n_truth (the number of truth particles saved per event).
    # TODO: Would be nice to accomplish this with references if possible, see: https://docs.h5py.org/en/stable/refs.html#refs .
    #       Not clear if that is actually feasible.
    if(separate_truth_particles):
        if(n_separate <= 0): n_separate = n_truth
        for i in range(n_separate):
            key = 'truth_Pmu_{}'.format(i)
            data[key] = np.zeros((nentries_per_chunk,4),dtype=np.dtype('f8'))
    return data

def PrepH5File(filename,nentries,data_buffer):
    dsets = {}
    with h5.File(filename, 'w') as f:
        for key, val in data_buffer.items():
            shape = list(val.shape)
            shape[0] = nentries
            shape = tuple(shape)
            dsets[key] = f.create_dataset(key, shape, val.dtype,compression='gzip')
    return dsets

def PrepIndexRanges(nentries,nentries_per_chunk):
    nchunks = int(np.ceil(nentries / nentries_per_chunk))
    start_idxs = np.zeros(nchunks,dtype = np.dtype('i8'))
    for i in range(1,start_idxs.shape[0]): start_idxs[i] = start_idxs[i-1] + nentries_per_chunk
    stop_idxs = start_idxs + nentries_per_chunk
    stop_idxs[-1] = nentries
    ranges = stop_idxs - start_idxs
    return start_idxs,stop_idxs,ranges