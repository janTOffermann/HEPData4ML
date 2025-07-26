import os,glob,uuid
import subprocess as sub
import numpy as np
import h5py as h5


def ConcatenateH5(input_file_patterns,output_file,cwd=None,delete_inputs=False,compression='gzip',copts=9,ignore_keys=None,verbose=False, silent_drop=False):
    # Determine what are the input files.
    input_files = []
    for pattern in input_file_patterns:

        if(cwd is not None):
            pattern = '{}/{}'.format(cwd,pattern)

        print('Checking pattern = ',pattern)

        if(',' in pattern):
            input_files += [x.strip() for x in pattern.split(',')]
        elif('*' in pattern or '[' in pattern):
            input_files += glob.glob(os.path.expanduser(pattern),recursive=True)
        else:
            input_files += [pattern]

    input_files.sort()
    print('Concatenating files:')
    for j,file in enumerate(input_files):
        print('\t{}: {}'.format(j,file))
    infiles = [h5.File(x,'r') for x in input_files]

    # Determine the keys -- we assume they are the same for all files!
    keys = list(infiles[0].keys())

    if(ignore_keys is not None):
        keys = [k for k in keys if k not in ignore_keys]
        if(verbose and not silent_drop):
            print('\tExcluding the following keys from concatenation, these will be dropped:')
            for k in ignore_keys:
                print('\t\t{}'.format(k))

    # Exclude keys that also appear in metadata -- for our use these have a special meaning,
    # and have to be handled separately. These keys correspond to scalar integer fields, which
    # are meant to reference an index in the associated metadata (that shares the same key).
    metakeys = list(infiles[0].attrs.keys())
    keys = [x for x in keys if x not in metakeys]

    nentries = [x[keys[0]].shape[0] for x in infiles]
    n = np.sum(nentries,dtype=int)
    shapes = {}

    for key in keys:
        shape = list(infiles[0][key].shape)
        shape[0] = n
        shapes[key] = tuple(shape)

    f = h5.File(output_file,'w')

    # Not sure how this method compares to just using numpy concatenate?
    for i,key in enumerate(keys):
        if(verbose):
            print('\tKey: {}\t({}/{})'.format(key,i+1,len(keys)))
        data = np.zeros(shapes[key],dtype=infiles[0][key].dtype)
        counter = 0
        for j,infile in enumerate(infiles):
            start = counter
            stop = counter+nentries[j]
            data[start:stop] = infile[key][:]
            counter += nentries[j]
        f.create_dataset(key,data=data,compression=compression,compression_opts=copts)
        del data

    # Now handle the metadata.
    #TODO: This currently assumes all the files have the same metadata keys!
    for key in metakeys:
        # Determine what all the different metadata values are for this key.
        metadata_values = np.unique(np.concatenate([x.attrs[key] for x in infiles]))
        data = np.zeros(n,dtype=np.dtype('i4')) # buffer where we will store the final indexes for output.

        counter = 0
        for j,infile in enumerate(infiles):
            start = counter
            stop = counter+nentries[j]
            this_metadata = infile.attrs[key]
            this_metadata_index_mapping = [np.where(metadata_values==x)[0][0] for x in this_metadata]

            original_data = infile[key][:] # indexing w.r.t. this infile's metadata
            data[start:stop] = [this_metadata_index_mapping[x] for x in original_data]
            counter += nentries[j]
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        f.attrs[key] = metadata_values

    f.close()

    for file in infiles:
        file.close()

    if(verbose):
        if(delete_inputs):
            print('\tDeleting input files.')
        print('\t' + 21*"#" + '\n')

    if(delete_inputs):
        for f in input_files:
            sub.check_call(['rm',f])

    return


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
    as those with a negative value in the "SignalFlag" dataset.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)

    fname_tmp = h5_file.replace('.h5','_{}.h5'.format(str(uuid.uuid4())))

    f = h5.File(h5_file,'r')
    keys = list(f.keys())
    shapes = [f[key].shape for key in keys]
    N = shapes[0][0] # number of events

    keep_indices = f['SignalFlag'][:] >= 0
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
def AddEventIndices(h5_file,cwd=None,copts=9,key='Event.Index',offset=0):
    """
    Add a dataset with event indices to an HDF5 file. This is a sequential index,
    which may be useful if some entries will later be dropped and one wants to
    keep track of which events these were.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    nevents = f['SignalFlag'].shape[0]
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
    nevents = f['SignalFlag'].shape[0]

    if(dtype is None): dtype = np.dtype('f8')
    if(dtype != str):
        data = np.full((nevents),value)
    else:
        data = nevents * [value]
    f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
    f.close()

# def AddMetaData(h5_file,cwd=None,value='',key='metadata'):
#     """
#     Generic function for adding a value to the HDF5 file
#     metadata container.
#     """
#     if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
#     f = h5.File(h5_file,'r+')
#     f.attrs[key] = value
#     f.close()

# def AddMetaDataWithReference(h5_file,cwd=None,value='',key='metadata',copts=9):
#     """
#     Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
#     Also creates a column in the dataset that will point to this metadata's index.
#     Somewhat redundant for file generation but this type of logic will be useful when concatenating files
#     with different entries in the metadata fields.
#     """
#     if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
#     f = h5.File(h5_file,'r+')
#     nevents = f['SignalFlag'].shape[0]
#     metadata = f.attrs
#     if(key not in metadata.keys()): f.attrs[key] = [value]
#     else: f.attrs[key] = list(f.attrs[key]) + [value] # I think the list <-> array stuff should be OK here
#     idx = len(f.attrs[key]) - 1
#     f.create_dataset(key,data=np.full(nevents,idx,dtype=np.dtype('i4')),compression='gzip',compression_opts=copts)
#     f.close()

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