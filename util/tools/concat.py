import sys,os,glob
import numpy as np
import h5py as h5
import argparse as ap
import time

def concatenate(input_patterns,output,copts,verbose=False):
    # Determine what are the input files.
    input_files = []
    for pattern in input_patterns:
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

    f = h5.File(output,'w')

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
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
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
    return

def concatenate_slow(input_patterns,output,copts,verbose=False):
    """
    More file I/O, but fewer files are open at once.
    Should decrease memory usage, and deal with issues
    when too many files are open.
    """
    # Determine what are the input files.
    input_files = []
    for pattern in input_patterns:
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

    infile_0 = h5.File(input_files[0],'r')
    keys = list(infile_0[0].keys())
    metakeys = list(infile_0[0].attrs.keys())
    original_shapes = {key:list(infile_0[key].shape) for key in keys}
    dtypes = {key:infile_0[key].dtype for key in keys}
    infile_0.close()

    nentries = []
    shapes = {}
    for i,input_file in enumerate(input_files):
        infile = h5.File(input_file,'r')
        nentries.append(infile[keys[0]].shape[0])
        infile.close()
    n = np.sum(nentries,dtype=int)

    shapes = {}

    for key in keys:
        shape = original_shapes[key]
        shape[0] = n
        shapes[key] = tuple(shape)

    # output file
    f = h5.File(output,'w')

    for i,key in enumerate(keys):
        if(verbose):
            print('\tKey: {}\t({}/{})'.format(key,i+1,len(keys)))
        data = np.zeros(shapes[key],dtype=dtypes[key].dtype)
        counter = 0
        for j,input_file in enumerate(input_files):
            start = counter
            stop = counter+nentries[j]
            infile = h5.File(input_file,'r')
            data[start:stop] = infile[key][:]
            counter += nentries[j]
            infile.close()
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        del data

    # Now handle the metadata.
    #TODO: This currently assumes all the files have the same metadata keys!
    for key in metakeys:
        values = []
        for input_file in input_files:
            infile = h5.File(input_file,'r')
            values.append(infile.attrs[key])
            infile.close()
        # Determine what all the different metadata values are for this key.
        metadata_values = np.unique(np.concatenate(values))
        data = np.zeros(n,dtype=np.dtype('i4')) # buffer where we will store the final indexes for output.

        counter = 0
        for j,input_file in enumerate(input_files):
            infile = h5.File(input_file,'r')
            start = counter
            stop = counter+nentries[j]
            this_metadata = infile.attrs[key]
            this_metadata_index_mapping = [np.where(metadata_values==x)[0][0] for x in this_metadata]

            original_data = infile[key][:] # indexing w.r.t. this infile's metadata
            infile.close()
            data[start:stop] = [this_metadata_index_mapping[x] for x in original_data]
            counter += nentries[j]
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        f.attrs[key] = metadata_values

    f.close()
    return



def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input',   nargs='+', help='Input file pattern (or a list of patterns). Each pattern can also be provided as a comma-separated string of filenames.', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file.', required=True)
    parser.add_argument('-c', '--compression',type=int,help='Compression level (0-10).',default=9)
    parser.add_argument('-v', '--verbose',type=int,help='Verbosity flag.',default=9)
    parser.add_argument('-s', '--slow',type=int,help='Slower mode, but may be necessary when dealing with large number of file.',default=9)

    args = vars(parser.parse_args())

    input_patterns = args['input']
    output = args['output']
    copts = args['compression']
    verbose = args['verbose'] > 0

    t1 = time.time()
    concatenate(input_patterns,output,copts,verbose)
    t2 = time.time()
    dt = t2 - t1
    print('\tElapsed time: {:.2f} seconds.'.format(dt))

if __name__ == '__main__':
    main(sys.argv)
