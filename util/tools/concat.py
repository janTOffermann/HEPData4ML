import sys,os,glob
import numpy as np
import h5py as h5
import argparse as ap

def concatenate(input_patterns,output,copts):
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
    for i,file in enumerate(input_files):
        print('\t{}: {}'.format(i,file))
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
    for key in keys:
        data = np.zeros(shapes[key],dtype=infiles[0][key].dtype)
        counter = 0
        for i,infile in enumerate(infiles):
            start = counter
            stop = counter+nentries[i]
            data[start:stop] = infile[key][:]
            counter += nentries[i]
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        del data

    # Now handle the metadata.
    #TODO: This currently assumes all the files have the same metadata keys!
    for key in metakeys:
        # Determine what all the different metadata values are for this key.
        metadata_values = np.unique(np.concatenate([x.attrs[key] for x in infiles]))
        data = np.zeros(n,dtype=np.dtype('i4')) # buffer where we will store the final indexes for output.

        counter = 0
        for i,infile in enumerate(infiles):
            start = counter
            stop = counter+nentries[i]
            this_metadata = infile.attrs[key]
            this_metadata_index_mapping = [np.where(metadata_values==x)[0][0] for x in this_metadata]

            original_data = infile[key][:] # indexing w.r.t. this infile's metadata
            data[start:stop] = [this_metadata_index_mapping[x] for x in original_data]
            counter += nentries[i]
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        f.attrs[key] = metadata_values

    f.close()

    for file in infiles:
        file.close()
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input',   nargs='+', help='Input file pattern (or a list of patterns). Each pattern can also be provided as a comma-separated string of filenames.', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file.', default=None)
    parser.add_argument('-c', '--compression',type=int,help='Compression level (0-10).',default=7)
    args = vars(parser.parse_args())

    input_patterns = args['input']
    output = args['output']
    copts = args['compression']

    concatenate(input_patterns,output,copts)

if __name__ == '__main__':
    main(sys.argv)