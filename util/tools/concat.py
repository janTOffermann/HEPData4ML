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
    nentries = [x[keys[0]].shape[0] for x in infiles]
    shapes = {}

    for key in keys:
        shape = list(infiles[0][key].shape)
        shape[0] = nentries
        shapes[key] = tuple(shape)

    f = h5.File(output,'w')

    for key in keys:
        data = np.concatenate([x[key][:] for x in infiles],axis=0)
        f.create_dataset(key,data=data,compression='gzip',compression_opts=copts)
        del data

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
