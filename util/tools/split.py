import sys,os,glob
import numpy as np
import h5py as h5
import argparse as ap
import pathlib

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file.',required=True)
    parser.add_argument('-o','--output', type=str, help='Output directory.',default=None)
    parser.add_argument('-n1', '--ntrain', type=int, help='Number of training events.',required=True)
    parser.add_argument('-n2', '--ntest', type=int, help='Number of testing events.',required=True)
    parser.add_argument('-s', '--shuffle',type=int,help='Whether or not to shuffle events.', default=0)
    parser.add_argument('-c','--compression_opts',type=int,help='Compression level for h5py.',default=7)

    args = vars(parser.parse_args())

    input_file = args['input']
    outdir = args['output']
    ntrain = args['ntrain']
    ntest = args['ntest']
    shuffle = args['shuffle']
    compression_opts = args['compression_opts']

    rng = np.random.default_rng()

    if(outdir is None): outdir = os.getcwd()
    else: os.makedirs(outdir,exist_ok=True)

    f = h5.File(input_file,'r')
    keys = list(f.keys())
    n_total = f[keys[0]].shape[0]

    if(ntrain + ntest > n_total):
        print('ERROR: Requested more events than are present.')
        return

    nvalid = n_total - ntrain - ntest
    indices = np.arange(n_total)

    if(shuffle):
        rng.shuffle(indices)

    indices = {
        'train':indices[:ntrain],
        'test':indices[-ntest:],
        'valid':indices[ntrain:ntrain+nvalid]
    }

    for i,(file_key,n) in enumerate(zip(['train','test','valid'],[ntrain,ntest,nvalid])):
        filename = '{}/{}.h5'.format(outdir,file_key)
        if(pathlib.Path(filename).exists()):
            print('WARNING: File {} exists. Skipping.'.format(filename))

        shapes = {}
        for key in keys:
            shape = list(f[key].shape)
            shape[0] = len(indices[file_key])
            shapes[key] = tuple(shape)
            print(key,shapes[key])

        g = h5.File(filename,'w')

        for key in keys:
            data = f[key][indices[file_key]]
            g.create_dataset(key,data=data,compression='gzip',compression_opts=compression_opts)
            del data
        g.close()
    f.close()

if __name__ == '__main__':
    main(sys.argv)