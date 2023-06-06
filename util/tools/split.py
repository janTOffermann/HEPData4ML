import sys,os
import argparse as ap
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../..')
from util.conversion import SplitH5

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file.',required=True)
    parser.add_argument('-o','--output', type=str, help='Output directory.',default=None)
    parser.add_argument('-f1', '--ftrain', type=float, help='Fraction of input events for training.',required=True)
    parser.add_argument('-f2', '--ftest', type=float, help='Fraction of input events for testing',required=True)
    parser.add_argument('-s', '--seed',type=int,help='RNG seed for shuffle before split', default=1)
    parser.add_argument('-c','--compression_opts',type=int,help='Compression level for h5py.',default=9)

    args = vars(parser.parse_args())

    input_file = args['input']
    outdir = args['output']
    ftrain = args['ftrain']
    ftest = args['ftest']
    seed = args['seed']
    compression_opts = args['compression_opts']

    if(outdir is None): outdir = os.getcwd()
    os.makedirs(outdir,exist_ok=True)

    if(ftrain + ftest > 1):
        print('ERROR: Requested more events than are present.')
        return
    fvalid = 1. - ftrain - ftest

    train_name = '{}/train.h5'.format(outdir)
    test_name = '{}/test.h5'.format(outdir)
    valid_name = '{}/valid.h5'.format(outdir)

    SplitH5(input_file,(ftrain,ftest,fvalid),train_name,valid_name,test_name,copts=compression_opts,verbose=True,seed=seed)

if __name__ == '__main__':
    main(sys.argv)