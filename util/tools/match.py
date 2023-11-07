# A simple tool for matching two files -- for use if you have a Delphes and non-Delphes version of the same dataset.
# The jet cuts might throw out different events for these two files but they are produced from the same underlying set
# of events, so this tool will match on truth-level data to produce subsets of the two datasets where the events are
# matched between each.
# TODO: This should be replaced with some new functionality in the generation pipeline itself.
import sys,os,glob,pathlib,uuid
import numpy as np
import h5py as h5
import argparse as ap
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../..' )
import util.qol_utils.qol_util as qu

def produce_match(input1,input2,output1,output2,match_key='event_idx',copts=9):
    f1 = h5.File(input1,'r')
    f2 = h5.File(input2,'r')

    match_1 = f1[match_key][:]
    match_2 = f2[match_key][:]
    intersection,idx_1,idx_2 = np.intersect1d(match_1,match_2,assume_unique=True,return_indices=True)
    # NOTE: sort in idx_1, idx_2 shouldn't be needed by construction

    g1 = h5.File(output1,'w')
    keys = sorted(list(f1.keys()))
    metakeys = sorted(list(f1.attrs.keys()))
    for metakey in metakeys:
        g1.attrs[metakey] = f1.attrs[metakey]
    for key in keys:
        g1.create_dataset(key,data=f1[key][idx_1],compression='gzip',compression_opts=copts)
    g1.close()
    f1.close()

    g2 = h5.File(output2,'w')
    keys = sorted(list(f2.keys()))
    metakeys = sorted(list(f2.attrs.keys()))
    for metakey in metakeys:
        g2.attrs[metakey] = f2.attrs[metakey]
    for key in keys:
        g2.create_dataset(key,data=f2[key][idx_2],compression='gzip',compression_opts=copts)
    g2.close()
    f2.close()
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i1','--input1',type=str,required=True)
    parser.add_argument('-i2','--input2',type=str,required=True)
    parser.add_argument('-o1','--output1',type=str,default=None)
    parser.add_argument('-o2','--output2',type=str,default=None)
    parser.add_argument('-k','--key',type=str,default='event_idx',help='Key of (scalar) data field upon which to match.')
    args = vars(parser.parse_args())

    input1 = args['input1']
    input2 = args['input2']
    output1 = args['output1']
    output2 = args['output2']
    match_key = args['key']

    if(output1 is None):
        output1 = '.'.join(input1.split('.')[:-1]) + '_matched.h5'

    if(output2 is None):
        output2 = '.'.join(input2.split('.')[:-1]) + '_matched.h5'

    outdir1 = '/'.join(output1.split('/')[:-1])
    if(outdir1 != ''):
        os.makedirs(outdir1,exist_ok=True)

    outdir2 = '/'.join(output2.split('/')[:-1])
    if(outdir2 != ''):
        os.makedirs(outdir2,exist_ok=True)

    produce_match(input1,input2,output1,output2,match_key,copts=9)

if(__name__ == '__main__'):
    main(sys.argv)