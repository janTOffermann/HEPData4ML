import sys,os,glob
import numpy as np
import h5py as h5
import argparse as ap

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--input',type=str,required=True)
    parser.add_argument('-o','--output',type=str,required=True)
    parser.add_argument('-n','--nkeys',type=int,required=True,help='Number of truth keys to keep.')
    args = vars(parser.parse_args())
    input = args['input']
    output = args['output']
    nkeys = args['nkeys']

    f = h5.File(input,"r")
    keys = list(f.keys())
    keys.sort()

    output_keys = []
    for key in keys:
        if('truth_Pmu_') not in key:
            output_keys.append(key)
            continue

        idx = int(key.replace('truth_Pmu_',''))
        if(idx < nkeys):
            output_keys.append(key)
            continue

    output_keys.sort()

    g = h5.File(output,"w")
    for key in output_keys:
        g.create_dataset(key,data=f[key][:],compression='gzip',compression_opts=9)
    g.close()
    f.close()
    return

if(__name__ == '__main__'):
    main(sys.argv)

