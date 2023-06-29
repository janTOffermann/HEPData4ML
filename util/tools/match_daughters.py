# A tool for computing "energy_ratio_truth" for truth-level datasets
# (where the ratio will always be 0 or 1). For advanced usage.
#

import sys,os,glob,pathlib
import numpy as np
import h5py as h5
import subprocess as sub
import argparse as ap
sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../..' )
import util.qol_utils.qol_util as qu

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputs',type=str,required=True,help='Input HDF5 file(s), can be glob-compatible string.')
    parser.add_argument('-n','--n_ignore',type=int,default=5,help='Ignore the first n truth_Pmu for matching purposes.')
    args = vars(parser.parse_args())

    infiles = glob.glob(args['inputs'])

    for infile in infiles:
        print('Input file = {}'.format(infile))
        n_ignore = args['n_ignore']
        outfile = infile.replace('.h5','_matching_indices.h5')
        assert(not pathlib.Path(outfile).exists())

        sub.check_call(['cp',infile,outfile])
        f = h5.File(outfile,'a')
        Pmu = f['Pmu'][:]
        Nobj = f['Nobj'][:]
        truth_Nobj = f['truth_Nobj'][:] - n_ignore
        truth_Pmu = f['truth_Pmu'][:,n_ignore:,:]
        nevents = Nobj.shape[0]

        flag = np.full(Pmu.shape[:2],0.,dtype=np.dtype('f8'))

        prefix = 'Matching Pmu and truth_Pmu'
        suffix = 'Complete'
        bl = 50
        qu.printProgressBarColor(0,nevents,prefix=prefix,suffix=suffix,decimals=2,length=bl)

        for i in range(nevents):
            nobj = Nobj[i]
            for j in range(nobj):
                diff = [Pmu[i,j] - x for x in truth_Pmu[i,:truth_Nobj[i]]]
                diff = np.array([np.sum(np.abs(x)) for x in diff])
                idxs = np.where(diff==0.)[0]
                if(len(idxs == 0)): continue
                flag[i,j] = 1.
                # flag[i,j] = idxs[0]
            qu.printProgressBarColor(i+1,nevents,prefix=prefix,suffix=suffix,decimals=2,length=bl)

        f.create_dataset('energy_ratio_truth',data=flag,compression='gzip',compression_opts=9)
        f.close()
        print('\tOutput file = {}\n'.format(outfile))

if(__name__=='__main__'):
    main(sys.argv)