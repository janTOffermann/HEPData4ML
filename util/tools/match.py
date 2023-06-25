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

class Comparator:
    def __init__(self,f1,f2,ntruth=5):
        self.f1 = h5.File(f1,'r')
        self.f2 = h5.File(f2,'r')
        self.ntruth = ntruth

        self.Initialize()

    def Initialize(self):
        self.nevents1 = self.f1['Nobj'].shape[0]
        self.nevents2 = self.f2['Nobj'].shape[0]

        self.truth_Nobj1 = self.f1['truth_Nobj'][:]
        self.truth_Nobj2 = self.f2['truth_Nobj'][:]

        self.truth_Pmu1 = self.f1['truth_Pmu'][:,:self.ntruth,:]
        self.truth_Pmu2 = self.f2['truth_Pmu'][:,:self.ntruth,:]

    def CompareEventsByIndices(self,i1,i2):
        nobj1 = self.truth_Nobj1[i1]
        nobj2 = self.truth_Nobj2[i2]
        if(nobj1 != nobj2): return False
        diff = self.truth_Pmu1[i1] - self.truth_Pmu2[i2]
        s = np.sum(np.abs(diff))
        if(s != 0.):
            # if(s < 100.): print('\t',i1,i2,s)
            return False
        return True

    def Search(self,i1,i2_start=0,i2_end=-1):
        if(i2_start < 0): i2_start = 0
        if(i2_end == -1): i2_end = self.nevents2 - 1
        for i2 in range(i2_start,i2_end+1):
            match = self.CompareEventsByIndices(i1,i2)
            if(match):
                return i2
        return -1

    def __del__(self):
        self.f1.close()
        self.f2.close()

def produce_match(input1,input2,output1,output2,ntruth,buffersize=-1):
    comparator = Comparator(input1,input2,ntruth)
    nevents_1 = comparator.nevents1
    n_success = 0
    n_fail = 0

    i2 = 0
    indices_1 = []
    indices_2 = []

    prefix = 'Matching events:'
    suffix = 'Complete'
    bl = 50
    qu.printProgressBarColor(0,nevents_1,prefix,suffix,1,bl)
    for i in range(0,nevents_1):
        i2_end = -1
        if(buffersize > 0): i2_end = i2 + buffersize
        i2 = comparator.Search(i,i2_start=i2,i2_end=i2_end)
        if(i2 > -1):
            n_success += 1
            indices_1.append(i)
            indices_2.append(i2)

        else: n_fail += 1
        qu.printProgressBarColor(i+1,nevents_1,prefix,suffix,1,bl)

    g1 = h5.File(output1,'w')
    g2 = h5.File(output2,'w')

    attribute_keys = list(comparator.f1.attrs.keys())

    for key in attribute_keys:
        g1.attrs[key] = comparator.f1.attrs[key]
        g2.attrs[key] = comparator.f2.attrs[key]

    keys1 = list(comparator.f1.keys())
    keys2 = list(comparator.f2.keys())

    for key in keys1:
        g1.create_dataset(key,data=comparator.f1[key][:][indices_1],compression='gzip',compression_opts=9)

    for key in keys2:
        g2.create_dataset(key,data=comparator.f2[key][:][indices_2],compression='gzip',compression_opts=9)

    g1.close()
    g2.close()

    print(n_success / (n_success + n_fail) * 100.)
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i1','--input1',type=str,required=True)
    parser.add_argument('-i2','--input2',type=str,required=True)
    parser.add_argument('-o1','--outdir1',type=str,default='match1')
    parser.add_argument('-o2','--outdir2',type=str,default='match2')
    parser.add_argument('-n','--ntruth',type=int,default=5,help='Number of truth 4-momenta to explicitly compare. (first n)')
    parser.add_argument('-buffersize','--buffersize',type=int,default=-1,help='Size of search buffer. If -1, will search whole file. (Advanced usage).')
    args = vars(parser.parse_args())

    ntruth = args['ntruth']
    outdir1 = args['outdir1']
    outdir2 = args['outdir2']
    buffersize = args['buffersize']

    os.makedirs(outdir1,exist_ok=True)
    os.makedirs(outdir2,exist_ok=True)

    input1 = glob.glob(args['input1'])
    input2 = glob.glob(args['input2'])

    input1.sort()
    input2.sort()

    output1 = [outdir1 + '/' + x.split('/')[-1].replace('.h5','_matched.h5') for x in input1]
    output2 = [outdir2 + '/' + x.split('/')[-1].replace('.h5','_matched.h5') for x in input2]

    nfiles = len(input1)
    for i in range(nfiles):
        print('{}/{}'.format(i+1,nfiles))
        print('\t',input1[i])
        print('\t',input2[i])
        produce_match(input1[i],input2[i],output1[i],output2[i],ntruth,buffersize)

if(__name__ == '__main__'):
    main(sys.argv)