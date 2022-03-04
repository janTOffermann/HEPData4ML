import sys,os
import numpy as np
import ROOT as rt
import time
import matplotlib.pyplot as plt

# custom imports
path_prefix = os.getcwd() + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util import qol_util as qu
from util.calcs import PtEtaPhiMToPxPyPzE_root, PtEtaPhiMToPxPyPzE_numpy
from util.vectorcalcs import LoadVectorCalcs

def generate_data(nvecs = 100):
    pt = 100. * np.abs(np.random.rand(nvecs))
    eta = 4. * (np.random.rand(nvecs) - 0.5)
    phi = np.random.rand(nvecs) * (2. * np.pi)
    m = np.zeros(nvecs) # zero masses, like in practice
    data = np.vstack((pt,eta,phi,m)).T
    return np.vstack((pt,eta,phi,m)).T

def run_test(data, ntests=20):
    pt = data[:,0]
    eta = data[:,1]
    phi = data[:,2]
    m = data[:,3]

    # Numpy test
    output = PtEtaPhiMToPxPyPzE_numpy(pt,eta,phi,m) # initial call, not sure if meaningful
    start = time.time()
    for i in range(ntests):
        output = PtEtaPhiMToPxPyPzE_numpy(pt,eta,phi,m)
    end = time.time()
    dt2 = (end - start) / ntests

    # ROOT test
    output = PtEtaPhiMToPxPyPzE_root(pt,eta,phi,m) # initial call, not sure if meaningful
    start = time.time()
    for i in range(ntests):
        output = PtEtaPhiMToPxPyPzE_root(pt,eta,phi,m)
    end = time.time()
    dt1 = (end - start) / ntests

    return dt1,dt2

def sem(x):
    return np.std(x) / np.sqrt(len(x))

def main(args):
    LoadVectorCalcs()
    print('Starting timing test for coordinate conversion.')

    # Prepare some data to test with.
    ntests = 25 # How many tests we do per "batch"
    nbatch = 5 # How many batches -- we average over these, also compute stderr
    dt1_avg = []
    dt1_err = []
    dt2_avg = []
    dt2_err = []

    nvecs = [10,20,50,100,200,500,1000,1500,2000,2500,5000,10000,20000,int(1.0e6)]
    for nvec in nvecs:
        data = generate_data(nvec)

        dt1_vals = []
        dt2_vals = []
        for i in range(nbatch):
            dt1_i,dt2_i = run_test(data, ntests)
            dt1_vals.append(dt1_i)
            dt2_vals.append(dt2_i)

        dt1_avg.append(np.mean(dt1_vals))
        dt2_avg.append(np.mean(dt2_vals))
        dt1_err.append(sem(np.array(dt1_vals)))
        dt2_err.append(sem(np.array(dt2_vals)))

    # Now plot the results -- average execution time vs. number of vectors.
    fig,ax = plt.subplots()
    ax.errorbar(nvecs,dt1_avg,yerr=dt1_err,color='xkcd:red',label='ROOT-based method')
    ax.errorbar(nvecs,dt2_avg,yerr=dt2_err,color='xkcd:blue',label='Pure Numpy method')

    plt.xscale('log')
    plt.yscale('log')
    fontsize = 14
    plt.grid(b=True, which='major', color='0.65', linestyle='--')
    plt.xlabel(r'Number of input vectors',fontsize=fontsize)
    plt.ylabel(r'Execution time [s]',fontsize=fontsize)

    plt.legend()
    plt.savefig('conversion_test.png')

if __name__ == '__main__':
    main(sys.argv)