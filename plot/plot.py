import sys,os
import numpy as np
import ROOT as rt
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import argparse as ap
from plot_util.plot_util import RN

# custom imports
path_prefix = os.path.dirname(os.path.abspath(__file__)) + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util.qol_utils.pdg import pdg_names
from util.calcs import Calculator

def DeltaR(vec1,vec2):
    n = vec1.shape[0]
    assert(vec2.shape[0] == n)
    rvec1 = rt.Math.PtEtaPhiMVector(0.,0.,0.,0.)
    rvec2 = rt.Math.PtEtaPhiMVector(0.,0.,0.,0.)

    dr = np.zeros(n)
    for i in range(n):
        rvec1.SetCoordinates(0.,*vec1[i],0.)
        rvec2.SetCoordinates(0.,*vec2[i],0.)
        dr[i] = rt.Math.VectorUtil.DeltaR(rvec1,rvec2)
    return dr

def SimpleHist(data,binning,title,normalize=True):
    h = rt.TH1D(RN(),title,*binning)
    for entry in data:
        h.Fill(entry)
    if(normalize):
        integral = h.Integral()
        if(integral != 0.):
            h.Scale(1./h.Integral())
    return h

def Root2Plt_hist1d(h,ax,grid=True):
    nbins = h.GetNbinsX()
    x_vals = [h.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins)]
    weights = [h.GetBinContent(i+1) for i in range(nbins)]
    ax.hist(x_vals,bins=len(x_vals),weights=weights)
    ax.set_xlabel(h.GetXaxis().GetTitle())
    ax.set_ylabel(h.GetYaxis().GetTitle())

    if(grid):
        ax.grid()
    PltStatBox(h,ax)
    return

def PltStatBox(h,ax):
    mean = h.GetMean()
    rms = h.GetRMS()
    text = '\n'.join(
        [
        'mean = {:.2e}'.format(mean),
        'RMS  = {:.2e}'.format(rms)
        ]
    )
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(0.05,0.95,text,fontsize=14,verticalalignment='top',bbox=props)
    anchored_text = AnchoredText(text,loc='upper right',prop=dict(family='monospace'))
    ax.add_artist(anchored_text)
    return

def SaveOutput(name,outdir):
    plt.savefig('{}/{}'.format(outdir,name))
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--infile',type=str,help='Input HDF5 file.',required=True)
    parser.add_argument('-o','--outdir',type=str,help='Output directory.',default=None)
    parser.add_argument('-n','--ntruth',type=int,help='Number of truth-level particles to plot for.',default=None)
    args = vars(parser.parse_args())

    infile = args['infile']
    outdir = args['outdir']
    if(outdir is None):
        outdir = os.getcwd()
    n_truth = args['ntruth']
    os.makedirs(outdir,exist_ok=True)

    calculator = Calculator(use_vectorcalcs=False)

    f = h5.File(infile,'r')
    keys = list(f.keys())

    # Pmu = f['Pmu'][:]
    Nobj = f['Nobj'][:]
    # Pmu = GetJagged(f['Pmu'][:],Nobj)
    jet_E = f['jet_Pmu'][:,0]
    jet_Pmu_cyl = f['jet_Pmu_cyl'][:]
    truth_Pmu = f['truth_Pmu'][:]
    truth_Pdg = f['truth_Pdg'][0]

    # Make a plot of the number of jet constituents.
    binning = (200,0.,200.)
    title = ';Number of jet constituents;Fractional Count'
    h = SimpleHist(Nobj,binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    # PltStatBox(h,ax)
    SaveOutput('n.png',outdir)
    ax.set_yscale('log')
    SaveOutput('n_log.png',outdir)
    # plt.savefig('n.png')
    plt.clf()

    # Make plots of the jet kinematics
    binning = (100,0.,1000.)
    title = ';jet $p_{T}$ [GeV];Fractional Count'
    h = SimpleHist(jet_Pmu_cyl[:,0],binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    SaveOutput('jet_pt.png',outdir)
    ax.set_yscale('log')
    SaveOutput('jet_pt_log.png',outdir)
    plt.clf()

    binning = (200,-2.5,2.5)
    title = ';jet $\eta$;Fractional Count'
    h = SimpleHist(jet_Pmu_cyl[:,1],binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    SaveOutput('jet_eta.png',outdir)
    ax.set_yscale('log')
    SaveOutput('jet_eta_log.png',outdir)
    plt.clf()

    binning = (100,-np.pi,np.pi)
    title = ';jet $\phi$;Fractional Count'
    h = SimpleHist(jet_Pmu_cyl[:,2],binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    SaveOutput('jet_phi.png',outdir)
    ax.set_yscale('log')
    SaveOutput('jet_phi_log.png',outdir)
    plt.clf()

    binning = (250,0.,250.)
    title = ';jet $m$ [GeV];Fractional Count'
    h = SimpleHist(jet_Pmu_cyl[:,3],binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    SaveOutput('jet_m.png',outdir)
    ax.set_yscale('log')
    SaveOutput('jet_m_log.png',outdir)
    plt.clf()

    binning = (200,0.,2000.)
    title = ';jet $E$ [GeV];Fractional Count'
    h = SimpleHist(jet_E,binning,title)
    fig,ax = plt.subplots(1,1)
    Root2Plt_hist1d(h,ax)
    SaveOutput('jet_e.png',outdir)
    ax.set_yscale('log')
    SaveOutput('jet_e_log.png',outdir)
    plt.clf()

    plt.close()

    # Make plots for the truth particle kinematics.
    # TODO: Assuming truth_Pdg entries are identical across events.
    #       This is generally a safe assumption, but may break down
    #       in some future use case!
    if(n_truth is None):
        n_truth = len(truth_Pdg)
    truth_names = [pdg_names[x] for x in truth_Pdg]

    quark_counter = 1
    for i,truth_name in enumerate(truth_names):
        if(truth_name in ['u','d','#bar{u}','#bar{d}']):
            truth_names[i] = 'q_{' + str(quark_counter) + '}'
            quark_counter += 1

    for i in range(n_truth):
        # if(i > 0): break
        particle_name = truth_names[i]
        has_subscript = False
        if('_' in particle_name): has_subscript = True

        if(not has_subscript):
            particle_name_fancy = '${}'.format(truth_names[i]) + '_{truth}' + '$'
        else:
            particle_name_fancy = '${}'.format(truth_names[i]).replace('}',', truth}') + '$'
            # print(particle_name_fancy)
        particle_name_fancy = particle_name_fancy.replace('#','\\')
        vec = truth_Pmu[:,i,:]
        print('Making plots for {}.'.format(particle_name))

        vec = truth_Pmu[:,i,:]
        vec_cyl = np.array([calculator.EPxPyPzToPtEtaPhiM_single(*x) for x in vec])
        # vec_cyl = np.array([EPxPyPzToPtEtaPhiM(*x) for x in vec])

        binning = (100,0.,1000.)
        title = ';{}'.format(particle_name_fancy) + ' $p_{T}$ [GeV];Fractional Count'
        h = SimpleHist(vec_cyl[:,0],binning,title)
        fig,ax = plt.subplots(1,1)
        ax.set_yscale('log')
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_pt.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        binning = (200,-2.5,2.5)
        title = ';{}'.format(particle_name_fancy) + ' $\eta$;Fractional Count'
        h = SimpleHist(vec_cyl[:,1],binning,title)
        fig,ax = plt.subplots(1,1)
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_eta.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        binning = (100,-np.pi,np.pi)
        title = ';{}'.format(particle_name_fancy) + ' $\phi$;Fractional Count'
        h = SimpleHist(vec_cyl[:,2],binning,title)
        fig,ax = plt.subplots(1,1)
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_phi.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        binning = (250,0.,250.)
        title = ';{}'.format(particle_name_fancy) + ' $m$ [GeV];Fractional Count'
        h = SimpleHist(vec_cyl[:,3],binning,title)
        fig,ax = plt.subplots(1,1)
        ax.set_yscale('log')
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_m.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        binning = (200,0.,2000.)
        title = ';{}'.format(particle_name_fancy) + ' $E$ [GeV];Fractional Count'
        h = SimpleHist(vec[:,0],binning,title)
        fig,ax = plt.subplots(1,1)
        ax.set_yscale('log')
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_e.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        # Also get the dR between the truth particle and the jet.
        binning = (100,0.,2.)
        title = ';$\Delta R$({}, jet)'.format(particle_name_fancy) + ';Fractional Count'
        dr = DeltaR(vec_cyl[:,1:3],jet_Pmu_cyl[:,1:3])
        h = SimpleHist(dr,binning,title)
        fig,ax = plt.subplots(1,1)
        ax.set_yscale('log')
        Root2Plt_hist1d(h,ax)
        SaveOutput('{}_jet_dr.png'.format(particle_name),outdir)
        plt.clf()
        plt.close()

        # plt.close()

    return

if(__name__ == '__main__'):
    main(sys.argv)