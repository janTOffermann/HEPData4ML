import sys,os
import numpy as np
import ROOT as rt
import h5py as h5
import matplotlib as mpl
mpl.use('Agg')
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

class Plotter:
    def __init__(self,outdir=None):
        if(outdir is None): outdir = os.getcwd()
        self.DefaultProps()
        self.SetOutputDirectory(outdir)

    def SetOutputDirectory(self,value):
        self.outdir = value

    def ClearProps(self):
        self.props = {}
        self.DefaultProps()

    def DefaultProps(self):
        self.props = {}
        self.SetAxisFont()
        self.SetPaveFont()
        self.SetPaveFillColor()
        self.SetFaceColor()
        self.SetBackgroundColor()

    def SetStyle(self,style_number=0):

        if(style_number == -1):
            # a little easter egg
            self.SetAxisFont('monospace')
            self.SetPaveFont('monospace')
            self.SetPaveFillColor('white')
            self.SetFaceColor('wheat')
            self.SetBackgroundColor('wheat')
            self.SetAxisFontColor('xkcd:white')

        else:
            self.SetAxisFont('monospace')
            self.SetPaveFont('monospace')
            self.SetPaveFillColor('white')
            self.SetFaceColor('white')
            self.SetBackgroundColor('white')
            self.SetAxisFontColor('black')

    def SetAxisFont(self,font='monospace'):
        self.props['axis_font'] = font
        return

    def SetPaveFont(self,font='monospace'):
        self.props['pave_font'] = font
        return

    def SetPaveFillColor(self,color='wheat'):
        self.props['pave_fill_color'] = color
        return

    def SetFaceColor(self,color='white'):
        self.props['facecolor'] = color
        return

    def SetBackgroundColor(self,color='white'):
        self.props['background_color'] = color
        return

    def SetAxisFontColor(self,color='black'):
        self.props['axis_font_color'] = color
        return

    def SimpleHist(self,data,binning,title,normalize=True):
        h = rt.TH1D(RN(),title,*binning)
        for entry in data:
            h.Fill(entry)
        if(normalize):
            integral = h.Integral()
            if(integral != 0.):
                h.Scale(1./h.Integral())
        return h

    def Root2Plt_hist1d(self,h,grid=True,ymin=None,ymax=None,ylog=False,name='plot'):
        """
        Given a ROOT histogram (TH1), draw it using matplotlib.
        """
        fig,ax = plt.subplots(1,1)
        nbins = h.GetNbinsX()
        x_vals = [h.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins)]
        weights = [h.GetBinContent(i+1) for i in range(nbins)]
        ax.hist(x_vals,bins=len(x_vals),weights=weights)
        ax.set_xlabel(h.GetXaxis().GetTitle(),family=self.props['axis_font'])
        ax.set_ylabel(h.GetYaxis().GetTitle(),family=self.props['axis_font'])

        if(ymin is None): ymin = 1.0e-3
        if(ymax is None): ymax = 1.0
        ax.set_ylim((ymin,ymax))

        if(ylog): ax.set_yscale('log')
        if(grid): ax.grid()
        self.PltStatBox(h,ax)

        ax.set_facecolor(self.props['facecolor'])
        fig.set_facecolor(self.props['background_color'])

        self.SaveOutput(name)
        return

    def PltStatBox(self,h,ax):
        mean = h.GetMean()
        rms = h.GetRMS()
        text = '\n'.join(
            [
            'mean = {:.2e}'.format(mean),
            'RMS  = {:.2e}'.format(rms)
            ]
        )
        props = {
            'family':self.props['pave_font'],
            'backgroundcolor':self.props['pave_fill_color']
        }
        anchored_text = AnchoredText(text,loc='upper right',prop=props)
        ax.add_artist(anchored_text)
        return

    def SaveOutput(self,name,outdir=None):
        if(outdir is None): outdir = self.outdir
        plt.savefig('{}/{}'.format(outdir,name))
        return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--infile',type=str,help='Input HDF5 file.',required=True)
    parser.add_argument('-o','--outdir',type=str,help='Output directory.',default=None)
    parser.add_argument('-n','--ntruth',type=int,help='Number of truth-level particles to plot for.',default=None)
    parser.add_argument('-d','--descriptor',type=str,help='String describing dataset. For multiple lines, use semicolon separation.')
    parser.add_argument('-N','--Nevents',type=int,default=-1,help='Number of events to plot (plots for the first N events). If <= 0, plots from whole dataset.')
    parser.add_argument('-s','--style',type=int,default=0)
    args = vars(parser.parse_args())

    infile = args['infile']
    outdir = args['outdir']
    if(outdir is None):
        outdir = os.getcwd()
    n_truth = args['ntruth']
    descriptor_strings = args['descriptor']
    if(descriptor_strings is not None):
        descriptor_strings = descriptor_strings.split(';')
    nevents = args['Nevents']
    style = args['style']

    os.makedirs(outdir,exist_ok=True)
    calculator = Calculator(use_vectorcalcs=False)

    f = h5.File(infile,'r')
    # keys = list(f.keys())

    print('Loading data into memory. The truth_Pmu arrays will be loaded as needed during plotting.')
    # Pmu = f['Pmu'][:]
    Nobj = f['Nobj'][:]
    if(nevents <= 0):
        nevents = Nobj.shape[0]
    Nobj = f['Nobj'][:nevents]

    # Pmu = GetJagged(f['Pmu'][:],Nobj)
    jet_E = f['jet_Pmu'][:nevents,0]
    jet_Pmu_cyl = f['jet_Pmu_cyl'][:nevents]
    truth_Pdg = f['truth_Pdg'][0]

    # Define a bunch of binning settings.
    binning = {
        'n':(200,0.,200.),
        'pt':(100,0.,1000.),
        'eta':(200,-2.5,2.5),
        'phi':(100,-np.pi,np.pi),
        'm':(250,0.,250.),
        'e':(200,0.,2000.),
        'dr':(100,0.,2.)
    }

    plotter = Plotter(outdir)
    plotter.SetStyle(style)

    print('Making plots for jets.')

    # Make a plot of the number of jet constituents.
    title = ';Number of jet constituents;Fractional Count'
    h = plotter.SimpleHist(Nobj,binning['n'],title)
    plotter.Root2Plt_hist1d(h,name='n.png')
    plotter.Root2Plt_hist1d(h,name='n_log.png',ylog=True)
    # plt.savefig('n.png')
    plt.clf()

    # Make plots of the jet kinematics
    title = ';jet $p_{T}$ [GeV];Fractional Count'
    h = plotter.SimpleHist(jet_Pmu_cyl[:,0],binning['pt'],title)
    plotter.Root2Plt_hist1d(h,name='jet_pt.png')
    plotter.Root2Plt_hist1d(h,name='jet_pt_log.png',ylog=True)
    plt.clf()

    title = ';jet $\eta$;Fractional Count'
    h = plotter.SimpleHist(jet_Pmu_cyl[:,1],binning['eta'],title)
    plotter.Root2Plt_hist1d(h,name='jet_eta.png')
    plotter.Root2Plt_hist1d(h,name='jet_eta_log.png',ylog=True)
    plt.clf()

    title = ';jet $\phi$;Fractional Count'
    h = plotter.SimpleHist(jet_Pmu_cyl[:,2],binning['phi'],title)
    plotter.Root2Plt_hist1d(h,name='jet_phi.png')
    plotter.Root2Plt_hist1d(h,name='jet_phi_log.png',ylog=True)
    plt.clf()

    title = ';jet $m$ [GeV];Fractional Count'
    h = plotter.SimpleHist(jet_Pmu_cyl[:,3],binning['m'],title)
    plotter.Root2Plt_hist1d(h,name='jet_m.png')
    plotter.Root2Plt_hist1d(h,name='jet_m_log.png',ylog=True)
    plt.clf()

    title = ';jet $E$ [GeV];Fractional Count'
    h = plotter.SimpleHist(jet_E,binning['e'],title)
    plotter.Root2Plt_hist1d(h,name='jet_e.png')
    plotter.Root2Plt_hist1d(h,name='jet_e_log.png',ylog=True)
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

        print('Making plots for {}.'.format(particle_name))
        vec = f['truth_Pmu'][:nevents,i,:]

        vec_cyl = np.array([calculator.EPxPyPzToPtEtaPhiM_single(*x) for x in vec])
        # vec_cyl = np.array([EPxPyPzToPtEtaPhiM(*x) for x in vec])

        title = ';{}'.format(particle_name_fancy) + ' $p_{T}$ [GeV];Fractional Count'
        h = plotter.SimpleHist(vec_cyl[:,0],binning['pt'],title)
        plotter.Root2Plt_hist1d(h,name='{}_pt.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        title = ';{}'.format(particle_name_fancy) + ' $\eta$;Fractional Count'
        h = plotter.SimpleHist(vec_cyl[:,1],binning['eta'],title)
        plotter.Root2Plt_hist1d(h,name='{}_eta.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        title = ';{}'.format(particle_name_fancy) + ' $\phi$;Fractional Count'
        h = plotter.SimpleHist(vec_cyl[:,2],binning['phi'],title)
        plotter.Root2Plt_hist1d(h,name='{}_phi.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        title = ';{}'.format(particle_name_fancy) + ' $m$ [GeV];Fractional Count'
        h = plotter.SimpleHist(vec_cyl[:,3],binning['m'],title)
        plotter.Root2Plt_hist1d(h,name='{}_m.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        title = ';{}'.format(particle_name_fancy) + ' $E$ [GeV];Fractional Count'
        h = plotter.SimpleHist(vec[:,0],binning['e'],title)
        plotter.Root2Plt_hist1d(h,name='{}_e.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        # Also get the dR between the truth particle and the jet.
        title = ';$\Delta R$({}, jet)'.format(particle_name_fancy) + ';Fractional Count'
        dr = DeltaR(vec_cyl[:,1:3],jet_Pmu_cyl[:,1:3])
        h = plotter.SimpleHist(dr,binning['dr'],title)
        plotter.Root2Plt_hist1d(h,name='{}_jet_dr.png'.format(particle_name),ylog=True)
        plt.clf()
        plt.close()

        # plt.close()

    return

if(__name__ == '__main__'):
    main(sys.argv)
