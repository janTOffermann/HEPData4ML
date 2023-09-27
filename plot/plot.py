import sys,os,datetime,glob
import numpy as np
import ROOT as rt
import h5py as h5
import matplotlib.pyplot as plt
import argparse as ap
from plot_util.plot_util import RN, MetaDataHandler,Plotter

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

class DataLoader:
    def __init__(self,infiles):
        self.infiles = infiles
        self.h5_files = []
        self.Initialize()

    def __del__(self):
        for hfile in self.h5_files:
            hfile.close()

    def Initialize(self):
        self.h5_files = [h5.File(x,'r') for x in self.infiles]

    def Get(self,key):
        return np.concatenate([x[key][:] for x in self.h5_files],axis=0)

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--infiles',type=str,help='Input HDF5 file(s). Can be a single filename, or a glob-compatible string',required=True)
    parser.add_argument('-o','--outdir',type=str,help='Output directory.',default=None)
    parser.add_argument('-n','--ntruth',type=int,help='Number of truth-level particles to plot for.',default=None)
    parser.add_argument('-d','--descriptor',type=str,help='String describing dataset. For multiple lines, use semicolon separation.')
    parser.add_argument('-N','--Nevents',type=int,default=-1,help='Number of events to plot (plots for the first N events). If <= 0, plots from whole dataset.')
    parser.add_argument('-s','--style',type=int,default=0)
    parser.add_argument('-smoothing','--smoothing',type=float,default=0,help='Smoothing coefficient for 2D matplotlib plots. If <= 0, no smoothing is applied.')
    args = vars(parser.parse_args())

    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptStat(0)

    infiles = glob.glob(args['infiles'])
    infiles.sort()

    outdir = args['outdir']
    if(outdir is None):
        outdir = os.getcwd()
    n_truth = args['ntruth']
    descriptor_strings = args['descriptor']
    if(descriptor_strings is not None):
        descriptor_strings = descriptor_strings.split(';')
    nevents = args['Nevents']
    style = args['style']
    smoothing = args['smoothing']

    os.makedirs(outdir,exist_ok=True)
    calculator = Calculator(use_vectorcalcs=False)
    dataloader = DataLoader(infiles)

    # print('Loading data into memory. The truth_Pmu arrays will be loaded as needed during plotting.')
    print('Loading data into memory.')

    Nobj = dataloader.Get('Nobj')
    if(nevents <= 0):
        nevents = Nobj.shape[0]
    Nobj = Nobj[:nevents]
    jet_E = dataloader.Get('jet_Pmu')[:nevents,0]
    jet_Pmu_cyl = dataloader.Get('jet_Pmu_cyl')[:nevents]
    truth_Pdg = dataloader.Get('truth_Pdg')[0]
    truth_Pmu = dataloader.Get('truth_Pmu')

    # Define a bunch of binning settings.
    binning = {
        'n':(100,0.,100.),
        'pt':(100,0.,1000.),
        'eta':(200,-2.5,2.5),
        'phi':(100,-np.pi,np.pi),
        'm':(250,0.,250.),
        'e':(200,0.,2000.),
        'jet_dr':(100,0.,2.)
    }

    metahandler = MetaDataHandler(infiles)
    plotter = Plotter(outdir)
    plotter.SetStyle(style)
    plotter.SetNevents(nevents) # for saving as metadata, useful for turning plots of fractional count into plots of count

    # Now we set a bunch of text boxes. A bit hacky since ROOT uses 2 corner positions,
    # whereas matplotlib just takes 1 set (and scales the size based on font size).

    # Some metadata.
    git_hashes = metahandler.GetMetaData()['git_hash']
    git_hashes.sort()
    nhash = len(git_hashes)

    if(nhash == 1):
        git_hash_text = [
            'Git hash: {}'.format(git_hashes[0])
        ]
    else:
        git_hash_text = [
            'Git hashes:',
        ]

        hash_text = ''
        for i in range(np.minimum(nhash,3)):
            if(i > 0): hash_text += ', '
            hash_text += '{}'.format(git_hashes[i])

        git_hash_text.append(hash_text)
        if(nhash > 3):
            git_hash_text.append('(+ {} more)'.format(nhash - 3))
    unique_ids = metahandler.GetMetaData()['unique_id_short']
    unique_ids.sort()
    unique_id_text = [
        'Dataset ID (short):',
        '{}, {}, {}'.format(*unique_ids[:3]),
        '(+ {} more)'.format(len(unique_ids)-3)
        ]

    timestamps = metahandler.GetMetaData()['timestamp']
    min_timestamp = np.min(timestamps)
    max_timestamp = np.max(timestamps)
    time_format = '%Y-%m-%d %H:%M:%S'
    min_timestamp = datetime.datetime.utcfromtimestamp(min_timestamp).strftime(time_format)
    max_timestamp = datetime.datetime.utcfromtimestamp(max_timestamp).strftime(time_format)
    time_stamp_text = [
        'Data Timestamps (UTC):',
        '{} (earliest)'.format(min_timestamp),
        '{} (latest)'.format(max_timestamp)
    ]

    plotter.AddTextBox(0.0,0.7,0.4,0.975,descriptor_strings,loc='lower left',alpha=0.7) # custom info in upper left.
    plotter.AddTextBox(0.7,0.7,1.0,0.975,git_hash_text + unique_id_text,loc='lower right',alpha=0.7) # custom info in upper left.
    plotter.AddTextBox(0.4,0.7,0.6,0.975,time_stamp_text,loc='lower left',alpha=0.7) # custom info in upper left.

    print('Making plots for jets.')

    jet_titles = {
        'n':'Number of jet constituents',
        'pt':'jet $p_{T}$ [GeV]',
        'eta':'jet $\eta$',
        'phi':'jet $\phi$',
        'm':'jet $m$ [GeV]',
        'e':'jet $E$ [GeV]'
    }

    jet_data = {
        'n'  : Nobj,
        'pt' : jet_Pmu_cyl[:,0],
        'eta': jet_Pmu_cyl[:,1],
        'phi': jet_Pmu_cyl[:,2],
        'm'  : jet_Pmu_cyl[:,3],
        'e'  : jet_E
    }

    # Some 1D jet kinematic plots.
    normalize = True
    for key in jet_titles.keys():
        if(normalize):
            title = ';{};Fractional Count'.format(jet_titles[key])
        else:
            title = ';{};Count'.format(jet_titles[key])
        h = plotter.SimpleHist(jet_data[key],binning[key],title,normalize=normalize)
        name = 'jet_{}'.format(key)
        plotter.Root2Plt_hist1d(h,name=name,ylog=True)
        plotter.Plt2RootRelabeling(h)
        plotter.SaveRootOutputCanvas(h,None,name,ndim=1,draw_option='HIST')
        plotter.SaveRootOutput(h,name)

    # Some 2D jet kinematic plots.
    drawn_keys = []
    for key_x in jet_titles.keys():
        drawn_keys.append(key_x)
        for key_y in jet_titles.keys():
            if(key_x == key_y or key_y in drawn_keys): continue
            if(normalize):
                title = ';{};{};Fractional Count'.format(jet_titles[key_x],jet_titles[key_y])
            else:
                title = ';{};{};Count'.format(jet_titles[key],jet_titles[key_y])
            h2 = plotter.SimpleHist2D(jet_data[key_x],jet_data[key_y],binning[key_x],binning[key_y],title,normalize=normalize)
            name = 'jet_{}_vs_{}'.format(key_y,key_x)
            plotter.Plt_hist2d(jet_data[key_x],jet_data[key_y],binning[key_x],binning[key_y],title,normalize=normalize,name=name,smoothing=smoothing)
            # plotter.Root2Plt_hist2d(h2,name=name,zlog=True)
            plotter.Plt2RootRelabeling(h2)
            plotter.SaveRootOutputCanvas(h2,None,name,ndim=2,draw_option='COLZ')
            plotter.SaveRootOutput(h2,name,ndim=2)

    plt.close()

    # Make plots for the truth particle kinematics.
    # TODO: Assuming truth_Pdg entries are identical across events.
    #       This is generally a safe assumption, but may break down
    #       in some future use case!
    if(n_truth is None):
        n_truth = len(truth_Pdg)
    truth_names = [pdg_names[x] for x in truth_Pdg]

    quark_counter = 0
    for i,truth_name in enumerate(truth_names):
        if(truth_name in ['u','d','c','s','#bar{u}','#bar{d}','#bar{c}','#bar{s}']):
            truth_names[i] = 'q_{' + str(quark_counter+1) + '}'
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
        vec = truth_Pmu[:nevents,i,:]

        vec_cyl = np.array([calculator.EPxPyPzToPtEtaPhiM_single(*x) for x in vec])

        titles = {
            'pt' :'{} '.format(particle_name_fancy) + ' $p_{T}$ [GeV]',
            'eta':'{} '.format(particle_name_fancy) +' $\eta$',
            'phi':'{} '.format(particle_name_fancy) +' $\phi$',
            'm'  :'{} '.format(particle_name_fancy) +' $m$ [GeV]',
            'e'  :'{} '.format(particle_name_fancy) +' $E$ [GeV]',
            'jet_dr': '$\Delta R$(' + '{}'.format(particle_name_fancy) + ', jet)'
        }

        data = {
            'pt' : vec_cyl[:,0],
            'eta': vec_cyl[:,1],
            'phi': vec_cyl[:,2],
            'm'  : vec_cyl[:,3],
            'e'  : vec[:,0],
            'jet_dr' : DeltaR(vec_cyl[:,1:3],jet_Pmu_cyl[:,1:3]) # TODO: Replace with Calculator usage
        }

        # Some 1D truth particle kinematic plots.
        normalize = True
        for key in titles.keys():
            if(normalize):
                title = ';{};Fractional Count'.format(titles[key])
            else:
                title = ';{};Count'.format(titles[key])
            h = plotter.SimpleHist(data[key],binning[key],title,normalize=normalize)
            name = '{}_{}'.format(particle_name,key)
            plotter.Root2Plt_hist1d(h,name=name,ylog=True)
            plotter.Plt2RootRelabeling(h)
            plotter.SaveRootOutputCanvas(h,None,name,ndim=1,draw_option='HIST')
            plotter.SaveRootOutput(h,name)
            plt.close()

        # Some 2D truth particle kinematic plots.
        normalize = True
        drawn_keys = []
        for key_x in titles.keys():
            drawn_keys.append(key_x)
            for key_y in titles.keys():
                if(key_x == key_y or key_y in drawn_keys): continue
                if(normalize):
                    title = ';{};{};Fractional Count'.format(titles[key_x],titles[key_y])
                else:
                    title = ';{};{};Count'.format(titles[key],titles[key_y])
                h2 = plotter.SimpleHist2D(data[key_x],data[key_y],binning[key_x],binning[key_y],title,normalize=normalize)
                name = '{}_{}_vs_{}'.format(particle_name,key_y,key_x)
                plotter.Plt_hist2d(data[key_x],data[key_y],binning[key_x],binning[key_y],title,normalize=normalize,name=name,smoothing=smoothing)
                plotter.Plt2RootRelabeling(h2)
                plotter.SaveRootOutputCanvas(h2,None,name,ndim=2,draw_option='COLZ')
                plotter.SaveRootOutput(h2,name,ndim=2)

        # Some 2D plots -- truth particle kinematics vs. jet kinematics.
        for key_x in jet_titles.keys():
            for key_y in titles.keys():
                if(normalize):
                    title = ';{};{};Fractional Count'.format(jet_titles[key_x],titles[key_y])
                else:
                    title = ';{};{};Count'.format(titles[key],titles[key_y])
                h2 = plotter.SimpleHist2D(jet_data[key_x],data[key_y],binning[key_x],binning[key_y],title,normalize=normalize)
                name = '{}_{}_vs_jet_{}'.format(particle_name,key_y,key_x)
                plotter.Plt_hist2d(jet_data[key_x],data[key_y],binning[key_x],binning[key_y],title,normalize=normalize,name=name,smoothing=smoothing)
                plotter.Plt2RootRelabeling(h2)
                plotter.SaveRootOutputCanvas(h2,None,name,ndim=2,draw_option='COLZ')
                plotter.SaveRootOutput(h2,name,ndim=2)

    # Save nevents to the ROOT files.
    plotter.SaveNevents()
    return

if(__name__ == '__main__'):
    main(sys.argv)
