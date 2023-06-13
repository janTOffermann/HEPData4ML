import sys,os,datetime,glob
import numpy as np
import ROOT as rt
import h5py as h5
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import argparse as ap
from scipy.ndimage import gaussian_filter
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

class MetaDataHandler:
    def __init__(self,filenames):
        files = [h5.File(x,'r') for x in filenames]
        attrs_list = [dict(x.attrs).copy() for x in files]

        # Now merge the dictionary is attrs_list
        self.attrs = {}
        keys = [list(x.keys()) for x in attrs_list]
        keys = list(set([j for i in keys for j in i]))

        for key in keys:
            self.attrs[key] = []
            for attrs in attrs_list:
                self.attrs[key].append(attrs[key])

            self.attrs[key] = np.unique(np.concatenate(self.attrs[key],axis=0))
            # print(key,self.attrs[key])

        for f in files:
            f.close()

    def GetMetaData(self):
        return self.attrs

class Plotter:
    def __init__(self,outdir=None):
        if(outdir is None): outdir = os.getcwd()
        self.DefaultProps()
        self.SetOutputDirectory(outdir)
        self.SetPltExtensions(['pdf'])
        self.SetRootOutputFiles('plots_1d.root','plots_2d.root','canvas_1d.root','canvas_2d.root')
        self.stat_box = False

        self.canvases = [] # for persistence
        self.paves = []

        for filename in [self.root_filename_1d,self.root_filename_2d,self.root_filename_canvas_1d,self.root_filename_canvas_2d]:
            try:
                os.unlink('{}/{}'.format(outdir,filename))
            except:
                pass

        self.text_boxes = []

    def AddTextBox(self,x1,y1,x2,y2,text_list,loc='upper left',alpha=0.2):
        self.text_boxes.append((x1,y1,x2,y2,text_list,loc,alpha)) # will be converted into textbox on-the-fly

    def ClearTextBoxes(self):
        self.text_boxes.clear()

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

    def SetPltExtensions(self,vals):
        self.plt_extensions = vals

    def SetRootOutputFile1D(self,filename):
        self.root_filename_1d = filename

    def SetRootOutputFile2D(self,filename):
        self.root_filename_2d = filename

    def SetRootOutputFiles(self,filename_1d, filename_2d,filename_canvas_1d,filename_canvas_2d):
        self.root_filename_1d = filename_1d
        self.root_filename_2d = filename_2d
        self.root_filename_canvas_1d = filename_canvas_1d
        self.root_filename_canvas_2d = filename_canvas_2d

    def SetStyle(self,style_number=0):

        if(style_number == -1):
            # a little easter egg
            self.SetAxisFont('monospace')
            self.SetPaveFont('monospace')
            self.SetPaveFillColor('white')
            self.SetFaceColor('wheat')
            self.SetBackgroundColor('xkcd:light grey')
            self.SetAxisFontColor('xkcd:white')

        elif(style_number == 1):
            # TODO: Some cool style. Dark mode?
            self.SetAxisFont('monospace')
            self.SetPaveFont('monospace')
            self.SetPaveFillColor('white')
            self.SetFaceColor('white')
            self.SetBackgroundColor('white')
            self.SetAxisFontColor('black')

        else:
            # Default style -- simple.
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

    def SimpleHist2D(self,data_x,data_y,binning_x,binning_y,title,normalize=True):
        h = rt.TH2D(RN(),title,*binning_x,*binning_y)
        assert(len(data_x) == len(data_y))
        for (x,y) in zip(data_x,data_y):
            h.Fill(x,y)
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
        if(self.stat_box): self.PltStatBox(h,ax)

        ax.set_facecolor(self.props['facecolor'])
        fig.set_facecolor(self.props['background_color'])

        self.PltTextBox(ax)

        self.SavePltOutput(name)
        return

    def Plt_hist2d(self,data_x,data_y,binning_x,binning_y,title,normalize=True,name='plot',smoothing=0):
        """
        Make a matplotlib 2D histogram. Would be nicer to convert from existing ROOT hist,
        but this seems ludicrously complicated or at least totally unclear from matplotlib's documentation.
        """

        fig,ax = plt.subplots(1,1)

        norm = mpl.colors.SymLogNorm(linthresh=1.0e-3,vmin=1.0e-4,vmax=1.0e-1)
        if(not normalize): norm = mpl.colors.SymLogNorm(linthresh=1.0e-1)

        if(smoothing > 0):
            heatmap, xedges, yedges = np.histogram2d(
                data_x,
                data_y,
                bins=(
                np.linspace(binning_x[1],binning_x[2],binning_x[0],endpoint=True),
                np.linspace(binning_y[1],binning_y[2],binning_y[0],endpoint=True),
                ),
                density=normalize
            )
            heatmap = gaussian_filter(heatmap,sigma=smoothing)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            c = ax.imshow(heatmap,extent=extent,origin='lower', cmap=plt.cm.jet, norm=norm)
            ax.set_aspect('auto')

        else:
            _,_,_,c = ax.hist2d(
                data_x,
                data_y,
                (
                    np.linspace(binning_x[1],binning_x[2],binning_x[0],endpoint=True),
                    np.linspace(binning_y[1],binning_y[2],binning_y[0],endpoint=True)
                ),
                density=normalize,
                norm=norm,
            )

        ax.set_xlabel(title.split(';')[1],family=self.props['axis_font'])
        ax.set_ylabel(title.split(';')[2],family=self.props['axis_font'])
        colorbar = fig.colorbar(c,ax=ax)
        colorbar.set_label(title.split(';')[3],family=self.props['axis_font'])

        ax.set_facecolor(self.props['facecolor'])
        fig.set_facecolor(self.props['background_color'])

        self.PltTextBox(ax)
        self.SavePltOutput(name)
        plt.close()
        return

    def Root2Plt_hist2d(self,h,grid=False,zmin=None,zmax=None,zlog=False,name='plot'):
        """
        TODO: Do not use this for now, it is broken. Matplotlib histogram interface is terrible!
        Given a ROOT histogram (TH2), draw it using matplotlib.
        Note: We will rebin things since mpl's hist2d seems incredibly clunky and memory-intensive, at least in this use.
        (Why is it so difficult in matplotlib to go from a 2D array of bin weights to a 2D histogram??????!!!)
        """
        fig,ax = plt.subplots(1,1)

        nbins_x = h.GetNbinsX()
        nbins_y = h.GetNbinsY()

        x_vals = np.array([[h.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins_x)] for j in range(nbins_y)])
        y_vals = np.array([[h.GetYaxis().GetBinLowEdge(i+1) for i in range(nbins_x)] for j in range(nbins_y)])
        weights = np.array([[h.GetBinContent(i+1,j+1) for i in range(nbins_x)] for j in range(nbins_y)])

        xbins = np.array([h.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins_x)] + [h.GetXaxis().GetBinLowEdge(nbins_x) + h.GetXaxis().GetBinWidth(nbins_x)])
        ybins = np.array([h.GetYaxis().GetBinLowEdge(i+1) for i in range(nbins_y)] + [h.GetYaxis().GetBinLowEdge(nbins_y) + h.GetYaxis().GetBinWidth(nbins_y)])
        # c = ax.imshow(weights,) # not what we want, doesn't do binning correctly (no use if bin edge info!)
        # _,_,_,c = ax.hist2d(x_vals.flatten(),y_vals.flatten(),bins=(xbins,ybins),weights=weights.flatten(),cmap='viridis') # doesn't work, not sure why not?
        ## pcolormesh *also* didn't work, gives no output
        # fig.colorbar(c,ax=ax)

        ax.set_xlabel(h.GetXaxis().GetTitle(),family=self.props['axis_font'])
        ax.set_ylabel(h.GetYaxis().GetTitle(),family=self.props['axis_font'])
        # ax.set_zlabel(h.GetZaxis().GetTitle(),family=self.props['axis_font']) # doesn't work -- why not???

        # if(zmin is None): zmin = 1.0e-3
        # if(zmax is None): zmax = 1.0
        # ax.set_zlim((zmin,zmax)) # doesn't work -- why doesn't z-axis follow same conventions as x and y????

        # if(zlog): ax.set_zscale('log') # doesn't work
        if(grid): ax.grid()

        ax.set_facecolor(self.props['facecolor'])
        fig.set_facecolor(self.props['background_color'])

        self.PltTextBox(ax)
        self.SavePltOutput(name)
        plt.close()
        return

    def Plt2RootRelabeling(self,h):
        if('TH1' in h.ClassName()):
            axes = [h.GetXaxis(),h.GetYaxis()]
        else:
            axes = [h.GetXaxis(),h.GetYaxis(),h.GetZaxis()]

        for axis in axes:
            axname = axis.GetTitle()
            rootname =  axname.replace('$','').replace('\\','#') # should work for most cases
            axis.SetTitle(rootname)
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

    def PltTextBox(self,ax):
        """
        Given a matplotlib axis object, draws text boxes.
        """
        for entry in self.text_boxes:
            x1,y1,x2,y2,text_lines,loc,alpha = entry
            if('right' in loc): x1 = x2
            text = '\n'.join(text_lines)
            text_box = AnchoredText(
                text,
                frameon=True,
                loc=loc,
                pad=0.5,prop={'fontsize':6},
                bbox_to_anchor=(x1,1.),
                bbox_transform=ax.transAxes
                )
            plt.setp(text_box.patch, facecolor=self.props['facecolor'], alpha=alpha)
            ax.add_artist(text_box)
        return

    def RootTextBox(self,canv):
        canv.cd()
        for entry in self.text_boxes:
            x1,y1,x2,y2,text_lines,loc,alpha = entry

            # The coordinates are optimized for matplotlib axis, we need to account
            # for them slightly differently here.
            # This is all a bit hacky for now...
            y1 = np.maximum(y2 - len(text_lines) * 0.025,0.9)
            x1 = np.maximum(x1,0.15)
            x2 = np.minimum(x2,0.9)

            pave = rt.TPaveText(x1,y1,x2,y2,'NDC')
            pave.SetBorderSize(0)
            pave.SetFillColorAlpha(rt.kWhite,alpha)
            pave.SetTextSize(0.075)
            pave.SetTextFont(82)
            pave.SetTextAlign(12)
            for line in text_lines:
                pave.AddText(line)
            pave.Draw()
            self.paves.append(pave)
        return

    def SavePltOutput(self,name,outdir=None):
        if(outdir is None): outdir = self.outdir
        for ext in self.plt_extensions:
            plt.savefig('{}/{}.{}'.format(outdir,name,ext))
        return

    def SaveRootOutputCanvas(self,h,c,name,ndim,draw_option=None):
        if(c is None): c = rt.TCanvas(RN(),'',800,600)
        self.canvases.append(c)
        c.cd()
        if(draw_option is None):
            if('TH1' in h.ClassName()): draw_option = 'HIST'
            else:  draw_option = 'COLZ'

        h.Draw(draw_option)
        self.RootTextBox(c)

        if('TH2' in h.ClassName()):
            rt.gPad.SetLeftMargin(0.1)
            rt.gPad.SetRightMargin(0.175)

        c.Draw()

        if('TH1' in h.ClassName()):
            h.SetFillColor(rt.kBlue)
            # rt.gPad.SetLogy()

        if('TH2' in h.ClassName()):
            rt.gPad.SetLogz()
            h.GetZaxis()
            h.GetZaxis().SetTitleOffset(0.95)

        if(ndim == 1): filename = self.root_filename_canvas_1d
        else: filename = self.root_filename_canvas_2d
        f = rt.TFile('{}/{}'.format(self.outdir,filename),'UPDATE')
        c.Write('canvas_' + name)
        f.Close()
        # self.SaveRootOutput(self.canvases[-1],'canvas_{}'.format(name),ndim=ndim) # TODO: Canvases are empty/undefined???
        return

    def SaveRootOutput(self,h,name,outdir=None,filename=None,ndim=1):
        if(outdir is None): outdir = self.outdir
        if(filename is None):
            if(ndim == 1): filename = self.root_filename_1d
            else: filename = self.root_filename_2d
        rootfile = '{}/{}'.format(outdir,filename)
        f = rt.TFile(rootfile,'UPDATE')
        h.Write(name)
        # print(h,name,h.ClassName())
        f.Close()

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
        'n':(200,0.,200.),
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
            'jet_dr' : DeltaR(vec_cyl[:,1:3],jet_Pmu_cyl[:,1:3])
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

    return

if(__name__ == '__main__'):
    main(sys.argv)
