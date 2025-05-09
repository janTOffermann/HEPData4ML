import uuid,os
import numpy as np
import h5py as h5
import ROOT as rt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.ndimage import gaussian_filter

def RN():
    return str(uuid.uuid4())

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
        self.nevents = -1

        self.canvases = [] # for persistence
        self.paves = []

        for filename in [self.root_filename_1d,self.root_filename_2d,self.root_filename_canvas_1d,self.root_filename_canvas_2d]:
            try:
                os.unlink('{}/{}'.format(outdir,filename))
            except:
                pass

        self.text_boxes = []

    def SetNevents(self,val):
        self.nevents = val

    def SaveNevents(self):
        for filename in [self.root_filename_1d,self.root_filename_2d,self.root_filename_canvas_1d,self.root_filename_canvas_2d]:
            filename = '{}/{}'.format(self.outdir,filename)
            f = rt.TFile(filename,'UPDATE')
            nevents_param = rt.TParameter(int)('nevents',self.nevents)
            nevents_param.Write('nevents')
            f.Close()

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
        if(self.nevents < 0): self.nevents = len(data)
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
        if(self.nevents < 0): self.nevents = len(data_x)
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

        if(ymin is None): ymin = 0.
        if(ymax is None): ymax = 1.0
        ax.set_ylim((ymin,ymax))

        if(ylog): ax.set_yscale('symlog',linthresh=1.0e-6)
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

        norm = mpl.colors.SymLogNorm(linthresh=1.0e-4,vmin=0.,vmax=3.0e-2)
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
            c = ax.imshow(heatmap.T,extent=extent,origin='lower', cmap=plt.cm.jet, norm=norm)
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

        # xbins = np.array([h.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins_x)] + [h.GetXaxis().GetBinLowEdge(nbins_x) + h.GetXaxis().GetBinWidth(nbins_x)])
        # ybins = np.array([h.GetYaxis().GetBinLowEdge(i+1) for i in range(nbins_y)] + [h.GetYaxis().GetBinLowEdge(nbins_y) + h.GetYaxis().GetBinWidth(nbins_y)])
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
