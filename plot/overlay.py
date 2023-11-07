import sys,os,glob
import numpy as np
import ROOT as rt
import argparse as ap
from plot_util.plot_util import RN

path_prefix = os.path.dirname(os.path.abspath(__file__)) + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util.qol_utils.qol_util import PlotStyle

def SetLogRange(hist):
    maximum = hist.GetMaximum()
    minimum = hist.GetMinimum()

    visual_maximum = np.power(10.,np.ceil(np.log10(maximum)))
    visual_minimum = np.power(10.,np.floor(np.log10(minimum)))

    hist.SetMinimum(visual_minimum)
    hist.SetMaximum(visual_maximum)


def main(args):

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputs',type=str,help='Semicolon-separated list of filename (or glob-compatible strings).')
    parser.add_argument('-t','--textfile',type=str,required=True,help='Input plaintext file listing canvas names to combine.')
    parser.add_argument('-o','--outdir',type=str,default='overlay_output',help='Output directory.')
    parser.add_argument('-m','--mode',type=str,default='light',help='Plot style mode. Options are currently "light" and "dark".')
    parser.add_argument('-norm','--normalize',type=int,default=0,help='Whether or not to normalize each distribution to integrate to 1.')
    parser.add_argument('-log','--log',type=int,default=1,help='Log for y axis.')
    parser.add_argument('-draw_option','--draw_option',type=str,default='NOSTACK HIST',help='Draw option.')
    parser.add_argument('-legend','--legend',type=str,default=None,help='Legend entries (semicolon-separated).')
    parser.add_argument('-legend_header','--legend_header',type=str,default=None,help='Header lines for legend (semicolon-separated).')
    parser.add_argument('-metadata','--metadata',type=int,default=1,help='Whether or not to write the metadata at the top.')
    parser.add_argument('-d','--descriptor',type=str,default=None,help='Descriptor string -- written in a TPaveText on top left of plot.')

    args = vars(parser.parse_args())

    inputnames = args['inputs'].split(';')
    inputs = []
    for name in inputnames: inputs += glob.glob(name)
    textfile = args['textfile']
    outdir = args['outdir']
    mode = args['mode']
    normalize = args['normalize'] > 0
    logy = args['log'] > 0
    draw_option = args['draw_option']
    legend_entries = args['legend']
    legend_header = args['legend_header']
    metadata_flag = args['metadata'] > 0
    descriptor = args['descriptor']
    inputs.sort()

    canvas_names = []
    with open(textfile,'r') as f:
        canvas_names = [y for y in [x.strip('\n').strip() for x in f.readlines()] if y != '']
    os.makedirs(outdir,exist_ok=True)

    infiles = {key:rt.TFile(key,"READ") for key in inputs}
    nfiles = len(inputs)

    outfile_root_name = '{}/overlay.root'.format(outdir)
    # outfile_root = rt.TFile(outfile_root_name,"RECREATE")

    canvases = []
    stacks = []
    hists = []
    paves = []

    dims = (800,600)

    # --- Setting style ---
    # pad_fill_color = (rt.kWhite,0.0)
    # color_palette = rt.kCool
    # text_color = rt.kBlack
    # line_color = rt.kBlack
    style = PlotStyle(mode=mode)
    style.SetStyle()
    # --- end style ---
    # rt.gStyle.SetPalette(color_palette)

    # Determine the space at the top of the plot, for all the metadata from the various files.
    top_margin = 0.05
    right_margin = 0.1
    left_margin = 0.1
    if(metadata_flag): top_margin = np.minimum(0.05 * nfiles,0.3)
    pave_vertical_width = top_margin / nfiles

    for i,canvas_name in enumerate(canvas_names):
        c = rt.TCanvas(RN(),'',*dims)
        stack = rt.THStack()
        canvas_paves = []

        for j,(key,file) in enumerate(infiles.items()):
            # print(key,file)
            c_tmp = file.Get(canvas_name)
            primitives = c_tmp.GetListOfPrimitives() # will assume a certain ordering
            h_tmp = primitives.At(1).Clone()

            try:
                nevents = file.Get('nevents').GetVal()
            except:
                print('{}:\tFailed to get nevents.'.format(key))
                nevents = 1

            if(not normalize):
                h_tmp.Scale(nevents)

            paves_tmp = [primitives.At(x) for x in (2,3,4)]
            for pave in paves_tmp:
                pave.SetY2NDC(1. - j * pave_vertical_width)
                pave.SetY1NDC(1. - (j+1) * pave_vertical_width)
                pave.SetTextSize(0.0175)
                pave.SetFillColorAlpha(rt.kWhite,0.)
                # pave.SetFillColorAlpha(rt.kGray, 0.2)
                canvas_paves.append(pave)

            h_tmp.SetFillColorAlpha(h_tmp.GetFillColor(),0.)
            h_tmp.SetLineWidth(2)
            stack.Add(h_tmp)
            hists.append(h_tmp)

        c.cd()
        if(logy): rt.gPad.SetLogy()
        rt.gPad.SetTopMargin(top_margin)
        rt.gPad.SetRightMargin(right_margin)
        rt.gPad.SetLeftMargin(left_margin)
        # rt.gPad.SetFillColorAlpha(*pad_fill_color)

        if('PLC' not in draw_option):
            draw_option += ' PLC'
        stack.Draw(draw_option)
        if(logy and normalize):
            SetLogRange(stack.GetHistogram())
            # if(stack.GetHistogram().GetMaximum() < 0.1):
            #     stack.GetHistogram().SetMaximum(0.1)
            #     stack.GetHistogram().SetMinimum(0.001)
            # else:
            #     stack.GetHistogram().SetMaximum(1.)
            #     stack.GetHistogram().SetMinimum(0.01)
        elif(not logy):
            stack.GetHistogram().SetMaximum(1.1 * stack.GetMaximum())
            stack.GetHistogram().SetMinimum(0.)

        stack.SetMaximum(1.1 * stack.GetMaximum())

        if(legend_entries is not None):
            legend_entries_list = legend_entries.split(';')
            nentries = len(legend_entries_list)
            if(legend_header is not None):
                legend_header_lines = legend_header.split(';')
                nentries += len(legend_header_lines)
            x2 = 1. - right_margin
            y2 = 1. - top_margin
            x1 = x2 - 0.05 - np.minimum(0.6,0.075 * nentries)
            y1 = y2 - np.minimum(0.6,0.075 * nentries)
            legend = rt.TLegend(x1,y1,x2,y2)

            if(legend_header is not None):
                legend_header_lines = legend_header.split(';')
                for line in legend_header_lines:
                    legend.AddEntry(0,line,'')

            for j,entry in enumerate(legend_entries_list):
                legend.AddEntry(stack.GetHists().At(j),entry,'l')
            legend.SetBorderSize(0)
            legend.SetFillColorAlpha(style.canv,0.1)
            legend.SetTextSize(0.04)
            legend.Draw()

        if(descriptor is not None):
            x1 = left_margin + 0.05
            y2 = 1.0 - top_margin - 0.01
            x2 = x1 + len(descriptor) * 0.015
            y1 = y2 - 0.05
            pave = rt.TPaveText(x1,y1,x2,y2,'NDC')
            pave.AddText(descriptor)
            pave.SetTextAlign(31)
            pave.SetTextFont(82)
            pave.SetBorderSize(0)
            pave.SetFillColorAlpha(style.canv,0.1)
            pave.Draw()

        color_indices = np.linspace(0,254,int(len(canvas_paves)/3),dtype=int,endpoint=True)
        # print(color_indices)
        colors = [rt.TColor.GetPalette().At(int(x)) for x in color_indices]

        if(metadata_flag):
            for j,pave in enumerate(canvas_paves):
                # print(j,len(colors), int(np.floor(j/2)))
                color = colors[int(np.floor(j/3))]
                pave.SetTextColor(color)
                pave.Draw()
        stack.GetXaxis().SetTitle(hists[-1].GetXaxis().GetTitle())
        if(normalize):
            stack.GetYaxis().SetTitle('Fractional Count')
        else:
            stack.GetYaxis().SetTitle('Count')
        # stacks.append(stack)
        c.Draw()
        extensions = ['pdf','root']
        for ext in extensions:
            c.SaveAs('{}/{}.{}'.format(outdir,canvas_name,ext))

    for key,file in infiles.items():
        file.Close()

if(__name__ == '__main__'):
    main(sys.argv)