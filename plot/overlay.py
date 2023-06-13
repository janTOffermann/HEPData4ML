import sys,os,glob
import numpy as np
import ROOT as rt
import argparse as ap
from plot_util.plot_util import RN

path_prefix = os.path.dirname(os.path.abspath(__file__)) + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util.qol_utils.qol_util import PlotStyle

def main(args):

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputs',type=str,help='Semicolon-separated list of filename (or glob-compatible strings).')
    parser.add_argument('-t','--textfile',type=str,required=True,help='Input plaintext file listing canvas names to combine.')
    parser.add_argument('-o','--outdir',type=str,default='overlay_output',help='Output directory.')
    parser.add_argument('-m','--mode',type=str,default='light',help='Plot style mode. Options are currently "light" and "dark".')
    args = vars(parser.parse_args())

    inputnames = args['inputs'].split(';')
    inputs = []
    for name in inputnames: inputs += glob.glob(name)
    textfile = args['textfile']
    outdir = args['outdir']
    inputs.sort()

    canvas_names = []
    with open(textfile,'r') as f:
        canvas_names = [y for y in [x.strip('\n').strip() for x in f.readlines()] if y != '']
    os.makedirs(outdir,exist_ok=True)

    infiles = {key:rt.TFile(key,"READ") for key in inputs}
    nfiles = len(inputs)

    outfile_root_name = '{}/overlay.root'.format(outdir)
    outfile_root = rt.TFile(outfile_root_name,"RECREATE")

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
    style = PlotStyle(mode='dark')
    style.SetStyle()
    # --- end style ---
    # rt.gStyle.SetPalette(color_palette)

    # Determine the space at the top of the plot, for all the metadata from the various files.
    top_margin = np.minimum(0.05 * nfiles,0.3)
    pave_vertical_width = top_margin / nfiles

    for i,canvas_name in enumerate(canvas_names):
        c = rt.TCanvas(RN(),'',*dims)
        stack = rt.THStack()
        canvas_paves = []

        for j,(key,file) in enumerate(infiles.items()):
            c_tmp = file.Get(canvas_name)
            primitives = c_tmp.GetListOfPrimitives() # will assume a certain ordering
            h_tmp = primitives.At(1).Clone()

            paves_tmp = [primitives.At(x) for x in (2,3,4)]
            for pave in paves_tmp:
                pave.SetY2NDC(1. - j * pave_vertical_width)
                pave.SetY1NDC(1. - (j+1) * pave_vertical_width)
                pave.SetTextSize(0.01)
                pave.SetFillColorAlpha(rt.kWhite,0.)
                # pave.SetFillColorAlpha(rt.kGray, 0.2)
                canvas_paves.append(pave)

            h_tmp.SetFillColorAlpha(h_tmp.GetFillColor(),0.)
            h_tmp.SetLineWidth(2)
            stack.Add(h_tmp)
            hists.append(h_tmp)

        c.cd()
        rt.gPad.SetLogy()
        rt.gPad.SetTopMargin(top_margin)
        # rt.gPad.SetFillColorAlpha(*pad_fill_color)
        stack.Draw('NOSTACK HIST PLC')

        color_indices = np.linspace(0,254,int(len(canvas_paves)/3),dtype=int,endpoint=True)
        # print(color_indices)
        colors = [rt.TColor.GetPalette().At(int(x)) for x in color_indices]

        for j,pave in enumerate(canvas_paves):
            # print(j,len(colors), int(np.floor(j/2)))
            color = colors[int(np.floor(j/3))]
            pave.SetTextColor(color)
            pave.Draw()
        stack.GetXaxis().SetTitle(hists[-1].GetXaxis().GetTitle())
        stack.GetYaxis().SetTitle(hists[-1].GetYaxis().GetTitle())
        # stacks.append(stack)
        c.Draw()
        c.SaveAs('{}.pdf'.format(canvas_name))

    for key,file in infiles.items():
        file.Close()







if(__name__ == '__main__'):
    main(sys.argv)