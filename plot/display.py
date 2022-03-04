import sys, os, glob
import ROOT as rt
import numpy as np
import argparse as ap
path_prefix = os.getcwd() + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util import qol_util as qu

def main(args):

    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input',   nargs='+', help='Input file pattern (or a list of patterns).', required=True)
    parser.add_argument('-o', '--outdir',  type=str,  help='Output directory.', default=None)
    parser.add_argument('-l', '--legend', nargs='+', help='Legend entries for data series. Must be one per input file.', default=None)
    parser.add_argument('-d', '--draw_option', type=str, help="Draw option for 1D histograms.", default='')
    parser.add_argument('-m', '--mode',type=str,help="Drawing mode. (Try \'light\' or \'dark\').",default='dark')
    args = vars(parser.parse_args())

    # Set some ROOT & drawing options.
    rt.gStyle.SetOptStat(0)
    rt.gROOT.SetBatch(True)
    style = args['mode']
    ps = qu.PlotStyle(style)
    ps.SetStyle()

    # Determine what are the input files.
    infiles = []
    input_patterns = args['input']
    for pattern in input_patterns:
        if(',') in pattern:
            infiles += [x.strip() for x in pattern.split(',')]
        elif('*') in pattern:
            infiles += glob.glob(os.path.expanduser(pattern),recursive=True) # The "recursive" argument allows us to use patterns containing things like "**", which is not supported in Python2's glob. The os.path.expanduser func deals with any tilde expansion.
        else:
            infiles += [pattern]

    leg_entries = args['legend']
    if(leg_entries is None): leg_entries = ['No MPI/ISR/FSR','With ISR/FSR','With MPI/ISR/FSR']
    assert(len(leg_entries) >= len(infiles))

    # Prepare output directory
    outdir = args['outdir']
    if(outdir is not None): os.makedirs(outdir, exist_ok=True)
    else: outdir = os.getcwd()

    # sort infiles to be alphabetic
    infiles.sort()
    infiles = [rt.TFile(x,'READ') for x in infiles]

    variables = ['nobj','pt','E','m','dr', 'dr2d'] # bookkeeping
    nvars = len(variables)
    d2 = np.full(nvars,False,dtype=bool)
    d2[-1] = True

    # colors = [rt.kRed, rt.kViolet, rt.kBlue]
    # if(style == 'dark'): colors = [rt.kRed, ps.curve, rt.kSpring]
    draw_option = args['draw_option']
    draw_option = 'SAME HIST {} PLC'.format(draw_option) #HIST removes error bars

    titles = []
    for i in range(nvars): titles.append(';' + ';'.join(infiles[0].Get('s_{}'.format(i)).GetTitle().split(';')[1:]))
    particles = ['constituents','top','bottom','W^{+}','jet']
    npars = len(particles)

    logx = [False, True, True, False, False, False]
    logy = [False, True, True, True, True, True]
    assert(len(logx) == len(logy) and len(logx) == nvars)

    plot_names = [
        'nconsts',
        'pt',
        'e',
        'm',
        'dr',
        'dr2d'
    ]

    # special naming for "unique variables" like nobj, dr (not plotted for each particle species)
    spec_names = {
        0: 'constituents' + titles[0], # nobj only plotted for constituents
        5: 'top' # dR(x,jet) only plotted for x==top (this is a 2D plot)
    }

    ranges = [
        [[0,100]], # number of constituents,
        [[1.0e-2, 1000.], # consts
        [90.,1000.],# top
        [10., 10000.], # bottom
        [1.,2.0e3], # W
        [1.,1000.]# jet
        ], # pT,
        [[1.0e-2, 1000.], # consts
        [90.,5000.],# top
        [10., 2e4], # bottom
        [10.,5000.], # W
        [1.,5000.]# jet
        ], # Energy,
        [[0, 10.], # consts
        [0.,400.], # top
        [0, 10.], # bottom
        [0.,200.], # W
        [0.,400.] # jet
        ], # mass
        [[0.,2.], # top
        [0, 2.], # bottom
        [0.,2.], # W
        ] # dR
    ]

    canvases = []
    hists = []
    legends = []
    for i in range(nvars):
        c = rt.TCanvas(qu.RN(),'c_{}'.format(i),1600,1200)
        leg = rt.TLegend(0.7,0.75,0.9,0.9)
        leg.SetFillColorAlpha(ps.canv, 0.5)
        hists = infiles[0].Get('s_{}'.format(i)).GetHists()

        K = len(hists)
        if(K in [3,4]): c.Divide(2,2)
        elif(K == 5): c.Divide(2,3)

        for j, rfile in enumerate(infiles):
            hists = rfile.Get('s_{}'.format(i)).GetHists()

            for k in range(K):
                if(K > 1): c.cd(k+1)
                else: c.cd()

                h = hists[k]
                if(K == 1): name = spec_names[i]
                elif(i == 4): name = particles[k+1] + titles[i] # for dR, skip constituents
                else: name = particles[k] + titles[i]
                h.SetTitle(name)

                if(not d2[i]):
                    h.Draw(draw_option)
                    h.SetLineWidth(2)
                    # h.SetLineColor(colors[j])

                else:
                    rt.gPad.SetRightMargin(0.125)
                    h.Draw('COLZ')
                    palette = h.GetListOfFunctions()[0]
                    palette.SetX1NDC(0.825)
                    palette.SetX2NDC(0.875)
                    #h.GetZaxis().SetLabelSize(0.02)
                    h.GetZaxis().SetMaxDigits(2)

                rt.gPad.SetLogx(logx[i])
                rt.gPad.SetLogy(logy[i])

                for ax in [h.GetXaxis(), h.GetYaxis(), h.GetZaxis()]:
                    ax.SetTitleColor(ps.text)
                    ax.SetLabelColor(ps.text)
                    ax.SetAxisColor(ps.main)
                    rt.gPad.SetFrameLineColor(ps.main)

                if(not d2[i]): h.GetXaxis().SetRangeUser(*(ranges[i][k]))
                hists.append(h)
                if(k==0): leg.AddEntry(h,leg_entries[j],'l')

        if(K == 1):
            if(not d2[i]): leg.Draw()
        else:
            for k in range(K):
                c.cd(k+1)
                leg.Draw()

        leg.SetTextColor(ps.text)
        legends.append(leg)
        canvases.append(c)

    for i,c in enumerate(canvases):
        c.Draw()
        c.SaveAs('{}/{}.png'.format(outdir,plot_names[i]))

    for x in infiles: x.Close()

if __name__ == '__main__':
    main(sys.argv)