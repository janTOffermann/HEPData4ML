# For plotting kinematics of the dataset.

import sys, os
import numpy as np
import h5py as h5
import ROOT as rt
import argparse as ap
# custom imports
path_prefix = os.getcwd() + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util import qol_util as qu
from util.calcs import PxPyPzEToPtEtaPhiM
from plot_util.plot_util import *

def main(args):

    # Parse arguments.
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input directory.', required=True)
    parser.add_argument('-I', '--infile',  type=str,  help='Input filename.', default='events.h5')
    parser.add_argument('-o', '--output', type=str, help='Output ROOT filename.', default = 'events.root')
    parser.add_argument('-l', '--legend', type=str, help="Name for legend.", default='')
    args = vars(parser.parse_args())
    indir = args['input']
    infile = args['infile']
    title = args['legend']
    outfile = args['output']

    # Set some ROOT options.
    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptStat(0)

    colors = {
        'constituents': rt.kRed,
        'top': rt.kGreen,
        'bottom': rt.kOrange,
        'W': rt.kBlue,
        'jet': rt.kViolet
    }

    filename = '{}/{}'.format(indir,infile)
    f = h5.File(filename,'r')

    top = f['truth_Pdg'][:] == 6
    bottom = f['truth_Pdg'][:] == 5
    W = f['truth_Pdg'][:] == 24

    # Number of constituents per event
    nobj = f['Nobj'][:]

    # Constituent momenta (throwing away the zero-padding).
    pmu = GetJagged(f['Pmu'][:],f['Nobj'][:])

    # particle energy -- component 0
    E = {
        'constituents': pmu[:,0].flatten(),
        'top':f['truth_Pmu'][:][top][:,0].flatten(),
        'bottom':f['truth_Pmu'][:][bottom][:,0].flatten(),
        'W':f['truth_Pmu'][:][W][:,0].flatten(),
        'jet': f['jet_Pmu'][:][:,0].flatten()
    }

    # particle px -- component 1
    px = {
        'constituents': pmu[:,1].flatten(),
        'top':f['truth_Pmu'][:][top][:,1].flatten(),
        'bottom':f['truth_Pmu'][:][bottom][:,1].flatten(),
        'W':f['truth_Pmu'][:][W][:,1].flatten(),
        'jet':f['jet_Pmu'][:][:,1].flatten()
    }

    # particle py -- component 2
    py = {
        'constituents': pmu[:,2].flatten(),
        'top':f['truth_Pmu'][:][top][:,2].flatten(),
        'bottom':f['truth_Pmu'][:][bottom][:,2].flatten(),
        'W':f['truth_Pmu'][:][W][:,2].flatten(),
        'jet':f['jet_Pmu'][:][:,2].flatten()
    }

    # particle pz -- component 3
    pz = {
        'constituents': pmu[:,3].flatten(),
        'top':f['truth_Pmu'][:][top][:,3].flatten(),
        'bottom':f['truth_Pmu'][:][bottom][:,3].flatten(),
        'W':f['truth_Pmu'][:][W][:,3].flatten(),
        'jet':f['jet_Pmu'][:][:,3].flatten()
    }

    # We can also get particles' momenta in (pt,eta,phi,m)
    p_cyl = {key: PxPyPzEToPtEtaPhiM(px[key],py[key],pz[key],E[key]) for key in E.keys()}
    pt = {key: p_cyl[key][:,0] for key in E.keys()}
    eta = {key: p_cyl[key][:,1] for key in E.keys()}
    phi = {key: p_cyl[key][:,2] for key in E.keys()}
    m = {key: p_cyl[key][:,3] for key in E.keys()}

    f.close()

    # Let's also get the ∆R between the jets and the truth-level particles (top is the most interesting).
    dr = {
        key: DeltaR(p_cyl[key][:,1:3], p_cyl['jet'][:,1:3])
        for key in ['top','bottom', 'W'] # TODO: Order does matter when referencing in other scripts!
    }

    canvases = []
    stacks = []
    hists = []
    legends = []

    f = rt.TFile(outfile,'RECREATE')
    c, leg, hlist, hstack = KinematicDraw('constituents',nobj, (200,0.,200.), title + ';Number of jet constituents;Fractional Count', False, False, colors=colors)
    canvases.append(c)
    legends.append(leg)
    stacks.append(hstack)
    for h in hlist: hists.append(h)

    c, leg, hlist, hstack = KinematicDraw(E.keys(),pt, (400,0.,4000.), title + ';p_{T} [GeV];Fractional Count', True, True, colors=colors)
    canvases.append(c)
    legends.append(leg)
    stacks.append(hstack)
    for h in hlist: hists.append(h)

    c, leg, hlist, hstack = KinematicDraw(E.keys(),E, (200,0.,4000.), title + ';E [GeV];Fractional Count', True, True, colors=colors)
    canvases.append(c)
    legends.append(leg)
    stacks.append(hstack)
    for h in hlist: hists.append(h)

    c, leg, hlist, hstack = KinematicDraw(E.keys(),m, (125,0.,500.), title + ';m [GeV];Fractional Count', False, True, colors=colors)
    canvases.append(c)
    legends.append(leg)
    stacks.append(hstack)
    for h in hlist: hists.append(h)

    c, leg, hlist, hstack = KinematicDraw(dr.keys(),dr, (200,0.,2.), title + ';#Delta R_{jet};Fractional Count', True, True, colors=colors)
    canvases.append(c)
    legends.append(leg)
    stacks.append(hstack)
    for h in hlist: hists.append(h)

    dr_title = title + ';p_{T} (top);#Delta R (top,jet);Fractional Count'
    c_dr, h_dr = KinematicDraw2D(pt['top'], dr['top'], (120,540.,660.,500,1e-3,.5), dr_title, False, True)

    s = rt.THStack(qu.RN(),dr_title)
    s.Add(h_dr)
    canvases.append(c_dr)
    stacks.append(s)

    for i,c in enumerate(canvases): c.Write('c_{}'.format(i))
    for i,s in enumerate(stacks): s.Write('s_{}'.format(i)) # TODO: Remove the histogram stacks, and just save the histograms.
    f.Close()

if __name__ == '__main__':
    main(sys.argv)