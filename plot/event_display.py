import sys,os,glob,pathlib,uuid
import ROOT as rt
import numpy as np
import h5py as h5
import argparse as ap

# custom imports
path_prefix = os.path.dirname(os.path.abspath(__file__)) + '/../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util.qol_utils.qol_util import RN
from util.calcs import Calculator

def main(args):

    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptStat(False)
    rt.gROOT.SetStyle("ATLAS")

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--input',type=str,required=True,help='Input HDF5 file.')
    parser.add_argument('-index','--index',type=int,default=0,help='Event index to draw from.')
    parser.add_argument('-O','--outdir',type=str,default=None)
    parser.add_argument('-o','--output',type=str,default='event_display.pdf')
    parser.add_argument('-r','--jet_radius',type=float,default=0.8)
    parser.add_argument('-c','--center',type=int,default=0,help='Whether or not to translate images to center on jet.')
    parser.add_argument('-rotation','--rotation',type=int,default=0,help='Whether or not to plot rotated events.')
    parser.add_argument('-theta','--theta',type=int,default=0,help='If > 0, use theta instead of eta for plotting.')
    args = vars(parser.parse_args())

    input = args['input']
    outdir = args['outdir']
    output = args['output']
    event_index = args['index']
    jet_radius = args['jet_radius']
    center = args['center'] > 0
    rotation = args['rotation'] > 0
    use_theta = args['theta'] > 0

    if(outdir is None): outdir = os.getcwd()
    os.makedirs(outdir,exist_ok=True)

    # Initialize our calculator
    calculator = Calculator(use_vectorcalcs=False)

    # Open the file.

    key_suffix = ''
    if(rotation): key_suffix = '_rot'

    f = h5.File(input,'r')
    Nobj = f['Nobj'][event_index]
    Pmu = f['Pmu{}'.format(key_suffix)][event_index,:Nobj]
    truth_Nobj = f['truth_Nobj'][event_index]
    truth_Pmu = f['truth_Pmu{}'.format(key_suffix)][event_index,:truth_Nobj]
    truth_Pdg = f['truth_Pdg'][event_index,:truth_Nobj]
    jet_Pmu = f['jet_Pmu{}'.format(key_suffix)][event_index]

    # TODO: This is a temporary test.
    energy_ratio_truth = f['energy_ratio_truth'][event_index,:Nobj]
    daughters = Pmu[np.where(energy_ratio_truth == 1)[0][0]]
    # print(energy_ratio_truth)

    # keys = sorted(list(f.keys()))
    # for k in keys: print(k)
    f.close()

    jet_Pmu = calculator.EPxPyPzToPtEtaPhiM_single(*jet_Pmu)
    if(use_theta):
        jet_Pmu[1] = calculator.EtaToTheta(jet_Pmu[1])

    offset = np.array([0.,0.])
    if(center):
        offset = jet_Pmu[1:3] # (eta,phi) of jet

    dims = (800,800)
    c = rt.TCanvas(RN(),'',*dims)

    binning_eta = (10,-np.pi,np.pi)
    binning_phi = (10,-np.pi,np.pi)
    binning_theta = (10,-0.5 * np.pi, 0.5 * np.pi)

    if(center):
        mult = 1.25
        binning_eta = (10,- mult *jet_radius, mult * jet_radius)
        binning_phi = (10,- mult *jet_radius, mult * jet_radius)
        binning_dtheta = (10,-0.5 * np.pi, 0.5 * np.pi)

    histogram_name = 'Event Display;#eta;#phi'
    binning_x = binning_eta
    binning_y = binning_phi # in case we want to change things later

    if(center):
        histogram_name = 'Event Display;#Delta#eta;#Delta#phi'

    if(use_theta):
        histogram_name = histogram_name.replace('#eta','#theta')
        binning_x = binning_theta

    if(rotation):
        histogram_name = histogram_name.replace('#eta','#eta\'')
        histogram_name = histogram_name.replace('#theta','#theta\'')
        histogram_name = histogram_name.replace('#phi','#phi\'')

    dummy_hist = rt.TH2D(RN(),histogram_name,*binning_x,*binning_y)

    # TODO: Use the new TScatter class! Need to see where this is available...
    scatter_Pmu = rt.TGraph()
    for Pmu_vec in Pmu:
        Pmu_cyl = calculator.EPxPyPzToPtEtaPhiM_single(*Pmu_vec)
        if(use_theta):
            Pmu_cyl[1] = calculator.EtaToTheta(Pmu_cyl[1])
        Pmu_cyl[1:3] -= offset
        scatter_Pmu.AddPoint(*Pmu_cyl[1:3])
    scatter_Pmu.SetMarkerColor(rt.kBlack)
    scatter_Pmu.SetMarkerStyle(rt.kCircle)
    scatter_Pmu.SetMarkerSize(0.6)

    # For truth_Pmu, we will make individual scatters for the different particle types,
    # based on the PDG ID. Above some particle index we will group them together -- can be used
    # for example for W boson daughter particles.
    max_unique = 5
    unique_pid = np.unique(truth_Pdg[:max_unique])

    truth_scatters = {}
    colors = [rt.kGreen, rt.kOrange, rt.kViolet, rt.kRed, rt.kBlue] # TODO can this be done more nicely?
    styles = [rt.kDiamond, rt.kStar, rt.kStar, rt.kCircle, rt.kCircle]

    for i,pid in enumerate(unique_pid):
        idxs = np.where(truth_Pdg == pid)[0]

        color = colors[idxs[0]%len(colors)]
        style = styles[idxs[0]%len(styles)]

        truth_scatters[pid] = rt.TGraph()
        truth_scatters[pid].SetMarkerColorAlpha(color,0.7)
        truth_scatters[pid].SetMarkerStyle(style)
        truth_scatters[pid].SetMarkerSize(1.)
        vecs = calculator.PxPyPzEToPtEtaPhiM(*np.roll(truth_Pmu[idxs],-1,axis=-1).T)
        if(use_theta):
            vecs[:,1] = calculator.EtaToTheta(vecs[:,1])
        for vec in vecs:
            vec[1:3] -= offset
            truth_scatters[pid].AddPoint(*vec[1:3])

    scatter_daughter = rt.TGraph()
    vecs = calculator.PxPyPzEToPtEtaPhiM(*np.roll(truth_Pmu[max_unique:],-1,axis=-1).T)
    if(use_theta):
        vecs[:,1] = calculator.EtaToTheta(vecs[:,1])
    for vec in vecs:
        vec[1:3] -= offset
        scatter_daughter.AddPoint(*vec[1:3])

    scatter_daughter.SetMarkerColor(rt.kCyan)
    scatter_daughter.SetMarkerStyle(rt.kCircle)
    scatter_daughter.SetMarkerSize(0.4)

    scatter_daughter_tmp = rt.TGraph()
    print(daughters.shape)
    vecs = calculator.PxPyPzEToPtEtaPhiM(*np.roll(daughters,-1,axis=-1).T)
    if(use_theta):
        vecs[:,1] = calculator.EtaToTheta(vecs[:,1])
    for vec in vecs:
        vec[1:3] -= offset
        scatter_daughter_tmp.AddPoint(*vec[1:3])

    scatter_daughter_tmp.SetMarkerColor(rt.kViolet)
    scatter_daughter_tmp.SetMarkerStyle(rt.kCircle)
    scatter_daughter_tmp.SetMarkerSize(0.6)

    dummy_hist.Draw('COL')
    scatter_Pmu.Draw('P SAME')

    scatter_daughter.Draw('P SAME')
    scatter_daughter_tmp.Draw('P SAME')

    for pid,scatter in truth_scatters.items():
        scatter.Draw('P SAME')

    # Draw the jet radius.
    if(not use_theta):
        jet_Pmu[1:3] -= offset
        jet_circle = rt.TEllipse(*jet_Pmu[1:3],jet_radius,jet_radius)
        jet_circle.SetLineColor(rt.kSpring)
        jet_circle.SetFillColorAlpha(rt.kWhite,0.)
        jet_circle.Draw()

    else: # if using theta, this gets a bit more complex to visualize
        pass # TODO: make a curve to plot?

    rt.gPad.SetGrid()
    c.Draw()
    c.SaveAs(output)


if(__name__=='__main__'):
    main(sys.argv)