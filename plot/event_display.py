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

def JetCurveInThetaEdges(eta0=0,radius=0.8):
    edge_1 = 2. * np.arctan(np.exp(-radius + eta0)) - 0.5 * np.pi + 1.0e-12
    edge_2 = 2. * np.arctan(np.exp( radius + eta0)) - 0.5 * np.pi - 1.0e-12
    return (edge_1,edge_2)

def JetCurveInEtaEdges(eta0=0,radius=0.8):
    edge_1 = -radius + eta0 + 1.0e-12 # to deal with numerical instability
    edge_2 = radius + eta0 - 1.0e-12 # to deal with numerical instability
    return (edge_1,edge_2)

def JetCurveInTheta(theta,phi0=0.,eta0=0.,radius=0.8):
    """
    While this is a circle in (eta,phi), it will be distorted
    in (theta,phi).
    """
    result = np.square(radius)
    b = np.square(np.log(-np.tan(0.5 * (-theta - 0.5 * np.pi))) - eta0)
    result -= b
    good_points = np.where(result >= 0)

    result_pos = np.full(result.shape,np.nan)
    result_neg = np.full(result.shape,np.nan)

    result_pos[good_points] =  np.sqrt(result[good_points]) + phi0
    result_neg[good_points] = -np.sqrt(result[good_points]) + phi0
    return (result_neg,result_pos)

def JetCurveinEta(eta,phi0=0.,eta0=0.,radius=0.8):
    """
    This is just a circle in (eta,phi).
    """
    result = np.square(radius)
    b = np.square(eta - eta0)
    result -= b
    good_points = np.where(result >= 0)

    result_pos = np.full(result.shape,np.nan)
    result_neg = np.full(result.shape,np.nan)

    result_pos[good_points] =  np.sqrt(result[good_points]) + phi0
    result_neg[good_points] = -np.sqrt(result[good_points]) + phi0
    return (result_neg,result_pos)

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
    output = '{}/{}'.format(outdir,args['output'])
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

    f.close()

    jet_Pmu = calculator.EPxPyPzToPtEtaPhiM_single(*jet_Pmu)
    jet_eta = jet_Pmu[1]
    jet_phi = jet_Pmu[2]
    jet_theta = np.nan

    if(use_theta): # TODO: Is this a good idea?
        jet_Pmu[1] = calculator.EtaToTheta(jet_Pmu[1])
        jet_theta = jet_Pmu[1]

    offset = np.array([0.,0.])
    offset_eta = 0.
    offset_phi = 0.
    offset_theta = 0.

    if(center):
        offset = jet_Pmu[1:3] # (eta,phi) of jet
        offset_eta = jet_eta
        offset_phi = jet_phi
        offset_theta = jet_theta

    dims = (800,800)
    c = rt.TCanvas(RN(),'',*dims)

    binning_eta = (10,-np.pi,np.pi)
    binning_phi = (10,-np.pi,np.pi)
    binning_theta = (10,-np.pi,np.pi)

    if(center):
        mult = 1.5
        binning_eta = (10,- mult *jet_radius, mult * jet_radius)
        binning_phi = (10,- mult *jet_radius, mult * jet_radius)
        binning_theta = (10,- mult *jet_radius, mult * jet_radius)

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
    if(center):
        scatter_Pmu.SetMarkerSize(1.2)

    # For truth_Pmu, we will make individual scatters for the different particle types,
    # based on the PDG ID. Above some particle index we will group them together -- can be used
    # for example for W boson daughter particles.
    max_unique = 5
    unique_pid = np.unique(truth_Pdg[:max_unique])

    truth_scatters = {}
    colors = [rt.kGreen + 1, rt.kOrange, rt.kViolet, rt.kRed, rt.kBlue] # TODO can this be done more nicely?
    styles = [rt.kPlus, rt.kPlus, rt.kPlus, rt.kCircle, rt.kCircle]

    for i,pid in enumerate(unique_pid):
        idxs = np.where(truth_Pdg == pid)[0]

        color = colors[idxs[0]%len(colors)]
        style = styles[idxs[0]%len(styles)]

        truth_scatters[pid] = rt.TGraph()
        truth_scatters[pid].SetMarkerColorAlpha(color,0.9)
        truth_scatters[pid].SetMarkerStyle(style)
        truth_scatters[pid].SetMarkerSize(1.)
        if(center):
            truth_scatters[pid].SetMarkerSize(2.)
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
    if(center):
        scatter_daughter.SetMarkerSize(1.0)

    dummy_hist.Draw('COL')
    scatter_Pmu.Draw('P SAME')

    scatter_daughter.Draw('P SAME')

    for pid,scatter in truth_scatters.items():
        scatter.Draw('P SAME')

    # Draw the jet radius.
    npoints = 1000
    jet_circle_color = rt.kGreen + 1
    jet_circle_top = rt.TGraph()
    jet_circle_bottom = rt.TGraph()

    if(not use_theta):
        # We use TGraph instead of TEllipse so we don't get weird overlap effects near axes.
        eta_jet_edges = np.array(JetCurveInEtaEdges(jet_eta,jet_radius))
        eta_jet_edges.sort()
        eta_carrier = np.zeros(npoints + 2)
        eta_carrier[0] = eta_jet_edges[0]
        eta_carrier[-1] = eta_jet_edges[1]
        eta_carrier[1:npoints+1] = np.linspace(-np.pi, np.pi,npoints)
        eta_carrier = np.sort(eta_carrier)

        curves = JetCurveinEta(eta_carrier,jet_phi,jet_eta,radius=jet_radius)
        curves = [c - offset_phi for c in curves]

        for i,circ in enumerate((jet_circle_bottom,jet_circle_top)):
            circ.SetLineColor(jet_circle_color)
            curve = curves[i]
            idxs = np.where(np.isnan(curve) == False)[0]
            x = eta_carrier[idxs] - offset_eta
            y = curve[idxs]

            for point in zip(x,y):
                circ.AddPoint(*point)

            circ.Draw('SAME L')

    else: # if using theta, this gets a bit more complex to visualize
        theta_jet_edges = np.array(JetCurveInThetaEdges(jet_eta,jet_radius))
        theta_jet_edges.sort()
        theta_carrier = np.zeros(npoints + 2)
        theta_carrier[0] = theta_jet_edges[0]
        theta_carrier[-1] = theta_jet_edges[1]
        theta_carrier[1:npoints+1] = np.linspace(binning_theta[1],binning_theta[2],npoints)
        theta_carrier = np.sort(theta_carrier)

        curves = JetCurveInTheta(theta_carrier,jet_phi,jet_eta,radius=jet_radius)
        curves = [c - offset_phi for c in curves]

        # for entry in zip(theta_carrier,curves[0]):
        #     print(entry)

        for i,circ in enumerate((jet_circle_bottom,jet_circle_top)):
            circ.SetLineColor(jet_circle_color)
            curve = curves[i]
            idxs = np.where(np.isnan(curve) == False)[0]
            x = theta_carrier[idxs] - offset_theta
            y = curve[idxs]

            for point in zip(x,y):
                circ.AddPoint(*point)

            circ.Draw('SAME L')

    rt.gPad.SetGrid()
    c.Draw()
    c.SaveAs(output)


if(__name__=='__main__'):
    main(sys.argv)