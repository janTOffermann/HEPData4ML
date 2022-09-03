import sys, os
import numpy as np
import ROOT as rt

# custom imports
path_prefix = os.getcwd() + '/../../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
from util import qol_util as qu
# from util.calcs import PxPyPzEToPtEtaPhiM, DeltaR2

# Helper function for fetching four-momenta from a jagged array.
def GetJagged(pmu, nobj):
    nentries_final = np.sum(nobj)
    result = np.zeros((nentries_final,pmu.shape[-1]))
    nentries = pmu.shape[0]
    counter = 0
    for i in range(nentries):
        result[counter:counter+nobj[i],:] = pmu[i,:nobj[i],:]
        counter += nobj[i]
    return result

# Helper function for coordinate conversion.
def PxPyPzEToPtEtaPhiM(px,py,pz,e):
    nvec = len(px)
    result = np.zeros((nvec,4))
    lvec = rt.Math.PxPyPzEVector()
    for i in range(nvec):
        lvec.SetCoordinates(px[i],py[i],pz[i],e[i])
        result[i,:] = np.array([lvec.Pt(),lvec.Eta(),lvec.Phi(),lvec.M()])
    return result

# Helper function for calculating dR between two sets pf vectors.
# Here, vecs are of form (eta, phi), one entry per event.
def DeltaR(vec1, vec2):
    result = np.zeros(vec1.shape[0])
    lvec1 = rt.Math.PtEtaPhiMVector()
    lvec2 = rt.Math.PtEtaPhiMVector()
    for i in range(len(result)):

        lvec1.SetCoordinates(0.,vec1[i,0],vec1[i,1],0.)
        lvec2.SetCoordinates(0.,vec2[i,0],vec2[i,1],0.)
        result[i] = np.sqrt(rt.Math.VectorUtil.DeltaR2(lvec1,lvec2))

        # result[i] = np.sqrt(DeltaR2(*vec1[i], *vec2[i]))
    return result

# Function for making a kinematic plot (1D).
def KinematicDraw(keys, data, bin_info, title, logx=True, logy=True, colors=None):
    c = rt.TCanvas(qu.RN(),'c0',800,600)
    legend = rt.TLegend(0.8,0.8,0.95,0.95)
    hstack = rt.THStack(qu.RN(),title)
    hlist = []
    nbins, xmin, xmax = bin_info

    dummy = {'a':'a'}

    if(type(keys) not in (list, type(dummy.keys()))):
        keys = [keys]
        data = {keys[0]:data}

    for key in keys:

        h = rt.TH1F(qu.RN(),'',nbins,xmin,xmax)
        for elem in data[key]: h.Fill(elem)
        if(colors is not None): h.SetLineColor(colors[key])
        integral = h.Integral()
        if(integral == 0): integral = 1.
        h.Scale(1./integral)
        hstack.Add(h)
        legend.AddEntry(h,key,'l')
        hlist.append(h)

    draw_option = "NOSTACK HIST"
    if(colors is None): draw_option += " PLC PMC"
    hstack.Draw(draw_option)
    legend.Draw()
    if(logx): c.SetLogx()
    if(logy): c.SetLogy()
    c.Draw()
    return c, legend, hlist, hstack

# Function for making a kinematic plot (2D).
def KinematicDraw2D(data_x, data_y, bin_info, title, logx=True, logy=True):
    c = rt.TCanvas(qu.RN(),'c0',800,600)
    # legend = rt.TLegend(0.8,0.8,0.95,0.95)
    nx, xmin, xmax, ny, ymin, ymax = bin_info
    # dummy = {'a':'a'}

    h = rt.TH2F(qu.RN(),title,nx, xmin, xmax, ny, ymin, ymax)

    for i in range(len(data_x)):
        h.Fill(data_x[i], data_y[i])

    integral = h.Integral()
    if(integral == 0): integral = 1.
    h.Scale(1./integral)
    h.Draw('COLZ')
    if(logx): c.SetLogx()
    if(logy): c.SetLogy()
    c.SetLogz()
    c.SetRightMargin(0.15)
    c.Draw()
    return c, h
