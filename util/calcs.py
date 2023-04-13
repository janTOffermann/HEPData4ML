# Note: Functions here use our custom VectorCalcs C++/ROOT library.
#       That library must be loaded, otherwise funcs will not work.
import numpy as np
import ROOT as rt
from util.vectorcalcs import BuildVectorCalcs, LoadVectorCalcs

# Always build/load library when loading this.
# Putting the commands outside of any func makes sure this is loaded.
BuildVectorCalcs()
LoadVectorCalcs()

def DeltaR2(eta1,phi1,eta2,phi2):
    return rt.VectorCalcs.DeltaR2(eta1,phi1,eta2,phi2)

# vec1 and vec2 are lists of (eta,phi) coordinates.
# If len(vec1) = n and len(vec2) = m, returns an array
# of shape n x m. Component [i][j] is the distance between
# vec1[i] and vec2[j].
def DeltaR2Vectorized(vec1, vec2):
    # TODO: Running into some runtime errors with DeltaR2Vectorized, preceded by "cling JIT session error: Cannot allocate memory". Hard to reproduce outside the routine.
    n,m = vec1.shape[0], vec2.shape[0]
    distances = np.array(rt.VectorCalcs.DeltaR2Vectorized(vec1[:,0].flatten(),vec1[:,1].flatten(),vec2[:,0].flatten(),vec2[:,1].flatten())) # testing shows this is faster than using DeltaR2 w/ list comprehension
    return np.reshape(distances,(n,m))

# Make sure that phi is in [-pi, pi].
# TODO: This is kind of ugly, is there a method that works for both arrays and scalars?
def AdjustPhi(phi):
    phi = np.asarray(phi)
    phi[phi > np.pi] -= 2.0 * np.pi
    phi[phi < -np.pi] += 2.0 * np.pi
    if(phi.shape == (1,)): return phi[0]
    return phi

def PxPyPzEToPtEtaPhiM(px,py,pz,e,transpose=True):
    nvecs = len(px)
    input_vecs = np.vstack((px,py,pz,e)).T
    output_vecs = np.array(rt.VectorCalcs.PxPyPzE2PtEtaPhiMflat(input_vecs.flatten())).reshape((nvecs,-1))
    if(transpose): return output_vecs
    return output_vecs.T

def PtEtaPhiMToPxPyPzE(pt,eta,phi,m,transpose=True):
    return PtEtaPhiMToPxPyPzE_root(pt,eta,phi,m,transpose)

def PtEtaPhiMToEPxPyPz(pt,eta,phi,m,transpose=True):
    return PtEtaPhiMToEPxPyPz_root(pt,eta,phi,m,transpose)

def EPxPyPzToM(e,px,py,pz):
    return np.sqrt(np.square(e) - np.square(px) -np.square(py) - np.square(pz))

def EPxPyPzToPtEtaPhiM(e,px,py,pz,transpose=True):
    vec = rt.Math.PxPyPzEVector()
    vec.SetCoordinates(px,py,pz,e)
    pt = vec.Pt()
    eta = vec.Eta()
    phi = vec.Phi()
    m = vec.M()
    v = np.array([pt,eta,phi,m],dtype=np.dtype('f8'))
    if(transpose): return v.T
    return v

def EPzToRap(e,pz):
    return 0.5 * np.log((e + pz)/(e - pz))

#--------------------------------------------------
# Below are some functions that may get used above -- currently experimenting with purely numpy-based functions,
# versus using my custom ROOT/C++ library. TODO: When possible we want to avoid conversion due to numerical precision
# issues. The ROOT-based method (e.g. using my "VectorCalcs" library) sometimes seem more precise, though also slower.
#--------------------------------------------------

def PtEtaPhiMToPxPyPzE_numpy(pt,eta,phi,m,transpose=True):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e  = np.sqrt(px * px + py * py + pz * pz + m * m)
    if(transpose): return np.array([px,py,pz,e],dtype=np.dtype('f8')).T
    else: return np.array([px,py,pz,e],dtype=np.dtype('f8'))

def PtEtaPhiMToPxPyPzE_root(pt,eta,phi,m,transpose=True):
    input_vecs = np.vstack((pt,eta,phi,m)).T
    rvecs = [rt.Math.PtEtaPhiMVector(*v) for v in input_vecs]
    output_vecs = np.array([[v.Px(),v.Py(),v.Pz(),v.E()] for v in rvecs],dtype=np.dtype('f8'))
    if(not transpose): output_vecs = output_vecs.T
    return output_vecs

def PtEtaPhiMToEPxPyPz_numpy(pt,eta,phi,m,transpose=True):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e  = np.sqrt(px * px + py * py + pz * pz + m * m)
    if(transpose): return np.array([e,px,py,pz],dtype=np.dtype('f8')).T
    else: return np.array([e,px,py,pz],dtype=np.dtype('f8'))

def PtEtaPhiMToEPxPyPz_root(pt,eta,phi,m,transpose=True):
    nvecs = len(pt)
    input_vecs = np.vstack((pt,eta,phi,m)).T
    output_vecs = np.array(rt.VectorCalcs.PtEtaPhiM2EPxPyPzflat(input_vecs.flatten())).reshape((nvecs,-1))
    if(transpose): return output_vecs
    return output_vecs.T