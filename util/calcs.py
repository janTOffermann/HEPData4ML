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
    return rt.VectorCalcs.DeltaR2(eta1,phi2,eta2,phi2)

# vec1 and vec2 are lists of (eta,phi) coordinates.
# If len(vec1) = n and len(vec2) = m, returns an array
# of shape n x m. Component [i][j] is the distance between
# vec1[i] and vec2[j].
def DeltaR2Vectorized(vec1, vec2):
    n,m = vec1.shape[0], vec2.shape[0]
    distances = np.array(rt.VectorCalcs.DeltaR2Vectorized(vec1[:,0].flatten(),vec1[:,1].flatten(),vec2[:,0].flatten(),vec2[:,1].flatten()))
    distances = np.reshape(distances,(n,m))
    return distances

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
    nvecs = len(pt)
    input_vecs = np.vstack((pt,eta,phi,m)).T
    output_vecs = np.array(rt.VectorCalcs.PtEtaPhiM2PxPyPzEflat(input_vecs.flatten())).reshape((nvecs,-1))
    if(transpose): return output_vecs
    return output_vecs.T

def PtEtaPhiMToEPxPyPz(pt,eta,phi,m,transpose=True):
    nvecs = len(pt)
    input_vecs = np.vstack((pt,eta,phi,m)).T
    output_vecs = np.array(rt.VectorCalcs.PtEtaPhiM2EPxPyPzflat(input_vecs.flatten())).reshape((nvecs,-1))
    if(transpose): return output_vecs
    return output_vecs.T

# def PxPyPzEToPtEtaPhiM(px,py,pz,e,transpose=True):
#     n = len(px)
#     vec = rt.Math.PxPyPzEVector()
#     result = np.zeros((n,4),dtype=np.dtype('f8'))
#     for i in range(n):
#         vec.SetCoordinates(px[i],py[i],pz[i],e[i])
#         result[i,:] = np.array([vec.Pt(),vec.Eta(),vec.Phi(),vec.M()])
#     if(transpose): return result
#     return result.T

# def PtEtaPhiMToPxPyPzE(pt,eta,phi,m,transpose=True):
#     n = len(pt)
#     vec = rt.Math.PtEtaPhiMVector()
#     result = np.zeros((n,4),dtype=np.dtype('f8'))
#     for i in range(n):
#         vec.SetCoordinates(pt[i],eta[i],phi[i],m[i])
#         result[i,:] = np.array([vec.Px(),vec.Py(),vec.Pz(),vec.E()])
#     if(transpose): return result
#     return result.T

# def PtEtaPhiMToEPxPyPz(pt,eta,phi,m,transpose=True):
#     n = len(pt)
#     vec = rt.Math.PtEtaPhiMVector()
#     result = np.zeros((n,4),dtype=np.dtype('f8'))
#     for i in range(n):
#         vec.SetCoordinates(pt[i],eta[i],phi[i],m[i])
#         result[i,:] = np.array([vec.E(),vec.Px(),vec.Py(),vec.Pz()])
#     if(transpose): return result
#     return result.T