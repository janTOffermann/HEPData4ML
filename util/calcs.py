# Some 4-vector code based on the ROOT (root.cern.ch) source code.
import numpy as np
from numba import jit

@jit
def DeltaPhi(phi1, phi2):
    dphi = phi2 - phi1
    if(dphi > np.pi): dphi -= 2.0 * np.pi
    elif(dphi < -np.pi): dphi += 2.0 * np.pi
    return dphi

@jit
def DeltaR2(eta1,phi1,eta2,phi2):
    dphi = DeltaPhi(phi1,phi2)
    deta = eta2 - eta1
    return dphi*dphi + deta*deta

# Calculate the arctangent of a/b, account for cases of b=0.
@jit
def SafeArctan(a,b):
    ans = np.full(a.shape, np.pi / 2.) # pi/2 corresponds to b=0, a>0
    neg_zero = (b==0.) * (a < 0.) # cases of b==0 w/ a < 0
    ans[neg_zero] = -np.pi / 2

    # Now deal with cases of b!=0.
    nzero = b!=0 # cases of b!=0
    ans[nzero] = np.arctan(a[nzero]/b[nzero])

    # Now adjust answers for quadrant II (b < 0, a > 0).
    ans[(b < 0) * (a > 0)] += np.pi
    ans[(b < 0) * (a <= 0)] -= np.pi
    return ans

# vec1 and vec2 are lists of (eta,phi) coordinates.
# If len(vec1) = n and len(vec2) = m, returns an array
# of shape n x m. Component [i][j] is the distance between
# vec1[i] and vec2[j].
@jit
def DeltaR2Vectorized(vec1, vec2):
    n,m = vec1.shape[0], vec2.shape[0]
    distances = np.full((n,m),-1.,dtype=np.dtype('f8'))
    for i in range(n):
        for j in range(m):
            distances[i,j] = DeltaR2(vec1[i][0],vec1[i][1],vec2[j][0],vec2[j][1])
    return distances

# Make sure that phi is in [-pi, pi].
# TODO: This is kind of ugly, is there a method that works for both arrays and scalars?
def AdjustPhi(phi):

    phi = np.asarray(phi)
    phi[phi > np.pi] -= 2.0 * np.pi
    phi[phi < -np.pi] += 2.0 * np.pi

    if(phi.shape == (1,)): return phi[0]
    return phi

def PtEtaPhiMToPxPyPzE(pt,eta,phi,m, transpose=True):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e  = np.sqrt(px * px + py * py + pz * pz + m * m)
    if(transpose): return np.array([px,py,pz,e],dtype=np.dtype('f8')).T
    else: return np.array([px,py,pz,e],dtype=np.dtype('f8'))

def PtEtaPhiMToEPxPyPz(pt,eta,phi,m, transpose=True):
    vecs = PtEtaPhiMToPxPyPzE(pt,eta,phi,m, transpose=True) # gives (px, py, pz, e)
    energies = vecs[:,3]
    vecs[:,1:3] = vecs[:,0:2]
    vecs[:,0] = energies
    if(transpose): return vecs # because vecs had transposition applied, already transpose -> don't do .T
    else: return vecs.T # vecs was transpose -> to undo transposition, do .T again

def PxPyPzEToPtEtaPhiM(px,py,pz,e, transpose=True):
    pt = np.sqrt(px * px + py * py)
    p2 = (px * px) + (py * py) + (pz * pz)

    eta = np.arctanh(pz / np.sqrt(p2))
    phi = SafeArctan(py,px)
    m = np.sqrt((e * e - p2).clip(min=0.))

    # Correct phi values, so that they are in [-pi,pi].
    phi = AdjustPhi(phi)

    if(transpose): return np.array([pt,eta,phi,m],dtype=np.dtype('f8')).T
    else: return np.array([pt,eta,phi,m],dtype=np.dtype('f8'))

