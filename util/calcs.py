# Note: Functions here use our custom VectorCalcs C++/ROOT library.
#       That library must be loaded, otherwise funcs will not work.
import numpy as np
import ROOT as rt
from util.vectorcalcs import VectorCalcsManager

class Calculator:
    def __init__(self):
        self.vc_manager = VectorCalcsManager()
        self.vc_manager.FullPreparation(force=False)
        self.calculator = rt.VectorCalcs.Calculator()

    def __del__(self): # TODO: Not sure if this is needed?
        del self.calculator

    def DeltaR2(self,eta1,phi1,eta2,phi2):
        return self.calculator.DeltaR2(eta1,phi1,eta2,phi2)

    def DeltaR2Vectorized(self,vec1,vec2):
        n,m = vec1.shape[0], vec2.shape[0]
        return np.reshape( self.calculator.DeltaR2Vectorized( vec1[:,0].flatten(),vec1[:,1].flatten(),vec2[:,0].flatten(),vec2[:,1].flatten()) , (n,m))

    def AdjustPhi(self,phi):
        phi = np.asarray(phi)
        phi[phi > np.pi] -= 2.0 * np.pi
        phi[phi < -np.pi] += 2.0 * np.pi
        if(phi.shape == (1,)): return phi[0]
        return phi

    def EPxPyPzToPtEtaPhiM(self,e,px,py,pz,transpose=True): # TODO: Some inconsistent naming/design here, some funcs act on single vectors, others act on many!
        """
        Acts on a single vector.
        """
        vec = rt.Math.PxPyPzEVector()
        vec.SetCoordinates(px,py,pz,e)
        pt = vec.Pt()
        eta = vec.Eta()
        phi = vec.Phi()
        m = vec.M()
        v = np.array([pt,eta,phi,m],dtype=np.dtype('f8'))
        if(transpose): return v.T
        return v

    def EPzToRap(self,e,pz):
        return 0.5 * np.log((e + pz)/(e - pz))

    def PxPyPzEToPtEtaPhiM(self,px,py,pz,e,transpose=True):
        """
        Acts on a set of vectors (separate lists/arrays of px, py, pz and e).
        """
        nvecs = len(px)
        input_vecs = np.vstack((px,py,pz,e)).T
        output_vecs = np.array(self.calculator.PxPyPzE_to_PtEtaPhiM_Multi(input_vecs.flatten())).reshape((nvecs,-1))
        if(transpose): return output_vecs
        return output_vecs.T

    def PtEtaPhiMToEPxPyPz(self,pt,eta,phi,m,transpose=True):
        """
        Acts on a set of vectors (separate lists/arrays of pt, eta, phi and m).
        """
        nvecs = len(pt)
        input_vecs = np.vstack((pt,eta,phi,m)).T
        output_vecs = np.array(self.calculator.PtEtaPhiM_to_EPxPyPz_Multi(input_vecs.flatten())).reshape((nvecs,-1))
        if(transpose): return output_vecs
        return output_vecs.T