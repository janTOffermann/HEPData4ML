# Note: Functions here use our custom VectorCalcs C++/ROOT library.
#       That library must be loaded, otherwise funcs will not work.
import numpy as np
import ROOT as rt
from util.vectorcalcs import VectorCalcsManager

def embed_array(array, target_shape,padding_value=0):
    """
    A generic function for embedding an input array into some target shape.
    The array will be truncated or padded as needed.

    """
    # Create the output array
    result = np.full(target_shape, padding_value, dtype=array.dtype)

    # Calculate the effective shape (minimum dimensions)
    effective_shape = tuple(min(s, t) for s, t in zip(array.shape, target_shape))

    slices = tuple(slice(0, dim) for dim in effective_shape)
    result[slices] = array[slices]
    return result

def embed_array_inplace(array, target, padding_value=0):
    """
    A generic function for embedding an input array into some target array.
    The array will be truncated or padded as needed.
    Modifies target in-place.
    """
    array_np = np.array(array)

    # Calculate the effective shape (minimum dimensions)
    effective_shape = tuple(min(s, t) for s, t in zip(array_np.shape, target.shape))

    # Create slices for both arrays
    slices = tuple(slice(0, dim) for dim in effective_shape)

    # Reset target to padding value
    target.fill(padding_value)

    # Copy data from array to target
    target[slices] = array_np[slices]

    return

# Note: Some functions have a funny-looking structure - they
#       will try to use our custom "VectorCalcs" C++/ROOT library,
#       but fall back on pure Python if there is an issue. This is
#       due to some kind of bug I have encountered when using this
#       code on the UChicago Analysis Facility, which is a bug I
#       have not been able to reproduce elsewhere. It involves
#       an error with "cling JIT" not being able to allocate
#       memory, so this is a bit of a workaround. - Jan

class Calculator:
    def __init__(self, use_vectorcalcs=True):
        self.vc_manager = VectorCalcsManager()
        self.use_vectorcalcs = use_vectorcalcs
        if(self.use_vectorcalcs):
            self.vc_manager.FullPreparation(force=False)
            self.calculator = rt.VectorCalcs.Calculator()

        # for non-C++/ROOT usage.
        self.v1 = rt.Math.PtEtaPhiMVector(0.,0.,0.,0.)
        self.v2 = rt.Math.PtEtaPhiMVector(0.,0.,0.,0.)
        self.v3 = rt.Math.PxPyPzEVector(0.,0.,0.,0.)

    def __del__(self): # TODO: Not sure if this is needed?
        if(self.use_vectorcalcs):
            del self.calculator

    def DeltaR2(self,eta1,phi1,eta2,phi2):
        if(self.use_vectorcalcs):
            try:
                return self._DeltaR2_root(eta1,phi1,eta2,phi2)
            except:
                self.use_vectorcalcs = False
                return self.DeltaR2(eta1,phi1,eta2,phi2)
        else:
            return self._DeltaR2_numpy(eta1,phi1,eta2,phi2)

    def _DeltaR2_root(self,eta1,phi1,eta2,phi2):
        return self.calculator.DeltaR2(eta1,phi1,eta2,phi2)

    def _DeltaR2_numpy(self,eta1,phi1,eta2,phi2):
        deta = eta2 - eta1
        dphi = self._dPhi(phi2,phi1)
        result = np.square(deta) + np.square(dphi)
        return result

    def DeltaR2Vectorized(self,vec1,vec2):
        if(self.use_vectorcalcs):
            try:
                return self._DeltaR2Vectorized_root(vec1,vec2)
            except:
                self.use_vectorcalcs = False
                return self.DeltaR2Vectorized(vec1,vec2)
        else:
            return self._DeltaR2Vectorized_numpy(vec1,vec2)

    def _DeltaR2Vectorized_root(self,vec1,vec2):
        n,m = vec1.shape[0], vec2.shape[0]
        return np.reshape( self.calculator.DeltaR2Vectorized( vec1[:,0].flatten(),vec1[:,1].flatten(),vec2[:,0].flatten(),vec2[:,1].flatten()) , (n,m))

    def _DeltaR2Vectorized_numpy(self,vec1,vec2):
        return np.array([[self._DeltaR2_numpy(*x,*y) for y in vec2] for x in vec1])

    def EPxPyPzToPtEtaPhiM_single(self,e,px,py,pz,transpose=True):
        """
        Acts on a single vector.
        """
        if(self.use_vectorcalcs):
            try:
                return self._EPxPyPzToPtEtaPhiM_single_root(e,px,py,pz,transpose)
            except:
                self.use_vectorcalcs = False
                return self.EPxPyPzToPtEtaPhiM_single(e,px,py,pz,transpose)
        else:
            return self._EPxPyPzToPtEtaPhiM_single_numpy(e,px,py,pz,transpose)

    def _EPxPyPzToPtEtaPhiM_single_root(self,e,px,py,pz,transpose=True):
        vec = np.array([px,py,pz,e])
        return np.array(self.calculator.PxPyPzE_to_PtEtaPhiM(vec.flatten()))

    def _EPxPyPzToPtEtaPhiM_single_numpy(self,e,px,py,pz,transpose=True):
        pt = np.linalg.norm([px,py])
        p = np.linalg.norm([px,py,pz])
        if(p == 0.): eta = 0.
        else:
            eta = np.arctanh(pz/p)
        phi = np.arctan2(py,px)
        m2 = np.square(e) - np.square(p)
        m = np.sign(m2) * np.sqrt(np.abs(m2))
        v = np.array([pt,eta,phi,m],dtype=np.dtype('f8'))
        if(transpose): return v.T
        return v

    def PxPyPzEToPtEtaPhiM(self,px,py,pz,e,transpose=True):
        """
        Acts on a set of vectors (separate lists/arrays of px, py, pz and e).
        """
        if(self.use_vectorcalcs):
            try:
                return self._PxPyPzEToPtEtaPhiM_root(px,py,pz,e,transpose)
            except:
                self.use_vectorcalcs = False
                return self.PxPyPzEToPtEtaPhiM(px,py,pz,e,transpose)
        else:
            return self._PxPyPzEToPtEtaPhiM_numpy(px,py,pz,e,transpose)

    def _PxPyPzEToPtEtaPhiM_root(self,px,py,pz,e,transpose=True):
        nvecs = len(px)
        input_vecs = np.vstack((px,py,pz,e)).T
        output_vecs = np.array(self.calculator.PxPyPzE_to_PtEtaPhiM_Multi(input_vecs.flatten())).reshape((nvecs,-1))
        if(transpose): return output_vecs
        return output_vecs.T

    def _PxPyPzEToPtEtaPhiM_numpy(self,px,py,pz,e,transpose=True):
        vecs = np.vstack((e,px,py,pz)).T
        output_vecs = np.array([self._EPxPyPzToPtEtaPhiM_single_numpy(*x) for x in vecs])
        if(transpose): return output_vecs
        return output_vecs.T

    def PtEtaPhiMToEPxPyPz(self,pt,eta,phi,m,transpose=True):
        """
        Acts on a set of vectors (separate lists/arrays of pt, eta, phi and m).
        """
        nvecs = len(pt)
        input_vecs = np.vstack((pt,eta,phi,m)).T
        if(self.use_vectorcalcs):
            try:
                output_vecs = np.array(self.calculator.PtEtaPhiM_to_EPxPyPz_Multi(input_vecs.flatten())).reshape((nvecs,-1))
            except:
                self.use_vectorcalcs = False
                return self.PtEtaPhiMToEPxPyPz(pt,eta,phi,m,transpose)
        else:
            output_vecs = np.zeros((nvecs,4))
            for i,ivec in enumerate(input_vecs):
                self.v1.SetCoordinates(*ivec)
                output_vecs[i,:] = np.array([self.v3.E(), self.v3.Px(), self.v3.Py(), self.v3.Pz()])
        if(transpose): return output_vecs
        return output_vecs.T

    def _PtEtaPhiMToEPxPyPz_root(self,pt,eta,phi,m,transpose=True):
        nvecs = len(pt)
        input_vecs = np.vstack((pt,eta,phi,m)).T
        output_vecs = np.array(self.calculator.PtEtaPhiM_to_EPxPyPz_Multi(input_vecs.flatten())).reshape((nvecs,-1))
        if(transpose): return output_vecs
        return output_vecs.T

    def EPzToRap(self,e,pz):
        return 0.5 * np.log((e + pz)/(e - pz))

    def AdjustPhi(self,phi):
        """
        Adjusts phi values so that they lie in [-pi,pi].
        """
        phi = np.asarray(phi)
        phi[phi > np.pi] -= 2.0 * np.pi
        phi[phi < -np.pi] += 2.0 * np.pi
        if(phi.shape == (1,)): return phi[0]
        return phi

    def _dPhi(self,phi1,phi2):
        """
        Delta phi. Based off of the ROOT::Math::VectorUtil source code.
        """
        # see https://root.cern/doc/v608/GenVector_2VectorUtil_8h_source.html#l00061
        dphi = phi2 - phi1
        if(dphi > np.pi):
            dphi -= 2. * np.pi
        elif(dphi <= -np.pi):
            dphi += 2. * np.pi
        return dphi

    def EtaToTheta(self,eta):
        return 2. * np.arctan(np.exp(eta)) - 0.5 * np.pi