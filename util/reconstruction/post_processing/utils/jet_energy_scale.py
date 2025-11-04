import ROOT as rt
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class JetEnergyScale:
    """
    Performs jet energy scale calibration.
    This will apply the calibration to the jet
    and constituent Pmu/Pmu_cyl branches, and will
    also store the per-jet calibration factor so
    that one can in principle recover the uncalibrated
    jet/constituent momenta.
    """

    def __init__(self,formula:str=None):

        self.formula = None
        self.tformula = None
        self._set_formula(formula)

        self.scale_factors = None

        self.name = 'JetEnergyScale'
        self.branch_name = None
        self.print_prefix = '\n\t{}'.format(self.name)
        self.citations = {}

    def GetCitations(self):
        return self.citations

    def ModifyInitialization(self,obj):
        self._initialize_func()
        return

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will compute and apply the jet energy scale calibration.
        """

        # Compute the actual jet energy scale factors and apply them to the jet.
        for i,jet in obj.jets_dict.items():

            self.scale_factors[i] = self.tformula.Eval(
                jet.pt(),
                jet.eta(),
                jet.phi(),
                jet.energy()
            )
            obj.jets_dict[i] *= self.scale_factors[i] # modify the jet itself!

        obj._jetsToVectors() # force update of vectors
        obj._ptSort() # re-sorts obj.jet_ordering

        return

    def ModifyWrite(self,obj : 'JetFinder'):
        self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
        self._addToBuffer(obj)
        return

    def ModifyConstituents(self, obj : 'JetFinder'):
        """
        Apply the jet energy scale calibration to the individual constituents.
        """
        # scale factors were computed by ModifyJets(), which is always called before.
        for i in obj.jets_dict.keys():
            obj.constituent_vectors[i] *= self.scale_factors[i] # scale (px, py, pz, E)
            obj.constituent_vectors_cyl[i,:,0] *= self.scale_factors[i] # scale pt
            obj.constituent_vectors_cyl[i,:,-1] *= self.scale_factors[i] # scale m
        return

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Adds a branch to the buffer corresponding to the calibration factor.
        """
        self._createBranchNames(obj)

        if(self.branch_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.branch_name,ndim=1,dtype=np.dtype('f8'))

    def _createBranchNames(self,obj : 'JetFinder'):
        self.branch_name = '{}.{}.CalibrationFactor'.format(obj.jet_name,self.name)
        return

    def _addToBuffer(self,obj : 'JetFinder'):
        """
        Adds the calibration factor to the buffer.
        """
        #NOTE: The embed is not needed, since we've constructed the inputs and the buffer to already match in size.
        #      The zero-padding is actually being handled within self.ModifyJets(), where the embed function is used.

        obj.buffer.set(self.branch_name,obj._i,np.vstack([self.scale_factors[i] for i in obj.jet_ordering]))

    def _set_formula(self,formula_str):

        if(formula_str is None):
            self._print('Warning: Defaulting to ATLAS JES formula.')
            formula_str='ATLAS' # default to ATLAS

        if(formula_str.lower()=='atlas'):
            # default formula from ATLAS detector card
            formula_str = 'sqrt( (3.0 - 0.2*(abs(eta)))^2 / pt + 1.0 )'

        elif(formula_str.lower()=='cms'):
            # default formula from CMS detector card
            formula_str = 'sqrt( (2.5 - 0.15*(abs(eta)))^2 / pt + 1.0 )'

        formula = formula_str.lower()
        formula = formula.replace('pt','x')
        formula = formula.replace('eta','y')
        formula = formula.replace('phi','z')
        formula = formula.replace('energy','t')
        self.formula = formula

    def _initialize_func(self):
        self.tformula = rt.TFormula('f_JES',self.formula)

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return