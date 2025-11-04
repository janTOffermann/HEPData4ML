import ROOT as rt
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class JetEnergyScale:
    """
    Performs jet energy scale calibration.

    """

    def __init__(self,formula:str=None):

        self.formula = None
        self.tformula = None
        self._set_formula(formula)

        self.name = 'JetEnergyScale'
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
        This function will apply the jet energy scale calibration.
        """

        for i,jet in obj.jets_dict.items():
            # apply the calibration
            scale_factor = self.tformula.Eval(jet.pt(),jet.eta(),jet.phi(),jet.energy())
            obj.jets_dict[i] *= scale_factor

        obj._jetsToVectors() # force update of vectors
        obj._ptSort() # re-sorts obj.jet_ordering

        return

    def ModifyWrite(self,obj : 'JetFinder'):
        return # does nothing

    def ModifyConstituents(self, obj : 'JetFinder'):
        """
        Apply the jet energy scale calibration to the individual constituents.
        """
        for i,jet in obj.jets_dict.items():
            scale_factor = self.tformula.Eval(jet.pt(),jet.eta(),jet.phi(),jet.energy())
            obj.constituent_vectors[i] *= scale_factor # scale (px, py, pz, E)
            obj.constituent_vectors_cyl[i,:,0] *= scale_factor # scale pt
            obj.constituent_vectors_cyl[i,:,-1] *= scale_factor # scale m
        return

    def _set_formula(self,formula_str):

        if(formula_str is None):
            # default formula from ATLAS detector card
            formula_str = 'sqrt( (3.0 - 0.2*(abs(eta)))^2 / pt + 1.0 )'

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