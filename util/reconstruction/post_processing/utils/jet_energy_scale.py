import ROOT as rt
import numpy as np
from typing import TYPE_CHECKING
from util.reconstruction.post_processing.utils.postprocessor_base import PostProcessorBase

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class JetEnergyScale(PostProcessorBase):
    """
    Performs jet energy scale calibration.
    This will apply the calibration to the jet
    and constituent Pmu/Pmu_cyl branches, and will
    also store the per-jet calibration factor so
    that one can in principle recover the uncalibrated
    jet/constituent momenta.
    """

    def __init__(self,formula:str=None,override=False):
        super().__init__()

        self.formula = formula
        self.tformula = None
        self.override = override # if True, override check on default ATLAS/CMS formulae

        self.buffers = {
            'CalibrationFactor':None,
            'Pmu':None,
            'Pmu_cyl':None,
            'Constituents.Pmu':None,
            'Constituents.Pmu_cyl':None
        }
        self.branch_names = {key:None for key in self.buffers.keys()}

        self.name = 'JetEnergyScale'
        self.print_prefix = '{}'.format(self.name)
        self.citations = {}

    def ModifyInitialization(self,obj : 'JetFinder'):
        if(self.obj_name_input is None):
            self.SetInputObjectName(obj.jet_name)
        self.obj_name_output = '.'.join([self.obj_name_input,'Calibrated'])
        self.obj_name_constituents_output = self.obj_name_output

        self._set_formula(self.formula, obj)
        self._initialize_func()
        self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will compute and apply the jet energy scale calibration.
        """
        self.buffers['CalibrationFactor'] = {}
        self.buffers['Pmu'] = {}
        self.buffers['Pmu_cyl'] = {}

        pmu_key = '{}.Pmu'.format(self.obj_name_input)
        pmu_cyl_key = '{}.Pmu_cyl'.format(self.obj_name_input)

        # Compute the actual jet energy scale factors and apply them to the jet.
        pmu_dict = obj.output_buffer.get(pmu_key,filter=obj.jet_ordering)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        for i,pmu_cyl in pmu_cyl_dict.items():
            pmu = pmu_dict[i]
            self.buffers['CalibrationFactor'][i] = self.tformula.Eval(
                float(pmu_cyl[0]),
                float(pmu_cyl[1]),
                float(pmu_cyl[2]),
                float(pmu[0])
            )

            self.buffers['Pmu'][i] = self.buffers['CalibrationFactor'][i] * pmu
            self.buffers['Pmu_cyl'][i] = pmu_cyl
            self.buffers['Pmu_cyl'][i][0] *= self.buffers['CalibrationFactor'][i] # pt
            self.buffers['Pmu_cyl'][i][-1] *= self.buffers['CalibrationFactor'][i] # mass
        return

    def ModifyConstituents(self, obj : 'JetFinder'):
        self.buffers['Constituents.Pmu'] = {}
        self.buffers['Constituents.Pmu_cyl'] = {}

        pmu_key = '{}.Constituents.Pmu'.format(self.obj_name_constituents_input)
        pmu_cyl_key = '{}.Constituents.Pmu_cyl'.format(self.obj_name_constituents_input)

        pmu_dict = obj.output_buffer.get(pmu_key,filter=obj.jet_ordering)
        pmu_cyl_dict = obj.output_buffer.get(pmu_cyl_key,filter=obj.jet_ordering)

        for i,pmu in pmu_dict.items():
            pmu_cyl = pmu_cyl_dict[i]
            self.buffers['Constituents.Pmu'][i] = self.buffers['CalibrationFactor'][i] * pmu
            self.buffers['Constituents.Pmu_cyl'][i] = pmu_cyl
            self.buffers['Constituents.Pmu_cyl'][i][0] *= self.buffers['CalibrationFactor'][i] # pt
            self.buffers['Constituents.Pmu_cyl'][i][-1] *= self.buffers['CalibrationFactor'][i] # mass
        return

    def ModifyWrite(self,obj : 'JetFinder'):
        self._addToBuffer(obj)
        return

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Adds a branch to the buffer corresponding to the calibration factor.
        """
        self._createBranchNames(obj)

        if(self.branch_names['CalibrationFactor'] not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_names['CalibrationFactor'],ndim=1,dtype=np.dtype('f8'))
        if(self.branch_names['Pmu'] not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_names['Pmu'],ndim=2,dtype=np.dtype('f8'))
        if(self.branch_names['Pmu_cyl'] not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_names['Pmu_cyl'],ndim=2,dtype=np.dtype('f8'))
        if(self.branch_names['Constituents.Pmu'] not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_names['Constituents.Pmu'],ndim=3,dtype=np.dtype('f8'))
        if(self.branch_names['Constituents.Pmu_cyl'] not in obj.output_buffer.keys()):
            obj.output_buffer.create_array(self.branch_names['Constituents.Pmu_cyl'],ndim=3,dtype=np.dtype('f8'))

    def _createBranchNames(self,obj : 'JetFinder'):
        for key in self.branch_names.keys():
            if('Constituents' in key):
                self.branch_names[key] = '{}.{}'.format(self.obj_name_constituents_output,key)
            else:
                self.branch_names[key] = '{}.{}'.format(self.obj_name_output,key)
        return

    def _addToBuffer(self,obj : 'JetFinder'):
        """
        Adds the calibrated momenta and calibration factor to the buffer.
        """
        for key in self.branch_names.keys():
            obj.output_buffer.set(self.branch_names[key],obj._i,self.buffers[key])

    def _set_formula(self,formula_str=None, obj:'JetFinder'=None):

        if(formula_str is None):
            self._print('Warning: Defaulting to ATLAS JES formula (in Delphes), for R=0.4 anti-kt jets.')
            formula_str='ATLAS' # default to ATLAS

        if(formula_str.lower()=='atlas'):
            if(obj.radius == 0.4 or self.override):
                # default formula from ATLAS detector card
                formula_str = 'sqrt( (3.0 - 0.2*(abs(eta)))^2 / pt + 1.0 )'
            else:
                self._print('Warning: Requested default ATLAS JES formula (from Delphes), but it is not available for R!=0.4 .')
                formula_str = '1'

        elif(formula_str.lower()=='cms'):
            if(obj.radius == 0.4 or self.override):
                # default formula from CMS detector card
                formula_str = 'sqrt( (2.5 - 0.15*(abs(eta)))^2 / pt + 1.0 )'
            else:
                self._print('Warning: Requested default CMS JES formula (from Delphes), but it is not available for R!=0.4 .')
                formula_str = '1'

        formula = formula_str.lower()
        formula = formula.replace('pt','x')
        formula = formula.replace('eta','y')
        formula = formula.replace('phi','z')
        formula = formula.replace('energy','t')
        self.formula = formula

    def _initialize_func(self):
        self.tformula = rt.TFormula('f_JES',self.formula)