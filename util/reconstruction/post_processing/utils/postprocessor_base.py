
from typing import TYPE_CHECKING, Optional
import datetime

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.reconstruction.post_processing.jets import JetFinder

class PostProcessorBase:
    """
    Base class for post-processors.
    """

    def __init__(self):

        self.name = 'PostProcessor'
        self.print_prefix = '{}'.format(self.name)
        self.citations = {}
        self.obj_name_input = None
        self.obj_name_output = None

        # Allow separate handling of input name for jet constituents;
        # this is because some post-processors write new jet-level outputs,
        # but *don't* produce their own constituent-level outputs.
        # Any subsequent post-processor will need to be pointed to the correct
        # constituent-level branch if constituent-level info is needed.
        self.obj_name_constituents_input = None
        self.obj_name_constituents_output = None

        # for timing purposes
        self.runtime = 0.

    def AddRuntime(self,value:float):
        self.runtime += value

    def GetCitations(self):
        return self.citations

    def SetInputObjectName(self,obj_name:str, obj_constituents_name:Optional[str]=None):
        self.obj_name_input = obj_name
        self.obj_name_constituents_input = obj_constituents_name
        if(self.obj_name_constituents_input is None):
            self.obj_name_constituents_input = self.obj_name_input

    def GetInputObjectName(self):
        return self.obj_name_input

    def GetOutputObjectName(self):
        return self.obj_name_output

    def GetOutputObjectConstituentsName(self):
        if(self.obj_name_constituents_output is None):
            return self.GetOutputObjectName()
        return self.obj_name_constituents_output

    def SummarizeIO(self):
        if(self.obj_name_input == self.obj_name_output):
            return # nothing to print
        self._print('     Input: {}'.format(self.obj_name_input))
        if(self.obj_name_input != self.obj_name_constituents_input):
            self._print('     \tConstituents: {}'.format(self.obj_name_constituents_input))
        self._print('    Output: {}'.format(self.obj_name_output))

    def SummarizeRuntime(self,level=0):
        runtime_readable = str(datetime.timedelta(seconds=self.runtime))
        self._print('     Runtime: {:.1f} seconds\t({})'.format(self.runtime,runtime_readable),level=level)
        return

    def _print(self,val,level=0):
        prefix = level * '\t' + self.print_prefix
        print('{}: {}'.format(prefix,val))
        return

    # ----------------------------------------

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def ModifyConstituents(self, obj : 'JetFinder'):
        return

    def ModifyWrite(self,obj:'JetFinder'):
        return # does nothing