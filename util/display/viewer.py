import ROOT as rt
import numpy as np
import os
from util.display.setup import DisplaySetup


class EventDisplay:

    def __init__(self):

        # Make sure that our underlying EventDisplay
        # library is built and ready.
        self.setup = DisplaySetup()
        self.setup.Prepare()

    def InitializeDisplay(self,delphes_card):

        self.display = rt.EventDisplay.DisplayInterface()
        self.display.DisplayEvent(delphes_card)

        # this_dir = os.path.dirname(os.path.abspath(__file__))
        # filename = this_dir + '/test_data/atlas.root'
        # self.display.Test(filename)