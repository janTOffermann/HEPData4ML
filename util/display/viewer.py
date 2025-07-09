import ROOT as rt
import numpy as np
from util.display.setup import DisplaySetup


class EventDisplay:

    def __init__(self):

        # Make sure that our underlying EventDisplay
        # library is built and ready.
        self.setup = DisplaySetup()
        self.setup.Prepare()

    def InitializeDisplay(self,delphes_card):

        self.display = rt.EventDisplay.Display()

        self.display.DisplayEvent(delphes_card)