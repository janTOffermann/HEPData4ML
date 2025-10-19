# Some functions for filtering out events during generation -- if an event fails to pass some condition
# it can be thrown out instead of being written to the HepMC file -- and it is not counted towards the
# total number of events we have generated.
import ROOT as rt
from util.pythia.setup import PythiaWrapperSetup

class BaseFilter:
    def __init__(self):

        self.setup = PythiaWrapperSetup()
        self.setup.FullPreparation()

        self.filter = None

    def GetFilter(self):
        return self.filter # actual access to the underlying ROOT.NtupleProducer event filter


class ParticleFilter(BaseFilter):
    def __init__(self, statusHepMC:int,pdgid:int,ptMin:float=-1.,ptMax:float=-1.,etaMin=-999.,etaMax=999.):
        super().__init__()
        self.filter = rt.PythiaGenerator.ParticleFilter(pdgid, statusHepMC, ptMin, ptMax, etaMin, etaMax)
