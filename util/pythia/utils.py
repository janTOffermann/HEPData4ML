import numpy as np
import pythia8 as pyth8

class PythiaWrapper:

    """
    A simple wrapper for Pythia8, performing much of the same functionality
    as the `numpythia` package but relying on an external Pythia8 installation.
    Also adds new functionality!
    """

    def __init__(self,verbose=False):
        self.pythia = pyth8.Pythia('',False)
        self.config_dict = {}

        self.event = None
        self.SetVerbose(verbose)
        self.initialized = False

    def GetPythia(self):
        return self.pythia

    def GetEvent(self):
        return self.pythia.event

    def SetVerbose(self,flag):
        self.verbose = flag

    def ReadString(self, string):
        self.pythia.readString(string)

    def ReadStrings(self,strings):
        for string in strings:
            string = string.strip('\n')
            string = string.split('#')[0]
            string = string.split('!')[0]
            if(string in ['','\n']): continue
            if(self.verbose): print('Reading string: {}'.format(string))
            self.ReadString(string)

    def ReadStringsFromFile(self,file):
        with open(file,'r') as f:
            strings = f.readlines()
            self.ReadStrings(strings)

    def ClearConfigDict(self):
        self.config_dict = {}

    def AddToConfigDict(self,cdict):
        for key,val in cdict.items():
            if(key in self.config_dict.keys() and self.verbose):
                if(self.config_dict[key] != val):
                    print('Warning: Overwriting configuration: {} .'.format(key))
                    print('\t[{}] -> [{}]'.format(self.config_dict[key],val))
            self.config_dict[key] = val

    def ReadConfigDict(self):
        strings = ['{} = {}'.format(key,val) for key,val in self.config_dict.items()]
        self.ReadStrings(strings)

    def PrintConfigDict(self):
        for key,val in self.config_dict.items():
            print('{} : {}'.format(key,val))

    # ---- A few convenient handles for configuration, mostly useful for testing ----
    def SetMPI(self, flag):
        self.AddToConfigDict({'PartonLevel:MPI':self._bool2string(flag)})

    def SetISR(self, flag):
        self.AddToConfigDict({'PartonLevel:ISR':self._bool2string(flag)})

    def SetFSR(self, flag):
        self.AddToConfigDict({'PartonLevel:FSR':self._bool2string(flag)})

    def SetPtHatMin(self, pt):
        self.AddToConfigDict({'PhaseSpace:pTHatMin':str(pt)})

    def SetPtHatMax(self, pt):
        self.AddToConfigDict({'PhaseSpace:pTHatMax':str(pt)})

    def SetPtHat(self, pt_min, pt_max):
        self.SetPtHatMin(pt_min)
        self.SetPtHatMax(pt_max)

    def SetQuiet(self,flag):
        self.AddToConfigDict({'Print:quiet':self._bool2string(flag)})

    def InitializePythia(self):
        self.ReadConfigDict()
        self.pythia.init()
        self.initialized = True

    # Generates an event, places it in self.event .
    def Generate(self):
        if(not self.initialized): self.InitializePythia()
        self.pythia.next()
        self.event = self.pythia.event

    # =============== Getters =============== #

    # Get the Pythia event object. Probably not super useful,
    # it will likely be better/easier to use other interfaces in this class.
    def GetEvent(self):
        return self.event

    # Get total number of particles in the event listing.
    def GetN(self):
        if(self.event is None): return -1
        return self.event.size() - 1 # drop entry 0, which represents "full event"
        # return len(self.GetPdgId())

    # Get PDG IDs.
    def GetPdgId(self, indices=None):
        if(self.event is None): return None
        pids = np.array([p.id() for p in self.event],dtype='i4')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): pids = pids[indices]
        return pids

    def GetStatus(self, indices=None, hepmc=False):
        if(self.event is None): return None
        if(not hepmc): status = np.array([p.status() for p in self.event],dtype='i4')[1:] # drop entry 0, which represents "full event"
        else: status = np.array([p.statusHepMC() for p in self.event],dtype='i4')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): status = status[indices]
        return status

    def _getMomentum(self, indices=None, format='e px py pz'):
        if(self.event is None): return None
        if(format == 'e px py pz'):
            pmu = np.array([[p.e(),p.px(),p.py(),p.pz()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        elif(format == 'px py pz e'):
            pmu = np.array([[p.px(),p.py(),p.pz(),p.e()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        elif(format == 'pt eta phi m'):
            pmu = np.array([[p.pT(),p.eta(),p.phi(),p.m()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        elif(format == 'pt eta phi e'):
            pmu = np.array([[p.pT(),p.eta(),p.phi(),p.e()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        elif(format == 'pt y phi e'):
            pmu = np.array([[p.pT(),p.y(),p.phi(),p.e()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        else: #'pt y phi m'
            pmu = np.array([[p.pT(),p.y(),p.phi(),p.m()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"

        if(indices is not None): pmu = pmu[indices]
        return pmu

    def GetEPxPyPz(self,indices=None):
        return self._getMomentum(indices,'e px py pz')

    def GetPxPyPzE(self,indices=None):
        return self._getMomentum(indices,'px py pz e')

    def GetPtEtaPhiM(self,indices=None):
        return self._getMomentum(indices,'pt eta phi m')

    def GetPtEtaPhiE(self,indices=None):
        return self._getMomentum(indices,'pt eta phi e')

    def GetPtYPhiM(self,indices=None):
        return self._getMomentum(indices,'pt y phi m')

    def GetPtYPhiE(self,indices=None):
        return self._getMomentum(indices,'pt y phi e')

    def GetEta(self,indices=None):
        if(self.event is None): return None
        eta = np.array([p.eta() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): eta = eta[indices]
        return eta

    def GetY(self,indices=None):
        if(self.event is None): return None
        y = np.array([p.y() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): y = y[indices]
        return y

    def GetTheta(self,indices=None):
        if(self.event is None): return None
        theta = np.array([p.theta() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): theta = theta[indices]
        return theta

    def GetPhi(self,indices=None):
        if(self.event is None): return None
        phi = np.array([p.phi() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): phi = phi[indices]
        return phi

    def GetM(self,indices=None):
        if(self.event is None): return None
        m = np.array([p.m() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): m = m[indices]
        return m

    def GetPt(self,indices=None):
        if(self.event is None): return None
        pt = np.array([p.pT() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): pt = pt[indices]
        return pt

    def GetE(self,indices=None):
        if(self.event is None): return None
        e = np.array([p.e() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): e = e[indices]
        return e

    def GetEt(self,indices=None):
        if(self.event is None): return None
        et = np.array([p.eT() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): et = et[indices]
        return et

    def GetXProd(self,indices=None):
        if(self.event is None): return None
        x = np.array([p.xProd() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): x = x[indices]
        return x

    def GetYProd(self,indices=None):
        if(self.event is None): return None
        x = np.array([p.yProd() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): x = x[indices]
        return x

    def GetZProd(self,indices=None):
        if(self.event is None): return None
        x = np.array([p.zProd() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): x = x[indices]
        return x

    def GetTProd(self,indices=None):
        if(self.event is None): return None
        x = np.array([p.tProd() for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): x = x[indices]
        return x

    def GetProd(self,indices=None):
        if(self.event is None): return None
        x = np.array([[p.tProd(),p.xProd(),p.yProd(),p.zProd()] for p in self.event],dtype='f8')[1:] # drop entry 0, which represents "full event"
        if(indices is not None): x = x[indices]
        return x

    def GetNames(self,indices=None):
        if(self.event is None): return None
        names = [p.name() for p in self.event][1:] # drop entry 0, which represents "full event"
        if(indices is not None): names = names[indices]
        return names

    def GetDaughter1(self,index):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        return np.array(self.event[idx].daughter1(),dtype='i4') - 1 # subtract 1 again to deal with Pythia's entry 0

    def GetDaughter2(self,index):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        return np.array(self.event[idx].daughter2(),dtype='i4') - 1 # subtract 1 again to deal with Pythia's entry 0

    def GetMother1(self,index):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        return np.array(self.event[idx].mother1(),dtype='i4') - 1 # subtract 1 again to deal with Pythia's entry 0

    def GetMother2(self,index):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        return np.array(self.event[idx].mother2(),dtype='i4') - 1 # subtract 1 again to deal with Pythia's entry 0

    def GetDaughters1(self):
        n = self.GetN()
        return [self.GetDaughter1(i) for i in range(n)]

    def GetDaughters2(self):
        n = self.GetN()
        return [self.GetDaughter2(i) for i in range(n)]

    def GetMothers1(self):
        n = self.GetN()
        return [self.GetMother1(i) for i in range(n)]

    def GetMothers2(self):
        n = self.GetN()
        return [self.GetMother2(i) for i in range(n)]

    # Return indices of particle's daughters.
    def GetDaughtersSingle(self,index, recursive=False):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        if(recursive): d = np.array(self.event[idx].daughterListRecursive(),dtype='i4')
        else: d = np.array(self.event[idx].daughterList(),dtype='i4')
        return d - 1 # subtract 1 again to deal with Pythia's entry 0

    # Return type is list, not array! This is because it is jagged -- could also consider awkward array.
    def GetDaughters(self, recursive=False):
        n = self.GetN()
        return [self.GetDaughtersSingle(i,recursive) for i in range(n)]

    # Return indices of a particle's mothers.
    def GetMothersSingle(self,index, recursive=False):
        idx = index + 1 # using index + 1 to effectively drop entry 0, which represents "full event"
        if(recursive): m = np.array(self.event[idx].motherListRecursive(),dtype='i4')
        else: m = np.array(self.event[idx].motherList(),dtype='i4')
        return m - 1 # subtract 1 again to deal with Pythia's entry 0

    def GetMothers(self, recursive=False):
        n = self.GetN()
        return [self.GetMothersSingle(i,recursive) for i in range(n)]

    # Get indices of a particle's stable daughters (hepmc status == 1).
    def GetStableDaughtersSingle(self, index, recursive=False):
        daughters = self.GetDaughtersSingle(index,recursive)
        status = self.GetStatus(hepmc=True)
        stable = np.where(status == 1)[0]
        return np.intersect1d(daughters,stable)

    # Get a single particle.
    # Format is (E, px, py, pz, pid, status)
    def GetParticle(self, index):
        status = self.GetStatus(index,hepmc=True)
        pid = self.GetPdgId(index)
        pmu = self._getMomentum(index,'e px py pz')
        return (*pmu, pid, status)
        # return np.array(
        #     [(*pmu,pid,status)],
        #     dtype=[('e','f8'),('px','f8'),('py','f8'),('pz','f8'),('pdgid','i4'),('status','i4')]
        # )

    # ================
    # Event-level info
    # ================

    # Get the code of the latest event process.
    def GetProcessCode(self):
        if(self.pythia is None): return None
        return self.pythia.infoPython().code()

    # Get the name of the latest event process.
    def GetProcessName(self):
        if(self.pythia is None): return None
        return self.pythia.infoPython().name()

    # Get the weight of the latest event.
    def GetEventWeight(self):
        if(self.pythia is None): return None
        return self.pythia.infoPython().weight()

    def GetWeightSum(self):
        return self.pythia.infoPython().weightSum()

    # Get the codes of all event processes that have been run.
    # Does not return a list of process code per event!
    def GetProcessCodes(self):
        return np.array(self.pythia.infoPython().codesHard(),dtype=int)

    # Get the estimated cross-section for a particular process (via its code).
    # Passing a code of 0 will give the cross-section for the sum of all active processes.
    def GetSigmaGen(self, i = 0):
        return self.pythia.infoPython().sigmaGen(i)

    # Get the uncertainty in a cross-section estimate.
    def GetSigmaErr(self, i = 0):
        return self.pythia.infoPython().sigmaErr(i)

    # Get a dictionary containing cross-sections (and their uncertainties) for all active processes,
    # arranged by the process codes.
    def GetSigmaDictionary(self):
        sigma_dict = {}
        codes = self.GetProcessCodes()
        for code in codes:
            sigma_dict[code] = (self.GetSigmaGen(code),self.GetSigmaErr(code))
        return sigma_dict

    def _bool2string(self,flag):
        if(flag): return 'on'
        return 'off'