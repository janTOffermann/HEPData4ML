import numpy as np
import pythia8 as pyth8
import ROOT as rt
import awkward as ak
from util.pythia.setup import PythiaWrapperSetup
from util.misc.timing import profile_method, profile_block
from typing import Any

class BasicWrapper:
    """
    Abstract parent class for the Pythia wrappers.
    """

    def __init__(self):
        self.config_dict = {}
        self.initialized = False
        self.verbose = False

    def IsInitialized(self):
        return self.initialized

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

    def _bool2string(self,flag):
        if(flag): return 'on'
        return 'off'


class PythiaPythonWrapper(BasicWrapper):

    """
    A simple wrapper for Pythia8, performing much of the same functionality
    as the `numpythia` package but relying on an external Pythia8 installation.
    Also adds some new functionality!
    """
    def __init__(self,verbose=False):
        super().__init__()
        self.pythia = pyth8.Pythia('',False)
        self.event = None
        self.events = None # for holding batches of events, as awkward arrays
        self.SetVerbose(verbose)

    def GetPythia(self):
        return self.pythia

    # Generates an event, places it in self.event .
    def Generate(self):
        if(not self.initialized): self.InitializePythia()
        self.pythia.next()
        self.event = self.pythia.event

    # Generates batches of events -- faster than doing one by one!
    # Leverages awkward arrays.
    # Places results in self.event
    def GenerateBatch(self,batch_size):
        if(not self.initialized): self.InitializePythia()
        self.events = self.pythia.nextBatch(batch_size, errorMode='skip')
        return

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

class PythiaWrapper(BasicWrapper):
    """
    An alternative wrapper for Pythia8, that utilizes
    a custom C++/ROOT library.
    The interface is relatively rudimentary, mostly
    just functions needed for converting the Pythia8
    event record to HepMC3 format.
    """

    def __init__(self,verbose=False, parallel=False):
        super().__init__()

        self.setup = PythiaWrapperSetup()

        # Make sure that the underlying C++/ROOT library is built and loaded,
        # if not it will build here.
        self.setup.FullPreparation()
        self.pythia = rt.PythiaGenerator.Generator(parallel)

        self.config_dict = {}

        self.use_hepmc_status = True

        self.SetVerbose(verbose)
        self.initialized = False
        self.n_events = None
        self.particles_per_event = None

        # By default, turn off array and HepMC3 modes
        self.array_mode = False
        self.hepmc_mode = False
        self.hepmc_filename = 'output.hepmc.root'
        self.print_prefix = 'PythiaWrapper'

    def UseHepMCStatus(self,val:bool):
        self.use_hepmc_status = val

    def Generate(self, n):
        """
        Generate events.
        """
        if(not(self.array_mode or self.hepmc_mode)):
            self._print('Warning: Neither array nor HepMC3 modes turned on. Skipping event generation.')
            return

        if(not self.initialized): self.InitializePythia()

        self.pythia.generate(int(n),True)
        self.n_events = n

    def SetArrayMode(self,val:bool):
        self.array_mode = val
        self.pythia.setArrayMode(self.array_mode)

    def SetHepMC3Mode(self,val:bool):
        self.hepmc_mode = val
        self.pythia.setHepMC3Mode(self.hepmc_mode)

    def SetOutputFilename(self,val:str):
        self.SetHepMC3Mode(True)
        self.hepmc_filename = val

        if(self.hepmc_filename.split('.')[-1].lower() == 'root'):
            self.SetHepMC3WriteMode('root')
        else:
            self.SetHepMC3WriteMode('ascii')

    def SetHepMC3WriteMode(self,val:str):
        if(val.lower() == 'root'):
            self.pythia.setHepMC3RootWriter(True)
            self.pythia.setHepMC3AsciiWriter(False)
        else:
            self.pythia.setHepMC3RootWriter(False)
            self.pythia.setHepMC3AsciiWriter(True)

    def AddEventFilter(self,filter):
        self.pythia.addEventFilter(filter)

    def GetNEventsInBuffer(self):
        return self.pythia.GetNEventsInBuffer()

    # First, we have a ton of methods for accessing the data
    # produced in "array" mode, whereby the wrapper will fill
    # awkward arrays with information from the Pythia8 event record.
    # This is similar to the Pythia Python interface's "nextBatch()"
    # functionality, albeit not quite as fast.

    @profile_method('PythiaWrapper.GetData')
    def GetData(self,status_hepmc=None):
        """
        Provides output in awkward array format,
        similar to Pythia8's "nextBatch()" function.
        """

        # Our methods for fetching flat arrays will give pairs/tuples,
        # the last entry will be the number of particles per event.
        # But we don't want to keep converting this to a list, so we'll
        # actually cache it when we fetch particle ID.
        pid = self.GetPdgId()

        prt = ak.Array(
            {
                'id':pid,
                'status':self._getStatus(status_hepmc),
                'm':self.GetMass(),
                'p':self.GetMomentum(),
                'vProd':self.GetProdVertex(),
                'vProdStatus':self.GetProdVertexStatus(),
                # 'mother1':self.GetMother1(), # NOTE: motherList is much more useful
                # 'mother2':self.GetMother2(),
                'motherList':self.GetMotherList(),
                # 'daughter1':self.GetDaughter1(), # NOTE: daughterList is much more useful
                # 'daughter2':self.GetDaughter2(),
                'daughterList':self.GetDaughterList(),
                'col':self.GetColor(),
                'acol':self.GetAntiColor()
            }
        )

        event_info = ak.Array(
            {
                'id1':self.GetId1Pdf(),
                'id2':self.GetId2Pdf(),
                'pdf1':self.GetPdf1(),
                'pdf2':self.GetPdf2(),
                'x1':self.GetX1Pdf(),
                'x2':self.GetX2Pdf(),
                'QFac':self.GetQFac(),
                'QRen':self.GetQRen(),
                'nMPI':self.GetNumMPI(),
                'code':self.GetCode(),
                'alphaS':self.GetAlphaS(),
                'alphaEM':self.GetAlphaEM(),
                'sigmaGen':self.GetSigmaGen(),
                'sigmaErr':self.GetSigmaErr(),
                'weights':self.GetWeights()
            }
        )

        return ak.Array({'prt':prt,'info':event_info})

    def _getStatus(self,status_hepmc):
        if(status_hepmc is None):
            status_hepmc = self.use_hepmc_status
        if(status_hepmc):
            return self.GetStatusHepMC()
        return self.GetStatus()

    def _getNumParticlesPerEvent(self):
        if(self.particles_per_event is None):
            # still need to cache this
            _ = self.GetPdgId()
        return self.particles_per_event

    def GetPdgId(self):
        """
        Fetches particle ID.
        Also performs the special function of caching
        the number of particles per event.
        """

        pair = self.pythia.getParticleIDFlat()
        self.particles_per_event = list(pair[1])
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetStatus(self):
        pair = self.pythia.getStatusFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetStatusHepMC(self):
        pair = self.pythia.getStatusHepMCFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetMass(self):
        pair = self.pythia.getMassFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetMomentum(self):
        # TODO: Is it better to structure this differently?
        #       nextBatch() gives a jagged array of "records", each
        #       represents one momentum vector (its fields are of length 1).
        pair_e = self.pythia.getEnergyFlat()
        pair_px = self.pythia.getPxFlat()
        pair_py = self.pythia.getPyFlat()
        pair_pz = self.pythia.getPzFlat()

        energy = ak.unflatten(list(pair_e[0]),self._getNumParticlesPerEvent())
        px = ak.unflatten(list(pair_px[0]),self._getNumParticlesPerEvent())
        py = ak.unflatten(list(pair_py[0]),self._getNumParticlesPerEvent())
        pz = ak.unflatten(list(pair_pz[0]),self._getNumParticlesPerEvent())

        return ak.Array(
            {
                'e':energy,
                'px':px,
                'py':py,
                'pz':pz
            }
        )

    def GetProdVertex(self):
        pair_t = self.pythia.getProdTFlat()
        pair_x = self.pythia.getProdXFlat()
        pair_y = self.pythia.getProdYFlat()
        pair_z = self.pythia.getProdZFlat()

        t = ak.unflatten(list(pair_t[0]),self._getNumParticlesPerEvent())
        x = ak.unflatten(list(pair_x[0]),self._getNumParticlesPerEvent())
        y = ak.unflatten(list(pair_y[0]),self._getNumParticlesPerEvent())
        z = ak.unflatten(list(pair_z[0]),self._getNumParticlesPerEvent())

        return ak.Array(
            {
                't':t,
                'x':x,
                'y':y,
                'z':z
            }
        )
    def GetProdVertexStatus(self):
        pair = self.pythia.getHasProdVertexFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetMother1(self):
        pair = self.pythia.getMother1Flat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetMother2(self):
        pair = self.pythia.getMother2Flat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetMotherList(self):
        triple = self.pythia.getMothersFlat()
        tmp = ak.unflatten(list(triple._0),list(triple._1))
        return ak.unflatten(tmp,self._getNumParticlesPerEvent())

    def GetDaughter1(self):
        pair = self.pythia.getDaughter1Flat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetDaughter2(self):
        pair = self.pythia.getDaughter1Flat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetDaughterList(self):
        triple = self.pythia.getDaughtersFlat()
        tmp = ak.unflatten(list(triple._0),list(triple._1))
        return ak.unflatten(tmp,self._getNumParticlesPerEvent())

    def GetColor(self):
        pair = self.pythia.getColorFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    def GetAntiColor(self):
        pair = self.pythia.getAntiColorFlat()
        return ak.unflatten(list(pair[0]), self._getNumParticlesPerEvent())

    # Event-level variable access
    def GetId1Pdf(self):
        return ak.Array(list(self.pythia.getId1PdfArray()))

    def GetId2Pdf(self):
        return ak.Array(list(self.pythia.getId2PdfArray()))

    def GetX1Pdf(self):
        return ak.Array(list(self.pythia.getX1PdfArray()))

    def GetX2Pdf(self):
        return ak.Array(list(self.pythia.getX2PdfArray()))

    def GetPdf1(self):
        return ak.Array(list(self.pythia.getPdf1Array()))

    def GetPdf2(self):
        return ak.Array(list(self.pythia.getPdf2Array()))

    def GetQFac(self):
        return ak.Array(list(self.pythia.getQFacArray()))

    def GetQRen(self):
        return ak.Array(list(self.pythia.getQRenArray()))

    def GetNumMPI(self):
        return ak.Array(list(self.pythia.getNMPIArray()))

    def GetCode(self):
        return ak.Array(list(self.pythia.getCodeArray()))

    def GetAlphaS(self):
        return ak.Array(list(self.pythia.getAlphaSArray()))

    def GetAlphaEM(self):
        return ak.Array(list(self.pythia.getAlphaEMArray()))

    def GetSigmaGen(self):
        return ak.Array(list(self.pythia.getSigmaGenArray()))

    def GetSigmaErr(self):
        return ak.Array(list(self.pythia.getSigmaErrArray()))

    def GetNWeights(self):
        return ak.Array(list(self.pythia.getNWeightsArray()))

    def GetWeights(self):
        return ak.Array(list(self.pythia.getWeightsArray()))

    # Now, methods for handling the HepMC3 mode.
    def WriteHepMC3File(self,filename:str=None):
        if(filename is None):
            filename = self.hepmc_filename

        self.pythia.writeHepMC3File(filename)

    def _print(self,val:Any):
        print('{}: {}'.format(self.print_prefix,val))
        return
