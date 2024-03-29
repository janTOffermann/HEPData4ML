import numpy as np
import ROOT as rt
import pyhepmc as hep

class PythiaROOTWriter:
    """
    A class for writing Pythia output (from our PythiaWrapper class)
    to a ROOT file.
    """
    def __init__(self,outfile='output.root'):
        self.outname = outfile

        self.file = rt.TFile(self.outname,'RECREATE')
        self.InitializeTree()
        self.SetMaxParticles(int(1e4))
        self.SetMaxNDaughters() # will inherit from argument to SetMaxParticles()

    def SetMaxParticles(self,n=10000):
        """
        Maximum number of particles per event that we can save.
        This should be set to be largr than necessary (at the cost
        of more memory usage).
        """
        self.nmax = n

    def SetMaxNDaughters(self,n=None):
        if(n is None): n = self.nmax
        self.ndaughter_max = n

    def InitializeTree(self,treename='CollectionTree'):
        self.file.cd()
        self.tree = rt.TTree(treename,'')
        print(self.tree)

        self.particle_buffer = rt.TObjArray()

        self.buffer = {
            'nparticles':np.zeros(1,dtype=int),
            'pdgid' :np.zeros(self.nmax,dtype=int),
            'status':np.zeros(self.nmax,dtype=int),
            'e'  :np.zeros(self.nmax,dtype=float),
            'px' :np.zeros(self.nmax,dtype=float),
            'py' :np.zeros(self.nmax,dtype=float),
            'pz' :np.zeros(self.nmax,dtype=float),
            'pt' :np.zeros(self.nmax,dtype=float),
            'eta':np.zeros(self.nmax,dtype=float),
            'phi':np.zeros(self.nmax,dtype=float),
            'm'  :np.zeros(self.nmax,dtype=float),
            'rap':np.zeros(self.nmax,dtype=float),
            'theta':np.zeros(self.nmax,dtype=float),
            'et':np.zeros(self.nmax,dtype=float),
            'ndaughters':np.zeros(self.nmax,dtype=int),
            'nmothers'  :np.zeros(self.nmax,dtype=int),
            'daughters':np.zeros((self.nmax,self.ndaughter_max),dtype=int),
            'mothers'  :np.zeros((self.nmax,self.ndaughter_max),dtype=int)
        }

        keys_0d = ['nparticles']
        keys_1d = [key for key,val in self.buffer.items() if val.ndim == 1 and val.shape[0]==self.nmax]
        keys_2d = [key for key,val in self.buffer.items() if val.ndim == 2]

        # Prep the tree -- we store the TObjArray of TParticles, plus flattened info
        # (the latter is easier to access via uproot).
        self.tree.Branch('ParticleArray',self.particle_buffer)

        for key in keys_0d:
            self.tree.Branch(key,self.buffer[key],'{}/I'.format(key)) # TODO: Make the last arg more general? have a func checking dtype

        for key in keys_1d:
            self.tree.Branch(key,self.buffer[key],'{}[nparticles]/D'.format(key))

        for key in keys_2d: # TODO: this is very messy, needs to be rethought. Does it need to be general?
            if('daughters') in key: branch_description = '{}[nparticles][ndaughters]/I'.format(key)
            else: branch_description = '{}[nparticles][nmothers]/I'.format(key)
            # print(key,branch_description)
            self.tree.Branch(key,self.buffer[key],branch_description)

    # def __del__(self): # TODO: Why can't one call the write/close in __del__?
    #     print(self.tree)
    #     self.tree.Write()
    #     self.file.Close()

    def WriteTree(self):
        self.tree.Write()
        self.file.Close()

    def WriteEvent(self,pythia_wrapper):
        """
        Takes in a PythiaWrapper instance.
        TODO: Currently writing a "flat ntuple".
        Would be nice to move to something a bit fancier, e.g. TObjArray of TParticles
        """
        if(pythia_wrapper.GetEvent() is None): pythia_wrapper.GenerateEvent()

        n_particles = pythia_wrapper.GetN()
        pmu         = pythia_wrapper.GetEPxPyPz() # note the signature!
        pmu_cyl     = pythia_wrapper.GetPtEtaPhiM()
        status      = pythia_wrapper.GetStatus()
        pdgid       = pythia_wrapper.GetPdgId()
        mothers     = pythia_wrapper.GetMothers()
        daughters   = pythia_wrapper.GetDaughters()

        # some extra kinematic stuff
        rapidity = pythia_wrapper.GetY()
        et = pythia_wrapper.GetEt()
        theta = pythia_wrapper.GetTheta()

        # mother/daughter info for TParticle -- preserving the Pythia8 conventions
        # see: https://pythia.org/latest-manual/ParticleProperties.html
        daughter1 = pythia_wrapper.GetDaughters1()
        daughter2 = pythia_wrapper.GetDaughters2()
        mother1 = pythia_wrapper.GetMothers1()
        mother2 = pythia_wrapper.GetMothers2()

        # Fill the TParticle buffer.
        self.particle_buffer.Clear()

        for i in range(n_particles):
            p = rt.TParticle()
            p.SetMomentum(*np.roll(pmu[i],-1))
            p.SetPdgCode(int(pdgid[i]))
            p.SetStatusCode(int(status[i]))

            # TODO: For now, we do not store mother/daughter info. TParticle only holds 2 of each, need to keep the Pythia8 style.
            # Note that these values have been adjusted for 0-indexing! To get back to original Pythia8 indices, one mast add 1.
            p.SetFirstDaughter(int(daughter1[i]))
            p.SetLastDaughter(int(daughter2[i]))
            p.SetFirstMother(int(mother1[i]))
            p.SetLastMother(int(mother2[i]))

        # Fill the rest of the buffer.
        self.buffer['nparticles'] = n_particles
        self.buffer['pdgid'] = pdgid
        self.buffer['status'] = status
        self.buffer['e'][:n_particles]  = pmu[:n_particles,0]
        self.buffer['px'][:n_particles] = pmu[:n_particles,1]
        self.buffer['py'][:n_particles] = pmu[:n_particles,2]
        self.buffer['pz'][:n_particles] = pmu[:n_particles,3]
        self.buffer['pt'][:n_particles]  = pmu_cyl[:n_particles,0]
        self.buffer['eta'][:n_particles] = pmu_cyl[:n_particles,1]
        self.buffer['phi'][:n_particles] = pmu_cyl[:n_particles,2]
        self.buffer['m'][:n_particles]   = pmu_cyl[:n_particles,3]
        self.buffer['rap'][:n_particles]   = rapidity[:n_particles]
        self.buffer['et'][:n_particles]    = et[:n_particles]
        self.buffer['theta'][:n_particles] = theta[:n_particles]

        self.buffer['ndaughters'][:n_particles] = [len(x) for x in daughters]
        self.buffer['nmothers'][:n_particles] = [len(x) for x in mothers]

        print('Calling fill. n_particles = ',n_particles)
        self.tree.Fill()

class PythiaHepMC3Writer:
    """
    A class for writing Pythia output (from our PythiaWrapper class)
    to a HepMC file. Slightly different structure than our ROOT writer,
    since the actual I/O is handled a little differently -- here we will
    have a buffer of HepMC events, which will be occasionally flushed
    to a file (so that we don't have too many write operations).
    """
    def __init__(self,outfile='output.root'):
        self.outname = outfile

    def _extract_particle_momenta(pythia_wrapper, particle_indices=None):
        if(particle_indices is None):
            n = pythia_wrapper.GetN()
            particle_indices = np.arange(n)
        momentum = pythia_wrapper.GetPxPyPzE(particle_indices)
        pdgid = pythia_wrapper.GetPdgId(particle_indices)
        status = pythia_wrapper.GetStatus(particle_indices, hepmc=True)
        return momentum, pdgid, status

    def CreateHepMCEvent(self,pythia_wrapper, particle_indices, event_number):
        hepev = hep.GenEvent()
        hepev.event_number = event_number
        momentum,pdgid,status = self._extract_particle_momenta(pythia_wrapper,particle_indices)
        for j, pmu in enumerate(momentum):
            par = hep.GenParticle(momentum=pmu, pid=pdgid[j], status=status[j])
            hepev.add_particle(par)
        return hepev

