import sys,os
import numpy as np
from utils import PythiaWrapper
import pyhepmc as pyhep

path_prefix = os.getcwd() + '/../../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
import util.particle_selection.truth_selection_nnp as ts_npp
import util.particle_selection.algos as algos

p = PythiaWrapper()
p.ReadString("HardQCD:gg2qqbar = on")
p.SetQuiet(True)
p.SetMPI(False)
p.SetISR(False)
p.SetFSR(False)
p.InitializePythia()

p.Generate()

event = p.GetEvent()
n = p.GetN()

pids = p.GetPdgId()
status = p.GetStatus(hepmc=True)

stats = np.vstack((pids,status)).T
names = p.GetNames()

for i in range(len(names)):
    print('{}: {}, {}'.format(names[i], pids[i], status[i]))

idx = 10

stable_daughters = p.GetStableDaughtersSingle(idx, recursive=True)
print(stable_daughters)

# PyHEP tests
vec = pyhep.FourVector(1.,0.,0.,10.) # px, py, pz, E
pid = 6
status = 1
particle = pyhep.GenParticle(vec,pid,status)

print('particle:',particle)

# Truth selection test.
selection_1 = ts_npp.FirstSelector(status=1,pdgid=22)
selection_2 = ts_npp.FirstSelector(status=4, pdgid=2212)
selection = ts_npp.BasicSelection([selection_1,selection_2])

# HepMC style
print('=== Selecting particle, with HEPMC return style ===')
particles = selection(p,return_hepmc=True)
print(particles)

# Test if particle is stable
for particle in particles:
    stable = algos.IsStable(particle)
    print('Stable:',stable)
    print('Quark:', algos.IsQuark(particle))
    print('Lepton:', algos.IsLepton(particle))
    print('Boson:', algos.IsBoson(particle))
    print('Photon:', algos.IsPhoton(particle))
    print()

    print(particle)

# Numpy style
print('=== Selecting particle, with numpy return style ===')
particles = selection(p,return_hepmc=False)
print(particles)

for particle in particles:
    stable = algos.IsStable(particle, use_hepmc=False)
    print('Stable:',stable)
    print('Quark:', algos.IsQuark(particle, use_hepmc=False))
    print('Lepton:', algos.IsLepton(particle, use_hepmc=False))
    print('Boson:', algos.IsBoson(particle, use_hepmc=False))
    print('Photon:', algos.IsPhoton(particle, use_hepmc=False))
    print()

