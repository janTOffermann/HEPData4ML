import sys,os
import numpy as np
from utils import PythiaWrapper
import pyhepmc as pyhep

path_prefix = os.getcwd() + '/../../'
if(path_prefix not in sys.path): sys.path.append(path_prefix)
import util.particle_selection.particle_selection as par_sel
import util.particle_selection.algos as algos
import util.particle_selection.selection_algos as salgos

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
daughters = p.GetDaughters()

for i in range(len(names)):
    print(i, '{}: {}, {}'.format(names[i], pids[i], status[i]),daughters[i])

idx = 10

stable_daughters = p.GetStableDaughtersSingle(idx, recursive=True)
print('For particle {}, the stable daughters:'.format(idx))
for d in stable_daughters:
    print('\t{}'.format(d))


print('\n=== PyHEP test. === \n')

# PyHEP tests
vec = pyhep.FourVector(1.,0.,0.,10.) # px, py, pz, E
pid = 6
status = 1
particle = pyhep.GenParticle(vec,pid,status)
print('particle:',particle)

# Particle selection tests.
print('\n\n Testing basic particle selection: FirstSelector & BasicSelection.')
selection_1 = par_sel.FirstSelector(status=1,pdgid=22)
selection_2 = par_sel.FirstSelector(status=4, pdgid=2212)
selection = par_sel.BasicSelection([selection_1,selection_2])

selected_indices = selection(p)

print('Selected indices:')
for entry in selected_indices:
    print(entry)

selection_1.Print()

print('\n\n Testing some selection algorithms directly.')
i = 1
print('Using particle {}, which is {} ({})'.format(i,pids[i],names[i]))

print('\n\nTesting GatherQuarks.')
gatherer = algos.GatherQuarks()
result = gatherer(p,i)
print('Quark daughters:')
for k in result:
    print('\t [{}] ({})'.format(k,names[k]))

print('\n\nTesting GatherStableDaughters.')
gatherer = algos.GatherStableDaughters()
result = gatherer(p,i)
print('Stable daughters:')
for k in result:
    print('\t [{}] ({})'.format(k,names[k]))

print('\n\nTesting AlgoSelection with SelectFinalStateDaughters.')
algo = salgos.SelectFinalStateDaughters(truth_selection=selection) # use selection from earlier
n = 10
selector = par_sel.AlgoSelection(algo,n)
result = selector(p)
print('Result of AlgoSelection:')
for k in result:
    print('\t [{}] ({})'.format(k,names[k]))

print('\n\nTesting MultiSelection with the above AlgoSelection + FirstSelector.')
selectors = [selector, selection_2]
selector = par_sel.MultiSelection(selectors)
result = selector(p)
print('Result of MultiSelection:')
for k in result:
    print('\t [{}] ({})'.format(k,names[k]))

print('\n\nTesting event-level info.')
process_id = p.GetProcessCode()
process_name = p.GetProcessName()

print('Process code: {}\n\t({})'.format(process_id, process_name))

