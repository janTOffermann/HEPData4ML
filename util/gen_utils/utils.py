import sys, glob, uuid
import ROOT as rt
import numpy as np
import subprocess as sub
import pyhepmc_ng as hep
from util.fastjet import BuildFastjet, ParticleInfo
from util.config import GetJetConfig
from util.calcs import DeltaR2Vectorized
from util.hepmc import RestructureParticleArray

# --- FASTJET IMPORT ---
# TODO: Can this be done more nicely?
fastjet_dir = BuildFastjet(j=8)
fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]
if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
import fastjet as fj
# ----------------------

def ClusterJets(arr,jetdef):
    pseudojets = []
    N = len(arr)
    for j in range(N):
        # Create a pseudojet object and add it to the list.
        # Signature is px, py, pz, E.
        vec = [arr[j][x] for x in ['px','py','pz','E']]
        pj = fj.PseudoJet(*vec)
        pjinfo = ParticleInfo(j, arr[j]['status'], arr[j]['pdgid'])
        pj.set_python_info(pjinfo)
        pseudojets.append(pj)
    jets = jetdef(pseudojets)
    del pseudojets # TODO: Is this okay? Is it useful to delete?
    return jets

def ApplyJetCuts(jets):
    jet_config = GetJetConfig()

    # Apply minimum jet pT cut.
    jet_pt = np.array([jet.pt() for jet in jets])
    jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[jet_pt >= jet_config['jet_min_pt']]
    jets = [jets[i] for i in jet_indices]

    # Apply optional maximum |eta| cut.
    jet_eta = np.array([jet.eta() for jet in jets])
    jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[np.abs(jet_eta) <= jet_config['jet_max_eta']]
    jets = [jets[i] for i in jet_indices]
    return jets

def CreateHepMCEvent(particle_array,event_number):
    particle_array = RestructureParticleArray(particle_array)
    N = len(particle_array)
    hepev = hep.GenEvent()
    hepev.event_number = event_number
    for j in range(N):
        momentum = particle_array[j][:4] # px, py, pz, E.
        pid = int(particle_array[j][4])
        status = int(particle_array[j][5])
        par = hep.GenParticle(momentum=momentum, pid=pid, status=status)
        hepev.add_particle(par)
    return hepev

def Write2HepMCBuffer(buffername,hepev):
    with hep.WriterAscii(buffername) as f:
        f.write_event(hepev)
    return

def CopyHepMCBuffer2File(buffername,filename,loop_number,i,i_real,nevents,n_fail):
    # For speed, we do this using Unix commands (though it seems a bit hacky).
    upper_trim = 3
    lower_trim = 3 # 2 if using Unix head
    if(loop_number == 0 and i_real == 1): upper_trim = 0
    elif(i == nevents-1 and n_fail == 0): lower_trim = 0
    comm1 = 'tail -n +{} {}'.format(upper_trim, buffername).split(' ')

    proc = sub.Popen(comm1,stdout=sub.PIPE,text=True)
    if(lower_trim == 0): proc_out = proc.stdout.read().split('\n')
    else: proc_out = proc.stdout.read().split('\n')[:-lower_trim]

    with open(filename,'a') as f:
        f.writelines(line + '\n' for line in proc_out)

def HepMCOutput(hepev,buffername,filename,loop_number,i,i_real,nevents,n_fail):
    Write2HepMCBuffer(buffername,hepev) # write the given event to a buffer f ile
    CopyHepMCBuffer2File(buffername,filename,loop_number,i,i_real,nevents,n_fail) # copy buffer file into the full file
    return

# TODO: Make this a callable
def TruthDistanceSelection(arr, arr_truth, alpha):
    # Keep final-state particles within a certain distance of the selected truth particles.
    # We should be generous with our selection to avoid throwing out things that
    # have any chance of making their way into our jets later on. This is determined
    # by the alpha parameter, which should always be greater than 1.
    jet_config = GetJetConfig()

    arr_truth = RestructureParticleArray(arr_truth)
    arr = RestructureParticleArray(arr)
    dr2_limit = np.square(alpha * jet_config['jet_radius'])

    distances = DeltaR2Vectorized(arr[:,-2:], arr_truth[:,-2:])
    min_distances = np.min(distances,axis=1)
    selected = (np.abs(min_distances) < np.abs(dr2_limit))
    arr = arr[selected]

    status = len(arr) > 0

    return status, arr,arr_truth

def JetConstituentSelection(arr, arr_truth, jetdef, max_dr):
    jet_config = GetJetConfig()

    # Perform jet clustering on final-state particles.
    jets = ClusterJets(arr,jetdef)

    # Apply optional minimum jet pT cut, and maximum |eta| cut.
    jets = ApplyJetCuts(jets)

    if(len(jets) == 0): return False,[],[]

    # Get our selected jet.
    selected_jet_idx = jet_config['jet_selection'](truth=arr_truth, jets=jets, max_dr=max_dr)
    if(selected_jet_idx < 0): # No jet passed selection.
        return False,[],[]

    jet = jets[selected_jet_idx]
    jet_constituents = np.array([
        [x.px(), x.py(), x.pz(), x.E(), x.python_info().pdg_id, x.python_info().status]
        for x in jet.constituents()]
        )
    arr = jet_constituents
    arr_truth = RestructureParticleArray(arr_truth)
    return True,arr,arr_truth

