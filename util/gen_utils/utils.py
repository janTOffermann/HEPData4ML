import numpy as np
import subprocess as sub
import pyhepmc as hep
from util.particle_selection.algos import IsNeutrino
from util.calcs import PxPyPzEToPtEtaPhiM
from util.qol_utils.pdg import FillPdgHist

# --- Various utility functions for the class --

def _extract_particle_momenta(pythia_wrapper, particle_indices):
    momentum = pythia_wrapper.GetPxPyPzE(particle_indices)
    pdgid = pythia_wrapper.GetPdgId(particle_indices)
    status = pythia_wrapper.GetStatus(particle_indices, hepmc=True)
    return momentum, pdgid, status

def CreateHepMCEvent(pythia_wrapper, particle_indices, event_number):
    hepev = hep.GenEvent()
    hepev.event_number = event_number
    momentum,pdgid,status = _extract_particle_momenta(pythia_wrapper,particle_indices)
    for j, pmu in enumerate(momentum):
        par = hep.GenParticle(momentum=pmu, pid=pdgid[j], status=status[j])
        hepev.add_particle(par)
    return hepev

def AddToHepMCEvent(pythia_wrapper, particle_indices, hepev):
    momentum,pdgid,status = _extract_particle_momenta(pythia_wrapper,particle_indices)
    for j, pmu in enumerate(momentum):
        par = hep.GenParticle(momentum=pmu, pid=pdgid[j], status=status[j])
        hepev.add_particle(par)
    return

def Write2HepMCBuffer(buffername,hepev_list):
    if(type(hepev_list) is not list): hepev_list = [hepev_list]
    with hep.io.WriterAscii(buffername) as f:
        for hepev in hepev_list:
            f.write_event(hepev)
    return

def CopyHepMCBuffer2File(buffername,filename,header=False,footer=False):
    # For speed, we do this using Unix commands (though it seems a bit hacky).

    # We need to determine whether or not to include the HepMC header/footer,
    # which will always be in the buffer file.
    upper_trim = 0 if header else 3
    lower_trim = 0 if footer else 3 # 2 if using Unix head (?)
    comm1 = 'tail -n +{} {}'.format(upper_trim, buffername).split(' ')

    proc = sub.Popen(comm1,stdout=sub.PIPE,text=True)
    if(lower_trim == 0): proc_out = proc.stdout.read().split('\n')
    else: proc_out = proc.stdout.read().split('\n')[:-lower_trim]

    with open(filename,'a') as f:
        f.writelines(line + '\n' for line in proc_out)

def HepMCOutput(hepev_list,buffername,filename,header=False,footer=False):
    Write2HepMCBuffer(buffername,hepev_list) # write the given event(s) to a buffer f ile
    CopyHepMCBuffer2File(buffername,filename,header,footer) # copy buffer file into the full file
    return

# -- plotting functionality below --
def HistogramFinalStateCodes(arr,hist):
    codes = []
    for entry in arr:
        codes.append(entry[-2])
    FillPdgHist(hist,codes)
    return

def HistogramFinalStateKinematics(arr,kinematic_hist_dict):
    # We will make histograms of sum E, E_T and p_T. For each, we consider a sum over all particles,
    # as well as separate sums over invisibles (i.e. neutrinos) and visibles.
    full_sum = np.zeros(4)
    inv_sum = np.zeros(4)
    vis_sum = np.zeros(4)

    for par in arr:
        par_array = np.array([par[0],par[1],par[2],par[3]] )
        full_sum += par_array
        invisible = IsNeutrino(par,use_hepmc=False)
        if(invisible): inv_sum += par_array
        else: vis_sum += par_array

    for vec,key in zip((full_sum,inv_sum,vis_sum),('all','invisible','visible')):
        e = vec[0]
        px = vec[1]
        py = vec[2]
        pz = vec[3]
        vec_cyl = PxPyPzEToPtEtaPhiM([px],[py],[pz],[e])[0]
        pt = vec_cyl[0]
        et = e / np.cosh(vec_cyl[1])
        kinematic_hist_dict[key]['e'].Fill(e)
        kinematic_hist_dict[key]['et'].Fill(et)
        kinematic_hist_dict[key]['pt'].Fill(pt)
