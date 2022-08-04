import numpy as np
import subprocess as sub
import pyhepmc_ng as hep
from util.hepmc import RestructureParticleArray
from util.particle_selection import IsNeutrino
from util.calcs import PxPyPzEToPtEtaPhiM
from util.qol_utils.pdg import pdg_names, pdg_plotcodes, FillPdgHist

# --- Various utility functions for the class --
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
    lower_trim = 3 # 2 if using Unix head (?)
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
