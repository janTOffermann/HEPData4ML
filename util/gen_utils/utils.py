import subprocess as sub
import pyhepmc_ng as hep
from util.hepmc import RestructureParticleArray

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

    # print('\n')
    # for line in proc_out:
    #     print(line)
    # print('\n')

    with open(filename,'a') as f:
        f.writelines(line + '\n' for line in proc_out)

def HepMCOutput(hepev,buffername,filename,loop_number,i,i_real,nevents,n_fail):
    Write2HepMCBuffer(buffername,hepev) # write the given event to a buffer f ile
    CopyHepMCBuffer2File(buffername,filename,loop_number,i,i_real,nevents,n_fail) # copy buffer file into the full file
    return