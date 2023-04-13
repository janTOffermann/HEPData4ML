import numpy as np
import pyhepmc as hep
import subprocess as sub

def RestructureParticleArray(arr, fields=None):
    if(fields is None):
        fields = ['px','py','pz','E','pdgid','status','eta','phi']
    arr = np.array([[x[y] for y in fields] for x in arr])
    return arr

def HepMCArrayToNumpy(arr, short=False):
    npar = len(arr)
    # Ordering mimics the ordering of numpythia return_hepmc=False. Matching all fields though in practice we don't use them all.
    dtypes = {
        'E':'f8',
        'px':'f8',
        'py':'f8',
        'pz':'f8',
        'pT':'f8',
        'mass':'f8',
        'rap':'f8',
        'eta':'f8',
        'theta':'f8',
        'phi':'f8',
        'prodx':'f8',
        'prody':'f8',
        'prodz':'f8',
        'prodt':'f8',
        'pdgid':'i4',
        'status':'i4'
    }

    if(short): # TODO: Keeping old functionality -- not sure if there is any benefit to this.
        dtypes = {
            'E':'f8',
            'px':'f8',
            'py':'f8',
            'pz':'f8',
            'eta':'f8',
            'phi':'f8',
            'pdgid':'i4',
            'status':'i4'
        }
    names = []
    types = []
    for key,val in dtypes.items():
        names.append(key)
        types.append(val)
    arr2 = np.zeros((npar,),dtype={'names':names,'formats':types})
    for i,p in enumerate(arr):
        v = [p.e,p.px,p.py,p.pz] # E,px,py,pz
        if(short): arr2[i] = (*v,p.eta,p.phi,p.pid,p.status)
        else:
            prod = np.zeros(4) # TODO: Setting prodx,prody,prodz,prodt = 0 (GenParticle doesn't have these attributes?). We currently are not using these anyway, but keep this in mind.
            arr2[i] = (*v,p.pt,p.mass,p.rap,p.eta,p.theta,p.phi,*prod,p.pid,p.status)
    return arr2

# Some utility functions -- partly for interfacing with our Pythia8 wrapper.

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

def CompressHepMC(files, delete=True, cwd=None):
    for file in files:
        compress_file = file.replace('.hepmc','.tar.bz2')
        if(cwd is not None): compress_file = '{}/{}'.format(cwd,compress_file)
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],file.split('/')[-1]]
        # if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)
        if(delete):
            sub.check_call(['rm',file.split('/')[-1]],shell=False,cwd=cwd)
    return
