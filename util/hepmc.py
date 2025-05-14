import numpy as np
import pyhepmc as hep
import subprocess as sub
import pathlib

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

def _extract_particle_momenta(pythia_wrapper, particle_indices=None):
    momentum = pythia_wrapper.GetPxPyPzE(particle_indices)
    pdgid = pythia_wrapper.GetPdgId(particle_indices)
    status = pythia_wrapper.GetStatus(particle_indices, hepmc=True)
    return momentum, pdgid, status

def CreateHepMCEvent(pythia_wrapper, particle_indices, event_number):
    """
    This creates a HepMC file containing GenParticle objects. It does not
    include any vertex information, which is how HepMC keeps track of particles'
    mothers/daughters. Thus the resulting file does not keep track of the particles'
    relationships with one another! This is partly necessitated as particle_indices
    may be filtering out some particles from the resulting file, and it would be
    strange (maybe even non-functional) to have vertices referencing particles not
    present in the final file.
    """
    hepev = hep.GenEvent()
    hepev.event_number = event_number
    momentum,pdgid,status = _extract_particle_momenta(pythia_wrapper,particle_indices)
    for j, pmu in enumerate(momentum):
        par = hep.GenParticle(momentum=pmu, pid=pdgid[j], status=status[j])
        hepev.add_particle(par)
    return hepev

def CreateFullHepMCEvent(pythia_wrapper, event_number):
    """
    Together with the Pythia8 wrapper (another part of this package),
    this basically gives a (Pythonic) Pythia8->HepMC3 (ascii) interface.
    """
    hepev = hep.GenEvent()
    hepev.event_number = event_number
    momentum,pdgid,status = _extract_particle_momenta(pythia_wrapper) # fetch all the momenta
    mother_indices = pythia_wrapper.GetMothers()
    prod = pythia_wrapper.GetProd() # production vertices' positions (t,x,y,z)

    # Fetch the particles.
    particles = [hep.GenParticle(momentum=momentum[i], pid=pdgid[i], status=status[i]) for i in range(len(momentum))]
    # for par in particles: hepev.add_particle(par) # this line broke things in the weirdest way, led to inconsistencies with # of vertices & disconnected graphs

    # Now add the vertices, based on mother/daughter information.
    # Note that we do not want to create duplicate vertices, which will
    # happen if we naively loop over all particles' production vertices.
    # TODO: Maybe this method is more memory/time-intensive than it needs to be?
    #       I might not be thinking about this right, and there's a simpler
    #       approach (with less looping?).
    vertices = []
    used_mothers = None # will store tuples of (vertex index, mother index)
    for i in range(len(momentum)):

        mother_idxs = mother_indices[i]
        nmothers = len(mother_idxs)
        if(nmothers == 0): continue

        new_vertex = True
        used_mother_idx = None
        if(used_mothers is not None):
            for j in range(nmothers):
                if(mother_idxs[j] in used_mothers[:,-1]):
                    new_vertex = False
                    used_mother_idx = mother_idxs[j]
                    break

        if(new_vertex):
            vertex = hep.GenVertex(np.roll(prod[i],-1)) # converting from (t,x,y,z) to (x,y,z,t) for constructor
            # add stuff
            vertex.add_particle_out(particles[i])
            for idx in mother_idxs:
                vertex.add_particle_in(particles[idx])
                mother_tuple = np.array([len(vertices),idx],dtype=int)
                if(used_mothers is None): used_mothers = np.array([mother_tuple])
                else: used_mothers = np.append(used_mothers,[mother_tuple],axis=0)
            vertices.append(vertex)
        else:
            # find the existing vertex
            vertex_idx = used_mothers[np.where(used_mothers[:,-1] == used_mother_idx)[0][0],0]
            vertices[vertex_idx].add_particle_out(particles[i])

        hepev.add_vertex(vertex)

    # print('Wrote {} vertices'.format(len(vertices)))
    # print('\nChecking vertices.')
    # print('{} vertices in hepev.'.format(len(hepev.vertices)))
    # for i,vertex in enumerate(hepev.vertices):
    #     print('[{}]: {}\tid = {}'.format(i,vertex,vertex.id))
    #     for p in vertex.particles_in:
    #         print('\t-> {}'.format(p))
    #     for p in vertex.particles_out:
    #         print('\t<- {}'.format(p))

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

def ExtractHepMCEvents(files,get_nevents=False, silent=False):
    events = []
    nevents = 0
    for file in files:
        # Check that the file exists.
        if(not pathlib.Path(file).exists()):
            if(not silent):
                print('Warning: Tried to access file {} but it does not exist!'.format(file))
            continue

        with hep.io.ReaderAscii(file) as f:
            for i,evt in enumerate(f):
                events.append(evt)
                if(get_nevents): nevents += 1

    if(get_nevents): return events, nevents
    return events
