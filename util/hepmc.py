import numpy as np

def RestructureParticleArray(arr):
    fields = ['px','py','pz','E','pdgid','status','eta','phi']
    arr = np.array([[x[y] for y in fields] for x in arr])
    return arr

def HepMCArrayToNumpy(arr):
    npar = len(arr)
    # Ordering mimics the ordering of numpythia return_hepmc=False, with some fields removed
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
    for i,particle in enumerate(arr):
        arr2[i] = (particle.e,particle.px,particle.py,particle.pz,particle.eta,particle.phi,particle.pid,particle.status)
    return arr2