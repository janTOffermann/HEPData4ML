import numpy as np

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