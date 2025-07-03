import h5py as h5
import numpy as np
import ROOT as rt
from typing import TYPE_CHECKING

if TYPE_CHECKING: # Only imported during type checking -- avoids circular imports we'd otherwise get, since jets imports this file
    from util.post_processing.jets import JetFinder

class TrackCountingBTagging:
    """
    Performs a simple b-tagging algorithm, looking for
    displaced tracks near a jet. Based on the
    "TrackCountingBTagging" module of the Delphes library,
    by M. Selvaggi .
    """

    def __init__(self,mode,track_key, track_pt_min=1.,dr=0.3,track_ip_max=2.,sig_min=6.5,ntracks=3., use_3d=False, tag_name=None):

        self.mode = mode
        self.track_key = track_key
        assert self.mode in ['tag','filter']
        self.tag_name = tag_name

        self.track_pt_min = track_pt_min
        self.dr = dr
        self.track_ip_max = track_ip_max
        self.sig_min = sig_min
        self.ntracks = ntracks
        self.use_3d = use_3d

        # Transient, per-jet variable
        self.tag_status = False

        self.tags = None

        self.print_prefix = '\n\t\TrackCountingBTagging'

    def ModifyInitialization(self,obj):
        """
        This function will modify the initialization so that
        the required track inputs are loaded into memory.
        """

        # Fetch the 4-vector key, and make sure its data is loaded.
        # Note the use of cylindrical coordinates!
        if(self.track_key not in obj.input_collections):
            # Read the necessary keys in.
            f = h5.File(obj.h5_file,'r')
            keys = [x for x in list(f.keys()) if self.track_key in x] # all keys related to the track

            # drop a couple keys we won't need -- probably some very small memory saving
            for keyword in ['charge','pdgid']:
                keys = [x for x in keys if keyword not in x.lower()]

            for key in keys:
                obj.input_collection_arrays[key] = f[key][:]

            f.close()

    def ModifyInputs(self,obj : 'JetFinder'):
        return

    def _tag(self,jet_vector,obj : 'JetFinder'):

        jet = rt.Math.PxPyPzEVector(np.roll(jet_vector,-1)) # np.roll to get from (e,px,py,pz) to (px,py,pz,e)
        count = 0;

        # loop over the track collection
        ntracks = obj.input_collection_arrays['{}.N'.format(self.track_key)][obj._i]

        for i in range(ntracks):

            track_momentum = rt.Math.PxPyPzEVector(np.roll(obj.input_collection_arrays['{}.Pmu'.format(self.track_key)[obj._i]],-1)) # np.roll to get from (e,px,py,pz) to (px,py,pz,e)

            tpt = track_momentum.Pt()
            if(tpt < self.track_pt_min): continue

            d0 = np.abs(obj.input_collection_arrays['{}.D0'.format(self.track_key)][obj._i])
            if(d0 > self.track_ip_max): continue

            dr = rt.Math.VectorUtil.DeltaR(jet,track_momentum)
            if(dr > self.dr): continue

            xd,yd,zd = obj.input_collection_arrays['{}.Xdi'.format(self.track_key)][obj._i]
            dd0 = np.abs(obj.input_collection_arrays['{}.D0.Error'.format(self.track_key)][obj._i])
            z0 = np.abs(obj.input_collection_arrays['{}.DZ'.format(self.track_key)][obj._i])
            dz0 = np.abs(obj.input_collection_arrays['{}.DZ.Error'.format(self.track_key)][obj._i])

            # NOTE: This is all copied quite verbatim from Delphes, but can't I just check if sign > 0, since if not then sip will be negative and always less than self.sig_min? (assuming sig is positive)
            if(self.use_3d):
                sign = np.sign(np.dot([track_momentum.Px(),track_momentum.Py(),track_momentum.Pz()],[xd,yd,zd]))
                # add transverse and longitudinal significances in quadrature
                sip = sign * np.sqrt( np.square(d0/dd0) + np.power(z0/dz0) )
            else:
                sign = np.sign(np.dot([track_momentum.Px(),track_momentum.Py()],[xd,yd]))
                sip = sign * d0 / dd0

            if(sip > self.sig_min): count += 1

        self.tag_status = (count >= self.ntracks)
        return

    def ModifyJets(self, obj : 'JetFinder'):
        """
        This function will tag jets, and fill the corresponding branches.
        """
        import fastjet as fj # NOTE: In practice, fastjet will have been initialized already by JetFinder. Can similarly do this in Softdrop

        self.tags = {i:False for i in obj.jets_dict.keys()}

        for key in obj.jets_dict.keys():
            self._tag(obj.jet_vectors[key],obj) # fills self.tag_status, self.w_candidate
            self.tags[key] = self.tag_status

        if(self.mode=='filter'):
            obj.jet_ordering = [key for key in obj.jet_ordering if self.tags[key]]
            obj._updateJetDictionary()
            # Refresh vectors and constituents -- always need to do this if we filter jets_dict.
            obj._jetsToVectors()
            obj._fetchJetConstituents()

    def ModifyConstituents(self, obj : 'JetFinder'):
            return

    def ModifyWrite(self,obj : 'JetFinder'):
        if(self.mode=='filter'):
            return # do nothing
        else:
            self._initializeBuffer(obj) # will initialize buffer if it doesn't already exist
            self._addFlagToBuffer(obj)
            self._addWToBuffer(obj)

    def _initializeBuffer(self,obj : 'JetFinder'):
        """
        Used if self.mode=='tag', in which case we're writing all jets,
        and including a new branch to indicate whether or not a jet is
        b-tagged.
        """
        self._createBranchNames(obj)

        if(self.tag_name not in obj.buffer.keys()):
            obj.buffer.create_array(self.tag_name,(obj.n_jets_max,),dtype=bool)
        return

    def _createBranchNames(self,obj : 'JetFinder'):
        if(self.tag_name is None):
            self.tag_name = '{}.TrackCountingBTag'.format(obj.jet_name)

    def _addFlagToBuffer(self,obj : 'JetFinder'):
        """
        Adds the b tags to the buffer, for writing.
        Note that the pT sorting of obj is applied,
        which will have been filled by obj._ptSort().
        """
        obj.buffer.set(self.tag_name,obj._i,[self.tags[i] for i in obj.jet_ordering])

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return
