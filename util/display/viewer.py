import json, itertools
import ROOT as rt
import numpy as np
import re,pathlib
import h5py as h5
from util.display.setup import DisplaySetup

def ExtractJetRadius(text):
    match = re.search(r'\d+', text)
    return int(match.group()) / 10 if match else None

class EventDisplay:

    def __init__(self):

        # Make sure that our underlying EventDisplay
        # library is built and ready.
        self.setup = DisplaySetup()
        self.setup.Prepare()
        self.delphes_card = None
        self.mode = 0

        self.data = None # dictionary where data from an event is loaded
        self.metadata = None
        self.object_names = []

        # expect to have multiple jet collections, thus have multiple colors to cycle through
        self.jet_colors = [rt.kYellow, rt.kRed, rt.kBlue-7, rt.kSpring, rt.kOrange+1, rt.kMagenta]

    def SetMode(self,mode=0):
        self.mode = mode

    def SetJetPtMin(self,val):
        self.display.GetEventDisplay().SetJetPtMin(val)

    def SetTrackPtMin(self,val):
        self.display.GetEventDisplay().SetTrackPtMin(val)

    def SetTruthParticlePtMin(self,val):
        self.display.GetEventDisplay().SetTruthParticlePtMin(val)

    def InitializeDisplay(self,delphes_card,mode=0):

        self.delphes_card = delphes_card
        self.display = rt.Display.DisplayInterface(self.delphes_card)
        self.display.GetEventDisplay().SetColorMode(mode)
        # self.display.SetColorMode(mode)

        # this_dir = os.path.dirname(os.path.abspath(__file__))
        # filename = this_dir + '/test_data/atlas.root'
        # self.display.Test(filename)

    def Display(self, input_file, event_index = 0, delphes_card = None):
        if(delphes_card is not None):
            self.delphes_card = delphes_card
        assert self.delphes_card is not None
        self.input_file = input_file

        # load in data -- for now we don't allow for event changing, just picking one event to view.
        self.LoadData(event_index)

        self.display.DisplayEvent(self.delphes_card)

    def FindObjectNames(self):
        self.object_names = sorted(list(set([key.split('.')[0] for key in self.data.keys()])))


    def LoadData(self,event_index):
        assert self.input_file is not None
        assert pathlib.Path(self.input_file).exists()

        f = h5.File(self.input_file,'r')
        self.data = {key: f[key][event_index] for key in f.keys()}
        self.metadata = {key:f.attrs[key] for key in f.attrs.keys()}
        f.close()

        self.FindObjectNames()

        self.LoadVertexData()
        self.LoadTrackData()
        self.LoadCaloData()
        self.LoadLeptonData()
        self.LoadPhotonData()
        self.LoadJetData()
        self.LoadMETData()
        self.LoadGenParticleData()
        self.LoadJetConstituentData()

    def LoadTrackData(self):
        for object_name in self.object_names:
            if('Track' in object_name):
                nobj = self.data['{}.N'.format(object_name)]
                self.display.GetEventDisplay().AddTrackData(
                    object_name,
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,3],
                    self.data['{}.Xdi'.format(object_name)][:nobj,0],
                    self.data['{}.Xdi'.format(object_name)][:nobj,1],
                    self.data['{}.Xdi'.format(object_name)][:nobj,2],
                    self.data['{}.PdgId'.format(object_name)][:nobj]
                )

    def LoadCaloData(self):
        """
        We dump all the calorimeter data into the same collection for the viewer,
        due to some limitations with how it currently functions.
        """
        name = 'UnifiedCaloData'
        for object_name in self.object_names:
            if('{}.Edges.Eta'.format(object_name) in self.data.keys()):
                nobj = self.data['{}.N'.format(object_name)]

                self.display.GetEventDisplay().AddCaloData(
                    name,
                    self.data['{}.Edges.Eta'.format(object_name)][:nobj,0],
                    self.data['{}.Edges.Eta'.format(object_name)][:nobj,1],
                    self.data['{}.Edges.Phi'.format(object_name)][:nobj,0],
                    self.data['{}.Edges.Phi'.format(object_name)][:nobj,1],
                    self.data['{}.E.EM'.format(object_name)][:nobj],
                    self.data['{}.E.Hadronic'.format(object_name)][:nobj]
                )

    def LoadLeptonData(self):
        for object_name in self.object_names:
            if('Electron' in object_name):
                nobj = self.data['{}.N'.format(object_name)]
                self.display.GetEventDisplay().AddElectronData(
                    object_name,
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,3],
                    self.data['{}.Xdi'.format(object_name)][:nobj,0],
                    self.data['{}.Xdi'.format(object_name)][:nobj,1],
                    self.data['{}.Xdi'.format(object_name)][:nobj,2],
                    self.data['{}.Charge'.format(object_name)][:nobj]
                )
            elif('Muon' in object_name):
                nobj = self.data['{}.N'.format(object_name)]
                self.display.GetEventDisplay().AddMuonData(
                    object_name,
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,3],
                    self.data['{}.Xdi'.format(object_name)][:nobj,0],
                    self.data['{}.Xdi'.format(object_name)][:nobj,1],
                    self.data['{}.Xdi'.format(object_name)][:nobj,2],
                    self.data['{}.Charge'.format(object_name)][:nobj]
                )

    def LoadJetData(self):
        jet_counter = 0
        for object_name in self.object_names:
            if('jet' in object_name.lower()):
                nobj = self.data['{}.N'.format(object_name)]
                radius = ExtractJetRadius(object_name)
                if(radius is None): radius = 0.4

                color = self.jet_colors[jet_counter % len(self.jet_colors)]

                self.display.GetEventDisplay().AddJetData_EPxPyPz(
                    object_name,
                    self.data['{}.Pmu'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu'.format(object_name)][:nobj,3],
                    radius,
                    color
                )
                jet_counter += 1

    def LoadPhotonData(self):
        for object_name in self.object_names:
            if('photon' in object_name.lower() and 'eflow' not in object_name.lower()):
                nobj = self.data['{}.N'.format(object_name)]
                self.display.GetEventDisplay().AddPhotonData(
                    object_name,
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,3],
                )

    def LoadMETData(self):
        # TODO: Our MissingET objects can apparently have non-zero z-components...?
        #       Enforcing eta=0 for the visualization
        for object_name in self.object_names:
            if('missinget' in object_name.lower()):
                nobj = 1
                self.display.GetEventDisplay().AddMETData(
                    object_name,
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,0],
                    np.zeros(nobj),
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu_cyl'.format(object_name)][:nobj,3],
                )

    def LoadGenParticleData(self, enforce_lifetimes=False):
        for object_name in self.object_names:
            if('truthparticle' in object_name.lower()):
                nobj = self.data['{}.N'.format(object_name)]

                # Some trickery: we always have a StableTruthParticles collection,
                # which -- as the name implies -- is stable, so it does not contain
                # branches like "Decay.Xmu" or "Stable" because these would be
                # redundant. Other truth particle collections will have these, though.
                # Whether or not we want to use these is another thing -- for short-lived
                # truth particles, limiting their visible tracks to their actual lifetime
                # might be unproductive as they will be effectively invisible in the event
                # display.

                if(object_name.lower() == 'stabletruthparticles' or not enforce_lifetimes):
                    decay_xmu = np.full((nobj,4),0.)
                    stable = np.full(nobj,True)
                else: # realistic lifetimes -- at the possible cost of infinitesimally short tracks
                    decay_xmu = self.data['{}.Decay.Xmu'.format(object_name)]
                    stable = self.data['{}.Stable'.format(object_name)]

                self.display.GetEventDisplay().AddGenParticleData(
                    object_name,
                    self.data['{}.Pmu'.format(object_name)][:nobj,0],
                    self.data['{}.Pmu'.format(object_name)][:nobj,1],
                    self.data['{}.Pmu'.format(object_name)][:nobj,2],
                    self.data['{}.Pmu'.format(object_name)][:nobj,3],
                    self.data['{}.Production.Xmu'.format(object_name)][:nobj,1],
                    self.data['{}.Production.Xmu'.format(object_name)][:nobj,2],
                    self.data['{}.Production.Xmu'.format(object_name)][:nobj,3],
                    decay_xmu[:nobj,1],
                    decay_xmu[:nobj,2],
                    decay_xmu[:nobj,3],
                    stable[:nobj],
                    self.data['{}.PdgId'.format(object_name)][:nobj]
                )

    def LoadVertexData(self):
        for object_name in self.object_names:
            if('vertex' in object_name.lower()):
                nobj = self.data['{}.N'.format(object_name)]
                self.display.GetEventDisplay().AddVertexData(
                    object_name,
                    self.data['{}.Xmu'.format(object_name)][:nobj,1],
                    self.data['{}.Xmu'.format(object_name)][:nobj,2],
                    self.data['{}.Xmu'.format(object_name)][:nobj,3]
                )


    def LoadJetConstituentData(self):
        jet_counter = 0

        metadata_index = self.data['Metadata.JetCollections.InputCollections']
        jet_name_dict = json.loads(self.metadata['Metadata.JetCollections.InputCollections'][metadata_index])
        print(jet_name_dict)

        for object_name in self.object_names:
            if('jet' in object_name.lower()):
                nobj = self.data['{}.N'.format(object_name)]
                radius = ExtractJetRadius(object_name)
                if(radius is None): radius = 0.4
                color = self.jet_colors[jet_counter % len(self.jet_colors)]

                # map the constituents
                constituent_collections = self.data['{}.Constituents.Collection'.format(object_name)][:nobj]
                constituent_collection_indices = self.data['{}.Constituents.Collection.Index'.format(object_name)][:nobj]

                # go from constituent_collections to collection names
                print('object_name = {}, constituent_collections = '.format(object_name))
                for entry in constituent_collections:
                    print('\t',entry)
                print('jet_name_dict[object_name] = ',jet_name_dict[object_name])

                constituent_collection_names = [
                    [jet_name_dict[object_name][x] for x in y]
                    for y in constituent_collections
                ]

                # Now we construct eta and phi edges, and energy deposits, which will
                # go into a Lego-style calorimeter display.
                # The caveat is that things like EFlowTrack don't have eta/phi edges!
                eta_edges = []
                phi_edges = []
                e_em = []
                e_had = []
                for i in range(nobj):
                    names = constituent_collection_names[i]
                    mask = ['{}.Edges.Eta'.format(name) in self.data.keys() for name in names]
                    names = list(itertools.compress(names,mask))
                    tower_indices = list(itertools.compress(constituent_collection_indices[i],mask))

                    eta_edges_i = np.array([self.data['{}.Edges.Eta'.format(name)][x] for name,x in zip(names,tower_indices)])
                    phi_edges_i = np.array([self.data['{}.Edges.Phi'.format(name)][x] for name,x in zip(names,tower_indices)])
                    eem_i = np.array([self.data['{}.E.EM'.format(name)][x] for name,x in zip(names,tower_indices)])
                    ehad_i = np.array([self.data['{}.E.Hadronic'.format(name)][x] for name,x in zip(names,tower_indices)])

                    # Add these with some new function
                    self.display.GetEventDisplay().AddJetConstituentCaloData(
                        object_name,
                        eta_edges_i[:nobj,0],
                        eta_edges_i[:nobj,1],
                        phi_edges_i[:nobj,0],
                        phi_edges_i[:nobj,1],
                        eem_i[:nobj],
                        ehad_i[:nobj]
                    )

    def LoadCaloData2(self):
        """
        We dump all the calorimeter data into the same collection for the viewer,
        due to some limitations with how it currently functions.
        """
        name = 'UnifiedCaloData'
        for object_name in self.object_names:
            if('{}.Edges.Eta'.format(object_name) in self.data.keys()):
                nobj = self.data['{}.N'.format(object_name)]

                self.display.GetEventDisplay().AddCaloData(
                    name,
                    self.data['{}.Edges.Eta'.format(object_name)][:nobj,0],
                    self.data['{}.Edges.Eta'.format(object_name)][:nobj,1],
                    self.data['{}.Edges.Phi'.format(object_name)][:nobj,0],
                    self.data['{}.Edges.Phi'.format(object_name)][:nobj,1],
                    self.data['{}.E.EM'.format(object_name)][:nobj],
                    self.data['{}.E.Hadronic'.format(object_name)][:nobj]
                )
