import numpy as np
import util.jet_selection as jetsel
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos
import util.event_selection as eventsel
import util.event_filter as eventfilter
import util.post_processing.tracing as tracing
import util.post_processing.truth_sum as truth_sum
import util.post_processing.rotate as rotate
import util.post_processing.containment as containment
import util.post_processing.jh_tagging as jh_tagging
import config.selectors as s

selections = s.selections

config = {
    'proc' : 'Top_Wqq', # Name must correspond to a card name (without extension) in util/pythia_templates
    'hadronization' : True, # Pythia8 hadronization flag
    'mpi' : False, # Pythia8 multi-parton interactions flag
    'isr' : False, # Pythia8 initial-state radiation flag
    'fsr' : False, # Pythia8 final-state radiation flag
    'rng' : 1, # Pythia8 RNG seed
    'delphes' : True, # Whether or not to use Delphes
    'delphes_card' : None, # path to the Delphes card to use. If None, will use the ATLAS Delphes card that ships with Delphes
    'delphes_dir' : '/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre05-775ca/x86_64-centos7-gcc11-opt' ,# directory containing the Delphes installation. If None, will be build in a local directory "delphes".
    'fastjet_dir' : '/home/jaofferm/HEPData4ML/fastjet',
    'jet_radius': 0.8, # jet radius in (eta,phi). Jets currently all use anti-kt algorithm.
    'jet_min_pt': 15., # minimum pT cut, GeV
    'jet_max_eta': 2., # absolute value eta cut
    'jet_n_par': 200, # max number of jet constituents to save per jet (pT-ordered)
    'n_truth' : 3 + 2 + 120, # max number of truth particles to save per jet
    # 'event_filter' : None, # if not using any event filter
    'event_filter' : eventfilter.MultiFilter( # multi-filter -- applies filters (cuts) sequentially, if an event fails a filter we generate another to replace it
        [
            eventfilter.PtMatchedJetFilter(0.8,selections['t->Wb'],pt_window_frac=0.005,pt_min_jet=15.,eta_max_jet=2.), # require jet pT to match top pT closely
            eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['bqq'],matching_radius=0.6,pt_min_jet=15.,eta_max_jet=2.), # require W(qq) & b within 0.6 of jet centroid
            eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['W daughters'],pt_threshold_sum=2.5,pt_min_jet=15.,eta_max_jet=2.), # require stable W daughters within jet -- allow sum of uncontained daughters to carry off up to 2.5 GeV in pT
            eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['b daughters'],pt_threshold_sum=2.5,pt_min_jet=15.,eta_max_jet=2.) # require stable b daughters within jet -- allow som of uncontained daughters to carry off up to 2.5 GeV in pT
        ]
    ),
    # 'event_filter_flag': None, # an event_filter algorithm, but instead of filtering out events it will simply add a boolean flag to the dataset
    'event_filter_flag' : eventfilter.FilterFlag(
        eventfilter.MultiFilter( # multi-filter -- applies filters (cuts) sequentially, if an event fails a filter we generate another to replace it
            [
                eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['W daughters'],pt_threshold_sum=2.5,pt_min_jet=15.,eta_max_jet=2.,matching_radius=0.8), # require stable W daughters within jet -- allow sum of uncontained daughters to carry off up to 2.5 GeV in pT
                eventfilter.ContainedJetFilter(0.8,selections['t->Wb'],selections['b daughters'],pt_threshold_sum=2.5,pt_min_jet=15.,eta_max_jet=2.,matching_radius=0.8) # require stable b daughters within jet -- allow som of uncontained daughters to carry off up to 2.5 GeV in pT
            ]
        ),
        name = 'jet_top_daughters_contained_dR08'
    ),
    'truth_selection' : selections['t->Wb w/ qq and W daughters'],
    'final_state_selection': parsel.AlgoSelection(algos.SelectFinalState(),200),
    'event_selection' : eventsel.TruthDistanceSelection(distance=2.4, n_truth=3), # filters an event to remove some final-state particles, primarily for lowering HepMC file-size. Ideally does not affect the final output.
    'invisibles' : False, # if False, invisible particles in the final-state selection (neutrinos) will be discarded. If true, they may be input to jet clustering!
    'jet_selection':jetsel.GetNearestJet(truth_code=6,max_dr=0.8),
    'signal_flag' : 1, # What to provide as the "signal_flag" for these events. (relevant if combining datasets). Must be >= 0.
    'post_processing': [
        tracing.Tracer(verbose=True), # Delphes tracing, computes the "W daughteriness" of Delphes detector towers
        jh_tagging.JHTagger(2,120,verbose=True), # Johns Hopkins top tagger -- W is in positiong 2 (used for computing prediction's resolution), gives up to 120 constituents of candidate W subjet
        truth_sum.TruthParticleSum(np.arange(5,125),0.8,verbose=True), # produce 4-vector sum of W daughters within the jet cone
        rotate.Rotation(2,1,['Pmu'],verbose=True), # create rotated copies of 4-vectors
        containment.Containment([3,4],  jet_distance=0.8,containment_key='jet_qq_contained_dR08',  max_dr_key='jet_qq_dr_max', verbose=True), # compute W containment, using jet radius
        containment.Containment([1,3,4],jet_distance=0.8,containment_key='jet_bqq_contained_dR08',max_dr_key='jet_bqq_dr_max',verbose=True), # compute top containment, using jet radius
        containment.Containment([3,4],  jet_distance=0.6,containment_key='jet_qq_contained_dR06',  max_dr_key=None, verbose=True), # compute W containment, using 0.75 * jet radius
        containment.Containment([1,3,4],jet_distance=0.6,containment_key='jet_bqq_contained_dR06',max_dr_key=None,verbose=True), # compute top containment, using 0.75 * jet radius
    ],
    'record_final_state_indices' : True, # Whether or not to record jet constituents' indices w.r.t. the order they were passed to jet clustering (order of particles in HepMC file, or order of Delphes objects if using Delphes).
    'split_seed' : 1, # seed to be used for the RNG when splitting dataset into train/test/validation samples
    'use_vectorcalcs' : False # whether or not to use the VectorCalcs C++/ROOT library (which is part of this repo). May speed up some computations, but can lead to issues if there's a problem with the ROOT build (e.g. CVMFS/LCG_103 seems to cause issues)
}

# Don't adjust the lines below.
try: config['truth_selection'].SetHadronization(config['hadronization'])
except: pass
