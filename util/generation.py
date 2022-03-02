import sys, os, glob
import subprocess as sub
import numpy as np
import pyhepmc_ng as hep
import awkward as ak
import numpythia as npyth # Pythia, hepmc_write
from util.config import GetPythiaConfig, GetTruthSelection, GetJetConfig
from util.fastjet import BuildFastjet, ParticleInfo
from util.calcs import DeltaR2Vectorized
import util.qol_util as qu

# Generate a bunch of events in the given pT range,
# and save them to a HepMC file.
# Note that we do not perform event selection: All generated events are saved.
def GenerateSimple(nevents, pt_min, pt_max, filename = 'events.hepmc'):

    pythia_config = GetPythiaConfig(pt_min = pt_min, pt_max = pt_max)
    pythia = npyth.Pythia(params = pythia_config)

    prefix = 'Generating events for pT bin [{},{}]:'.format(pt_min,pt_max)
    suffix = 'Complete'
    i = 0
    qu.printProgressBarColor(i,nevents, prefix=prefix, suffix=suffix, length=50)

    for i,event in enumerate(npyth.hepmc_write(filename, pythia(events=nevents))):
        qu.printProgressBarColor(i+1,nevents, prefix=prefix, suffix=suffix, length=50)
    return

def GenerationLoop(pythia, nevents,
                   sel1, sel_truth,
                   filename,
                   prefix, suffix, bl=50,
                   i_real = 1, nevents_disp=None, fastjet_dir=None,
                   loop_number=0, chunk_size=100,
                   alpha = 2., max_dr = 0.2):

    use_delphes = True

    if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
    import fastjet as fj
    jet_config = GetJetConfig()
    jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])

    n_fail = 0
    if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

    # The way that pyhepmc_ng's WriterAscii works, writing an event will overwrite the whole file.
    # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
    # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
    # down things, so we ultimately want to find some way to do a write with "append" functionality.

    buffername = filename.replace('.hepmc','_buffer.hepmc')

    for i,event in enumerate(pythia(events=nevents)):
        success = True

        # Get the selected truth particles.
        arr_truth = [event.all(selection=x, return_hepmc=False) for x in sel_truth]

        for truth in arr_truth:
            if(len(truth) == 0): # missing some truth particle -> potential trouble?
                n_fail += 1
                success = False
                break

        if(not success): continue

        # Get the truth-level particles.
        # TODO: In some cases of t -> b W, using GenEvent.all() yields two b-quarks (throwing off some later code).
        #       This happens without MPI/ISR/FSR/ and both of them have the same generator status. How?
        arr_truth = np.concatenate([event.first(selection=x, return_hepmc=False) for x in sel_truth],axis=0)

        # Get the final-state particles, as a numpy array.
        arr = event.all(selection=sel1, return_hepmc=False)

        # Do some event filtering using hadronization-level jets.
        # Note that if doing detector-level simulation (i.e. Delphes),
        # an event that has a jet passing cuts here is not guaranteed to
        # have a detector-level jet that also passes cuts, but it is likely.

        # ----- JET CLUSTERING -----
        # Perform jet clustering on final-state particles.
        pseudojets = []
        N = len(arr)
        for j in range(N):
            # Create a pseudojet object and add it to the list.
            # Signature is px, py, pz, E.
            vec = [arr[j][x] for x in ['px','py','pz','E']]
            pj = fj.PseudoJet(*vec)
            pjinfo = ParticleInfo(j, arr[j]['status'], arr[j]['pdgid'])
            pj.set_python_info(pjinfo)
            pseudojets.append(pj)

        jets = jetdef(pseudojets)
        # Apply optional minimum jet pT cut.
        jet_pt = np.array([jet.pt() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[jet_pt >= jet_config['jet_min_pt']]
        jets = [jets[i] for i in jet_indices]

        # Apply optional maximum |eta| cut.
        jet_eta = np.array([jet.eta() for jet in jets])
        jet_indices = np.linspace(0,len(jets)-1,len(jets),dtype=np.dtype('i8'))[np.abs(jet_eta) <= jet_config['jet_max_eta']]
        jets = [jets[i] for i in jet_indices]

        if(len(jets) == 0): # No jets passing cuts -> We can already toss out this event.
                n_fail += 1
                success = False
                break

#         # Now select a particular jet from this event -- we will only keep its constituents.
#         selected_jet_idx = jet_config['jet_selection'](truth = arr_truth, jets = jets)
#         jet = jets[selected_jet_idx]

#         # Get the constituents of our selected jet.
#         jet_constituents = np.array([
#             [x.px(), x.py(), x.pz(), x.E(), x.python_info().pdg_id, x.python_info().status]
#             for x in jet.constituents()
#         ])
#         # ----- END JET CLUSTERING -----
#         arr = jet_constituents

        if(use_delphes):
            # ----- BASIC FINAL-STATE PARTICLE SELECTION -----
            # Now in order to slim down events, we will only keep final-state particles within a certain
            # distance of the selected truth particles. We should be generous with our selection to avoid
            # throwing out things that have any chance of making their way into our jets later on.

            fields = ['px','py','pz','E','pdgid','status','eta','phi']
            arr_truth = np.array([[x[y] for y in fields] for x in arr_truth])
            arr = np.array([[x[y] for y in fields] for x in arr])
            dr_limit = alpha * jet_config['jet_radius']

            N = len(arr)
            min_distances = np.min(DeltaR2Vectorized(arr[:,-2:], arr_truth[:,-2:]),axis=1)
            selected = (np.abs(min_distances) < np.abs(dr_limit))
            arr = arr[selected]

        else:
            # ----- JET-BASED FINAL-STATE PARTICLE SELECTION -----
            # If we're not using Delphes, then the pre-detector-level jets we used above for filtering
            # are the very same jets that we're interested in, so we simply save the constituents of
            # the selected jet.
            selected_jet_idx = jet_config['jet_selection'](truth = arr_truth, jets = jets, max_dr = max_dr)
            if(selected_jet_idx < 0): # No jet passed selection.
                n_fail += 1
                success = False
                break

            jet = jets[selected_jet_idx]
            jet_constituents = np.array([
                [x.px(), x.py(), x.pz(), x.E(), x.python_info().pdg_id, x.python_info().status]
                for x in jet.constituents()]
                )
            arr = jet_constituents

        # Combine the particle arrays.
        arr = np.row_stack((arr_truth, arr))
        del arr_truth

        # Create a GenEvent (pyhepmc_ng) containing these selected particles.
        # Note that the event comes from pyhepmc_ng, *not* numpythia.
        # Thus we cannot just use GenParticles from numpythia (e.g. using return_hepmc=True above).
        N = len(arr)
        hepev = hep.GenEvent()
        hepev.event_number = i_real
        for j in range(N):
            momentum = arr[j][:4] # px, py, pz, E.
            pid = int(arr[j][4])
            status = int(arr[j][5])
            par = hep.GenParticle(momentum=momentum, pid=pid, status=status)
            hepev.add_particle(par)

        # ----- File I/O -----
        # Write this event to the HepMC buffer file.
        with hep.WriterAscii(buffername) as f:
            f.write_event(hepev)

        # Now copy the buffer file contents into the full HepMC file.
        # For speed, we do this using Unix commands (though it seems a bit hacky).
        upper_trim = 3
        lower_trim = 3 # 2 if using Unix head
        if(loop_number == 0 and i_real == 1): upper_trim = 0
        elif(i == nevents-1 and n_fail == 0): lower_trim = 0
        comm1 = 'tail -n +{} {}'.format(upper_trim, buffername).split(' ')

        proc = sub.Popen(comm1,stdout=sub.PIPE,text=True)
        if(lower_trim == 0): proc_out = proc.stdout.read().split('\n')
        else: proc_out = proc.stdout.read().split('\n')[:-lower_trim]

        with open(filename,'a') as f:
            f.writelines(line + '\n' for line in proc_out)
        # ----- End File I/O -----

        qu.printProgressBarColor(i_real,nevents_disp, prefix=prefix, suffix=suffix, length=bl)
        i_real += 1

    # delete the buffer file
    comm = ['rm', buffername]
    sub.check_call(comm)

    return i_real-1, n_fail

# Generate a bunch of events in the given pT range,
# and save them to a HepMC file.
# We do perform event selection: Only certain particles are saved to the file to begin with.
def Generate(nevents, pt_min, pt_max, filename = 'events.hepmc'): # TODO: implement file chunking
    pythia_config = GetPythiaConfig(pt_min = pt_min, pt_max = pt_max)
    pythia = npyth.Pythia(params = pythia_config)

    # Get our (Pythonic) Fastjet. # TODO: Make this optional (only need if *not* using Delphes)
    fastjet_dir = BuildFastjet(j=8)
    fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]

    # Get the Fastjet banner out of the way
    if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
    import fastjet as fj
    fj.ClusterSequence.print_banner()

    prefix = 'Generating events for pT bin [{},{}]:'.format(pt_min,pt_max)
    suffix = 'Complete'
    bl = 50
    qu.printProgressBarColor(0,nevents, prefix=prefix, suffix=suffix, length=bl)

    sel1 = (npyth.STATUS == 1) & (npyth.ABS_PDG_ID > -1) # Selection for all final-state particles.
    sel_truth = GetTruthSelection()

    # Loop in such a way as to guarantee that we get as many events as requested.
    # This logic is required as events could technically fail selections, e.g. not have the
    # requested truth particles (depends on requested truth particles & processes).
    n_success = 0
    n_fail = nevents
    nloops = 0
    while(n_fail > 0):
        n_success, n_fail = GenerationLoop(pythia, nevents-n_success, sel1, sel_truth, filename, prefix, suffix, bl=50,
                                           i_real=n_success+1, nevents_disp = nevents, fastjet_dir = fastjet_dir, loop_number = nloops)
        nloops = nloops + 1
    return

# Create a copy of a HepMC file containing only truth-level particles (i.e. remove final-state info).
def CopyTruth(hepfile, outfile=None):

    npars = []
    npar = 0

    if(outfile is None): outfile = hepfile.replace('.hepmc','_truth.hepmc')
    with open(hepfile,'r') as f, open(outfile,'w') as g:
        for line in f:
            if(':') in line: g.write(line)
            else:
                line = line.replace('\n','').split(' ')
                if(line[-1] != '1'):
                    if(line[0] == 'E'):
                        npars.append(npar)
                        npar = 0
                        g.write(' '.join(line) + '\n')
                    elif(line[0] == 'P'):
                        npar +=1
                        g.write(' '.join(line) + '\n')
    npars.append(npar)
    npars = npars[1:]

    # Now we need to correct the event headers, to have the right number of particles.
    # While it might not be the most efficient, this is easy to accomplish on a 2nd pass.
    # The file is also likely so small now that it's very easy to store in memory.
    event_counter = 0
    lines = []

    with open(outfile,'r') as f:
        for line in f:
            if(':' in line): lines.append(line)
            else:
                line = line.replace('\n','').split(' ')
                if(line[0] == 'E'):
                    line[3] = str(npars[event_counter])
                    event_counter += 1
                lines.append(' '.join(line) + '\n')

    with open(outfile,'w') as f:
        for line in lines:
            f.write(line)

    return outfile