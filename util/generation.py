import glob
import subprocess as sub
import numpythia as npyth # Pythia, hepmc_write
from util.config import GetEventSelection, GetFinalStateSelection, GetPythiaConfig, GetPythiaConfigFile, GetTruthSelection
from util.fastjet import BuildFastjet
from util.gen_utils.utils import HepMCOutput, CreateHepMCEvent, JetConstituentSelection, TruthDistanceSelection
from util.conv_utils.utils import InitFastJet
from util.particle_selection import SelectFinalState, SelectSimplestHadronic
import util.qol_util as qu

# TODO: This is a compatibility check for the pythia configuration.
# Numpythia v. 1.2 uses Pythia 8.244, but Numpythia v. 1.1 uses Pythia 8.226.
# There are some noticeable differences in available process configs, and the latest numpythia
# may not be available on macOS/arm yet.
# The check is not very efficient, it will make many nympythia objects.
def PythiaCompatibilityCheck(pythia_config):
    pythia_config_safe = {}
    safe_keys = ['Beam','Random','Stat','Print','PhaseSpace','PartonLevel','HadronLevel']
    for key,val in pythia_config.items():
        safe = False
        for safe_key in safe_keys:
            if(safe_key in key):
                safe = True
                break
        if(safe):
            pythia_config_safe[key] = val
            continue

        tmp_config = {key_safe:val_safe for key_safe,val_safe in pythia_config_safe.items()}
        tmp_config[key] = val

        with qu.stdout_redirected():
            try:
                pythia = npyth.Pythia(params=tmp_config)
                pythia_config_safe[key] = val
                del pythia # TODO: Is this helpful?
            except: pass
    return pythia_config_safe

# # Generate a bunch of events in the given pT range,
# # and save them to a HepMC file.
# # Note that we do not perform event selection: All generated events are saved.
# def GenerateSimple(nevents, pt_min, pt_max, filename = 'events.hepmc'):

#     pythia_config = GetPythiaConfig(pt_min = pt_min, pt_max = pt_max)
#     pythia_config = PythiaCompatibilityCheck(pythia_config)
#     pythia = npyth.Pythia(params = pythia_config)

#     prefix = 'Generating events for pT bin [{},{}]:'.format(pt_min,pt_max)
#     suffix = 'Complete'
#     i = 0
#     qu.printProgressBarColor(i,nevents, prefix=prefix, suffix=suffix, length=50)

#     for i,event in enumerate(npyth.hepmc_write(filename, pythia(events=nevents))):
#         qu.printProgressBarColor(i+1,nevents, prefix=prefix, suffix=suffix, length=50)
#     return

def GenerationLoop(pythia, nevents,
                   filename,
                   prefix, suffix, bl=50,
                   i_real = 1, nevents_disp=None,
                   loop_number=0, chunk_size=100,
                   max_dr = 0.2):

    jet_config, jetdef = InitFastJet()
    # if(fastjet_dir not in sys.path): sys.path.append(fastjet_dir)
    # import fastjet as fj
    # jet_config = GetJetConfig()
    # jetdef = fj.JetDefinition(fj.antikt_algorithm, jet_config['jet_radius'])

    n_fail = 0
    if(nevents_disp is None): nevents_disp = nevents # number of events to display in progress bar

    # The way that pyhepmc_ng's WriterAscii works, writing an event will overwrite the whole file.
    # Thus for the time being, we will circumvent this limitation by making a buffer file where each event
    # is written, and then copied to the "main" file before the next event is generated. This I/O might slow
    # down things, so we ultimately want to find some way to do a write with "append" functionality.

    filename_truth = filename.replace('.hepmc','_truth.hepmc')
    buffername = filename.replace('.hepmc','_buffer.hepmc')
    buffername_truth = buffername.replace('.hepmc','_truth.hepmc')

    for i,event in enumerate(pythia(events=nevents)):
        success = True

        # Get the truth-level particles.
        # TODO: In some cases of t -> b W, using GenEvent.all() yields two b-quarks (throwing off some later code).
        #       This happens without MPI/ISR/FSR/ and both of them have the same generator status. How?
        # arr_truth = np.concatenate([event.first(selection=x, return_hepmc=False) for x in sel_truth],axis=0)

        arr_truth = GetTruthSelection()(event=event)

        for truth in arr_truth:
            # print(truth)
            if(len(truth) == 0): # missing some truth particle -> potential trouble?
                n_fail += 1
                success = False
                break

        if(not success): continue

        # Get the final-state particles, as a numpy array.
        status, arr = GetFinalStateSelection()(event=event)

        # If we didn't pick up any final-state particles, discard this event.
        if(not status):
            n_fail += 1
            success = False
            continue

        # print('\n',len(arr),arr)

        # Optionally filter down events (e.g. throw out particles too far from selected truth particles).
        # This is typically some sort of filtering done just to reduce the HepMC file size (which can be rather large).
        event_selection = GetEventSelection()
        if(event_selection is not None):
            status, arr, arr_truth = event_selection(arr,arr_truth)

        # print()

        if(not status):
            n_fail += 1
            success = False
            continue

        # Create a GenEvent (pyhepmc_ng) containing these selected particles.
        # Note that the event comes from pyhepmc_ng, *not* numpythia.
        # Thus we cannot just use GenParticles from numpythia (e.g. using return_hepmc=True above).
        hepev = CreateHepMCEvent(arr,i_real)
        hepev_truth = CreateHepMCEvent(arr_truth,i_real)

        # ----- File I/O -----
        # Write this event to the HepMC buffer file, then copy the buffer contents to the full HepMC file.
        HepMCOutput(hepev,buffername,filename,loop_number,i,i_real,nevents,n_fail)
        HepMCOutput(hepev_truth,buffername_truth,filename_truth,loop_number,i,i_real,nevents,n_fail)
        # ----- End File I/O -----

        qu.printProgressBarColor(i_real,nevents_disp, prefix=prefix, suffix=suffix, length=bl)
        i_real += 1

    # Delete the buffer files.
    for fname in [buffername,buffername_truth]:
        comm = ['rm', fname]
        try: sub.check_call(comm,stderr=sub.DEVNULL)
        except: pass

    return i_real-1, n_fail

# Generate a bunch of events in the given pT range,
# and save them to a HepMC file.
# We do perform event selection: Only certain particles are saved to the file to begin with.
def Generate(nevents, pt_min, pt_max, filename = 'events.hepmc'): # TODO: implement file chunking
    pythia_config = GetPythiaConfig(pt_min = pt_min, pt_max = pt_max)
    pythia_config = PythiaCompatibilityCheck(pythia_config)
    pythia_config_file = GetPythiaConfigFile()
    pythia = npyth.Pythia(config=pythia_config_file,params=pythia_config)

    # Get our (Pythonic) Fastjet. # TODO: Make this optional (only need if *not* using Delphes)
    fastjet_dir = BuildFastjet(j=8)
    fastjet_dir = glob.glob('{}/**/site-packages'.format(fastjet_dir),recursive=True)[0]

    # Get the Fastjet banner out of the way
    tmp = InitFastJet()
    del tmp

    prefix = 'Generating events for pT bin [{},{}]:'.format(pt_min,pt_max)
    suffix = 'Complete'
    bl = 50
    qu.printProgressBarColor(0,nevents, prefix=prefix, suffix=suffix, length=bl)

    # TODO: Final-state selection is chosen within GenerationLoop. Should we do the same for truth selection?
    # sel_truth = GetTruthSelection()

    # Loop in such a way as to guarantee that we get as many events as requested.
    # This logic is required as events could technically fail selections, e.g. not have the
    # requested truth particles (depends on requested truth particles & processes).
    n_success = 0
    n_fail = nevents
    nloops = 0
    while(n_fail > 0):
        n_success, n_fail = GenerationLoop(pythia, nevents-n_success, filename, prefix, suffix, bl=50,
                                           i_real=n_success+1, nevents_disp = nevents, loop_number = nloops)
        nloops = nloops + 1
    return