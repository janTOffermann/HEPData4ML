import sys,os,pathlib,time,datetime,re,shlex,uuid
import numpy as np, h5py as h5
import argparse as ap
import subprocess as sub
from util.generation import Generator
from util.delphes import DelphesWrapper
from util.conversion import Processor, RemoveFailedFromHDF5, SplitH5, AddEventIndices, ConcatenateH5, MergeH5, AddConstantValue, AddMetaDataWithReference
from util.hepmc import CompressHepMC
from util.config import Configurator,GetConfigFileContent
import config.config as config

def none_or_str(value): # see https://stackoverflow.com/a/48295546
    if value == 'None':
        return None
    return value

def get_git_revision_short_hash(): # see https://stackoverflow.com/a/21901260
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        result = sub.check_output(['git', 'rev-parse', '--short', 'HEAD'],cwd=cwd).decode('ascii').strip()
    except:
        result = 'NO_GIT_HASH'
    return result

# Convenience function for file naming
def float_to_str(value):
    value_str = str(value)
    value_str = value_str.replace('.',',')
    return re.sub(',0$','',value_str)

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n',            '--nevents',           type=int,          required=True,            help='Number of events per pt bin.')
    parser.add_argument('-p',            '--ptbins',            type=float,        default=None, nargs='+',  help='Transverse momentum bin edges, in GeV. (The \hat{p_T} variable for Pythia8 configuration.)')
    parser.add_argument('-p_s',          '--ptbins_string',     type=none_or_str,  default=None,             help='Transverse momentum bin edges, in GeV, as comma-separated string. An alternative to the -p argument.')
    parser.add_argument('-o',            '--outfile',           type=str,          default='events.h5',      help='Output HDF5 file name.')
    parser.add_argument('-O',            '--outdir',            type=none_or_str,  default=None,             help='Output directory.')
    parser.add_argument('-g',            '--generation',        type=int,          default=True,             help='Whether or not to do event generation.')
    parser.add_argument('-s',            '--sep_truth',         type=int,          default=True,             help='Whether or not to store truth-level particles in separate arrays.')
    parser.add_argument('-ns',           '--n_sep_truth',       type=int,          default=-1,               help='How many truth particles to save in separate arrays -- will save the first n as given by the truth selection.')
    parser.add_argument('-d',            '--diagnostic_plots',  type=int,          default=False,            help='Whether or not to make diagnostic plots.')
    parser.add_argument('-v',            '--verbose',           type=int,          default=0,                help='Verbosity.')
    parser.add_argument('-h5',           '--hdf5',              type=int,          default=1,                help='Whether or not to produce final HDF5 files. If false, stops after HepMC or Delphes/ROOT file production.')
    parser.add_argument('-f',            '--force',             type=int,          default=0,                help='Whether or not to force generation -- if true, will possibly overwrite existing HepMC files in output directory.')
    parser.add_argument('-c',            '--compress',          type=int,          default=0,                help='Whether or not to compress HepMC files.')
    parser.add_argument('-cd',           '--clean_delphes',     type=int,          default=0,                help='Whether or not to clean up DELPHES/ROOT files.')
    parser.add_argument('-rng',          '--rng',               type=int,          default=None,             help='Pythia RNG seed. Will override the one provided in the config file.')
    parser.add_argument('-npc',          '--nentries_per_chunk',type=int,          default=int(1e4),         help='Number of entries to process per chunk, for jet clustering & conversion to HDF5.')
    parser.add_argument('-pb',           '--progress_bar',      type=int,          default=1,                help='Whether or not to print progress bar during event generation')
    parser.add_argument('-sp',           '--split',             type=int,          default=1,                help='Whether or not to split HDF5 file into training/validation/testing files.')
    parser.add_argument('-tf',           '--train_fraction',    type=float,        default=0.7,              help='Fraction of events to place in the training file.')
    parser.add_argument('-vf',           '--val_fraction',      type=float,        default=0.2,              help='Fraction of events to place in the validation file.')
    parser.add_argument('-df',           '--delete_full',       type=int,          default=0,                help='Whether or not to delete the full HDF5 file after splitting into train/validation/testing files.')
    parser.add_argument('-co',           '--compression_opts',  type=int,          default=7,                help='Compression option for final HDF5 file (0-9). Higher value means more compression.')
    parser.add_argument('-separate_h5',  '--separate_h5',       type=int,          default=1,                help='Whether or not to make separate HDF5 files for each pT bin.')
    parser.add_argument('-pc',           '--pythia_config',     type=none_or_str,  default=None,             help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-index_offset', '--index_offset',      type=int,          default=0,                help='Offset for event_idx.')
    parser.add_argument('-debug',        '--debug',             type=int,          default=0,                help='If > 0, will record the full final-state (i.e. before jet clustering/selection) in a separate key.')
    args = vars(parser.parse_args())

    start_time = time.time()

    # Produce a random identifier for this dataset.
    unique_id = str(uuid.uuid4())
    unique_id_short = str(uuid.uuid4())[:5] # a second, shorter random string -- probably more convenient to use, at the risk of a higher collision rate

    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
    pt_bin_edges_string = args['ptbins_string']
    h5_file = args['outfile']
    outdir = args['outdir']
    do_generation =args['generation'] > 0 # TODO: Find nicer way to handle Boolean -- argparse is weird.
    separate_truth_particles = args['sep_truth'] > 0
    n_separate_truth_particles = args['n_sep_truth']
    diagnostic_plots = args['diagnostic_plots'] > 0
    verbose = args['verbose'] > 0
    do_h5 = args['hdf5']
    compress_hepmc = args['compress']
    delete_delphes = args['clean_delphes']
    force = args['force']
    pythia_rng = args['rng']
    nentries_per_chunk = args['nentries_per_chunk']
    progress_bar = args['progress_bar']
    compression_opts = args['compression_opts']
    separate_h5 = args['separate_h5'] > 0
    pythia_config = args['pythia_config']
    debug = args['debug'] > 0
    index_offset = args['index_offset']

    if(pt_bin_edges is not None and pt_bin_edges_string is not None):
        print('Warning: Both "ptbins" and "ptbins_string" arguments were given. Using the "ptbins" argument.')

    if(pt_bin_edges is None and pt_bin_edges_string is None):
        print('Error: pT bin edges not specified.')
        assert(False)

    if(pt_bin_edges is None and pt_bin_edges_string is not None):
        pt_bin_edges = pt_bin_edges_string.split(',')
        pt_bin_edges = [float(x) for x in pt_bin_edges]

    nbins = len(pt_bin_edges) - 1

    split_files = args['split'] > 0
    train_frac = args['train_fraction']
    val_frac = args['val_fraction']
    test_frac = 1. - train_frac - val_frac
    delete_full = args['delete_full'] > 0
    if(not split_files): delete_full = False # otherwise we are throwing out all the final files

    if(test_frac < 0. and split_files):
        print('Error: Requested training fraction and validation fraction sum to more than 1, this leaves no events for the test file.')
        assert(False)

    # Configurator class, used for fetching information from our config file.
    config_dictionary = config.config
    configurator = Configurator(config_dictionary=config_dictionary)

    # Set up FastJet -- we will need this later on (except for the special use case of no jet clustering!).
    # To keep our printouts clean, we are initializing FastJet here instead of later on in a loop, so that
    # we can get the FastJet banner printout out of the way. We remove the banner with some ANSI printing
    # hackery, since it's really not useful and clutters up our printout (we acknowledge the use in the
    # documentation, having this unavoidable printout for a single package's use is quite gratuitous).
    print()
    dummy_processor = Processor(configurator)
    line_up = '\033[1A'
    line_clear = '\x1b[2K'
    for i in range(13):
        print(line_up, end=line_clear)

    # Setting the verbosity for the HDF5 conversion.
    # If there are many events it might take a bit, so some printout
    # is helpful to monitor the progress.
    h5_conversion_verbosity = 0
    if(nevents_per_bin >= 100): h5_conversion_verbosity = 1
    elif(nevents_per_bin >= 10000): h5_conversion_verbosity = 2

    use_delphes = configurator.GetDelphesConfig()

    # Keep track of some files we create.
    jet_files = []
    truth_files = []
    hist_files = []
    stat_files = []
    final_state_truth_overlap_indices_files = [] # only used in certain cases, these files record indices of particles that show up in both truth and final-state selections (indices w.r.t. final-state HepMC files)
    filter_flag_files = [] # only used in certain cases, these files hold a boolean flag given by the "event_filter_flag"

    # Prepare the output directory.
    if(outdir is None): outdir = os.getcwd()
    else: os.makedirs(outdir,exist_ok=True)

    # Create a log file in the output directory, detailing the command line arguments.
    logfile = outdir + '/command_options.txt'
    with open(logfile,'w') as f:
        for key,val in args.items():
            if(val == True): val = 1
            elif(val == False): val = 0
            f.write('{} = {}\n'.format(key,val))

    # Also copy the configuration file, it is currently always config/config.py.
    this_dir = os.path.dirname(os.path.abspath(__file__))
    comm = ['cp','{}/config/config.py'.format(this_dir),'{}/config.py'.format(outdir)]
    sub.check_call(comm)

    if(do_generation):
        if(verbose):
            print('\n=================================')
            print('Running Pythia8 event generation.')
            print('=================================\n')

        if(pythia_rng is not None):
            print('\tSetting Pythia RNG seed to {}. (overriding config)'.format(pythia_rng))
        else:
            pythia_rng = configurator.GetPythiaRNGSeed()
        if(pythia_config is not None):
            print('\tSetting Pythia process configuration from {}. (overriding config)'.format(pythia_config))

        if(verbose):
            print('\tGenerating {} events per {} bin, with the following bin edges (in GeV):'.format(nevents_per_bin,'\\hat{p_T}'))
            for bin_edge in pt_bin_edges:
                print('\t\t{}'.format(bin_edge))
            print()

    print()
    for i in range(nbins):
        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.
        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]
        hep_file = 'events_{}-{}.hepmc'.format(float_to_str(pt_min),float_to_str(pt_max))

        generator = Generator(pt_min,pt_max, configurator, pythia_rng,pythia_config_file=pythia_config)
        generator.SetEventSelection(configurator.GetEventSelection())
        generator.SetTruthSelection(configurator.GetTruthSelection())
        generator.SetFinalStateSelection(configurator.GetFinalStateSelection())
        generator.SetEventFilter(configurator.GetEventFilter())
        generator.SetEventFilterFlag(configurator.GetEventFilterFlag())
        generator.SetJetConfig(configurator.GetJetConfig())
        generator.SetNTruth(configurator.GetNPars()['n_truth'])

        generator.SetOutputDirectory(outdir)

        # hist_filename = 'hists_{}.root'.format(i)
        hist_filename = hep_file.replace('.hepmc','_hists.root')
        generator.SetHistFilename(hist_filename)

        generator.SetFilename(hep_file, rename_extra_files=False)

        # We will use one file to hold additional event data -- cross-sections,
        # information on any "event filter flags", and on any overlapping listing
        # of particles between our truth- and final-state selections. These can
        # be handled by separate files but combining them will limit clutter.
        extra_data_file = hep_file.replace('.hepmc','_data.h5')
        generator.SetStatsFilename(extra_data_file)
        generator.SetEventFilterFlagFilename(extra_data_file)
        generator.SetIndexOverlapFilename(extra_data_file)
        generator.SetDiagnosticPlots(diagnostic_plots)
        generator.SetProgressBar(progress_bar)

        hepfile_exists = pathlib.Path('{}/{}'.format(outdir,hep_file)).exists()
        generate = True
        if(not do_generation):
            generate = False
        elif(hepfile_exists and not force):
            print('\tHepMC3 file {}/{} already found, skipping its generation.'.format(outdir,hep_file))
            generate = False

        if(generate): generator.Generate(nevents_per_bin)

        # Extract the truth-level particles from the full HepMC file.
        truthfile = generator.GetTruthFilename()
        truth_files.append(truthfile)

        stat_files.append(generator.GetStatsFilename())
        print('Stats filename = {}'.format(generator.GetStatsFilename()))
        final_state_truth_overlap_indices_files.append(generator.GetIndexOverlapFilename())
        filter_flag_files.append(generator.GetEventFilterFlagFilename())
        hist_files.append(generator.GetHistFilename())

        if(use_delphes): # Case 1: Using Delphes
            # Pass the HepMC file to Delphes. Will output a ROOT file.
            delphes_card = configurator.GetDelphesCard() # will default to the ATLAS card that is shipped with Delphes
            delphes_file = hep_file.replace('.hepmc','.root')
            delphes_wrapper = DelphesWrapper(configurator.GetDelphesDirectory())
            delphes_wrapper.PrepDelphes() # will download/build Delphes if necessary
            print('Running DelphesHepMC3: {} -> {}.'.format(hep_file,delphes_file))
            if(i == 0):
                print('\tDelphes executable: {}'.format(delphes_wrapper.GetExecutable()))

            delphes_file = delphes_wrapper.HepMC3ToDelphes(hepmc_file=hep_file, output_file=delphes_file, cwd=outdir, delphes_card=delphes_card)
            jet_files.append(delphes_file)

            # We can now compress the HepMC file. Once it's compressed, we delete the original file to recover space.
            if(compress_hepmc): CompressHepMC(['{}/{}'.format(outdir,hep_file)],True)

        else: # Case 2: No Delphes
            jet_files.append(hep_file)

    # Join together the histogram ROOT files.
    hist_filename = 'hists.root'
    if(do_generation and diagnostic_plots):
        comm = ['hadd', hist_filename]
        comm += hist_files
        sub.check_call(comm,cwd=outdir,stderr=sub.DEVNULL,stdout=sub.DEVNULL)
        comm = ['rm'] + hist_files
        sub.check_call(comm,cwd=outdir)

    # Optionally stop here, if we're not interested in progressing beyond HepMC or Delphes/ROOT.
    if(not do_h5):
        return

    # Now put everything into an HDF5 file.
    # We can optionally package things into separate HDF5 files, one per pT-hat bin, and then concatenate those later.
    # This could be useful for certain use cases.
    if(verbose): print('\nRunning jet clustering and producing final HDF5 output.\n')
    processor = Processor(configurator,use_delphes)
    processor.SetOutputDirectory(outdir)
    print("processor.outdir = {}".format(outdir))
    processor.SetHistFilename(hist_filename)
    processor.SetDiagnosticPlots(diagnostic_plots)
    processor.SetSeparateTruthParticles(separate_truth_particles)
    processor.SetNSeparateTruthParticles(n_separate_truth_particles)
    processor.SetRecordFullFinalState(debug)

    if(not separate_h5):
        # The simple way to run this code -- it immediately makes a single HDF5 file with all events, from all pT-binned HepMC/Delphes files,
        # and this file will optionally be split (in a later step) into training, testing and validation samples.
        print('\tProducing a single HDF5 file, from all the HepMC or Delphes/ROOT files.')
        print('\tWarning: This will skip any requested post-processing steps. If you want these, re-run with "-separate_h5 1".')

        processor.Process(jet_files,truth_files,h5_file,verbosity=h5_conversion_verbosity,nentries_per_chunk=nentries_per_chunk)
        processor.MergeEventFilterFlag(h5_file,filter_flag_files,copts=compression_opts)
    else:
        # Slightly more complex way of running things -- we first create pT-binned HDF5 files, corresponding with the binning of the
        # HepMC/Delphes files. These binned files are then concatenatd into a single HDF5 file with all events.
        # Making the binned files first may allow us to (more easily) run some algorithms where we may want to compare across
        # the HepMC/Delphes and HDF5 files, since they are all binned it'll be easy to match events across files.
        h5_files = []
        print('\tProducing separate HDF5 files for each pT bin, and then concatenating these.')
        delete_individual_h5 = True
        nentries_per_chunk = int(nentries_per_chunk/nbins)
        for i in range(nbins):
            pt_min = pt_bin_edges[i]
            pt_max = pt_bin_edges[i+1]
            jet_file = [jet_files[i]]
            truth_file = [truth_files[i]]
            h5_file_individual = h5_file.replace('.h5','_{}-{}.h5'.format(float_to_str(pt_min),float_to_str(pt_max)))
            processor.SetProgressBarPrefix('\n\tClustering jets & preparing data for pT bin [{},{}]:'.format(pt_min,pt_max))

            processor.Process(jet_file,truth_file,h5_file_individual,verbosity=h5_conversion_verbosity,nentries_per_chunk=nentries_per_chunk)

            # To each HDF5 event file, we will add event indices. These may be useful/necessary for the post-processing step.
            # We will remove these indices when concatenating files, since as one of our last steps we'll add indices again
            # but with respect to the full event listing (not just w.r.t. events in each generation pT bin).
            AddEventIndices(h5_file_individual,cwd=outdir,copts=compression_opts)

            # Optional post-processing. Any post-processing steps have been configured in the config file, config/config.py.
            processor.PostProcess(jet_file,[h5_file_individual],[final_state_truth_overlap_indices_files[i]])

            # # Optionally add any "event_filter_flags" that were set in the configuration. If none were set, this doesn't do anything.
            # processor.MergeEventFilterFlag(h5_file_individual,filter_flag_files[i],copts=compression_opts)

            h5_file_individual = '/'.join((outdir,h5_file_individual))
            h5_files.append(h5_file_individual)

        print('\n\tConcatenating HDF5 files. Will drop the "event_idx" key, \n\tthis was used internally for any post-processing steps.\n\tIndices will be recomputed and added at the end.')
        ConcatenateH5(h5_files,'/'.join((outdir,h5_file)),copts=compression_opts,delete_inputs=delete_individual_h5,ignore_keys=['event_idx'],verbose=True,silent_drop=True)

    if(use_delphes):
        if(delete_delphes):
            # Cleanup: Delete the jet files -- which are Delphes/ROOT files, since they can always be recreated from the (compressed) HepMC files.
            jet_files = ['{}/{}'.format(outdir,x) for x in jet_files]
            comm = ['rm'] + jet_files
            sub.check_call(comm)

    else:
        #Cleanup: Compress the HepMC files.
        if(compress_hepmc): CompressHepMC(jet_files,True,cwd=outdir)

    # Combine the stats files.
    # TODO: Can we handle some of the stats file stuff under-the-hood? Or just access all the files
    # without making the aggregate stats file.
    stats_filename = 'stats.h5'
    delete_individual_stats = True
    delete_full_stats = True # file is now just used for merging purposes, it is temporary
    try: ConcatenateH5(stat_files,stats_filename,cwd=outdir,delete_inputs=delete_individual_stats, copts=9)
    except: pass # as long as the full stats file exists, it's okay if the individual ones were deleted already

    try: MergeH5(h5_file,stats_filename,cwd=outdir,delete_stats_file=delete_full_stats, copts=9)
    except: print('Warning: Stats information not found!')

    # Now also compress the truth files.
    if(compress_hepmc): CompressHepMC(truth_files,True,cwd=outdir)

    # Add some event indices to our dataset.
    if(index_offset < 0): index_offset = 0
    print('\tAdding event indices to file {}.'.format('/'.join((outdir,h5_file))))
    print('\t\tStarting at event_idx = {}.'.format(index_offset))
    AddEventIndices(h5_file,cwd=outdir,copts=compression_opts,offset=index_offset)

    # Remove any failed events (e.g. detector-level events with no jets passing cuts).
    print('\tRemoving any failed events from file {}.\n\tThese may be events where there weren\'t any jets passing the requested cuts.'.format('/'.join((outdir,h5_file))))
    RemoveFailedFromHDF5(h5_file,cwd=outdir)

    # Now, add some metadata to the file -- we use the HDF5 file attributes to store lists of metadata, and create columns that reference these lists.
    # This is handled correctly by metadata.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=pythia_rng,                                                 key='pythia_random_seed'    ) # Add the Pythia8 RNG seed. Storing this way doesn't really save space -- it's just an int -- but we'll do this for consistency with how metadata is handled.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=" ".join(map(shlex.quote, sys.argv[1:])),                   key='command_line_arguments') # Add the command line arguments.
    AddMetaDataWithReference(h5_file,cwd=outdir,value='\n'.join(GetConfigFileContent()),                          key='config_file'           ) # Add the full config file a string.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=start_time,                                                 key='timestamp'             ) # Add the epoch time for the start of generation.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(start_time)), key='timestamp_string_utc'  ) # Add the epoch time for the start of generation, as a string.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=unique_id,                                                  key='unique_id'             ) # A unique random string to identify this dataset.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=unique_id_short,                                            key='unique_id_short'       ) # A unique random string to identify this dataset. (A shorter one, at the risk of increased likelihood of collisions)
    AddMetaDataWithReference(h5_file,cwd=outdir,value=get_git_revision_short_hash(),                              key='git_hash'              ) # Git hash for the data generation code.
    AddMetaDataWithReference(h5_file,cwd=outdir,value=configurator.GetPythiaConfigFileContents(pythia_config),    key='pythia_config'         ) # Pythia configuration (except for "\hat{p_T}", which is handled externally)

    # TODO: Might want to think about offering the ability to split the HepMC3 and Delphes files too?
    #       Could be useful for certain post-processing where we need to access those files, and we
    #       want to work on the split files (insead of the full events.h5).
    #       Before splitting those, we'd want to effectively join them together first, they are currently
    #       split by pT bin.
    if(split_files):
        # Now split the HDF5 file into training, testing and validation samples.
        split_ratio = (train_frac,val_frac,test_frac)
        print("\tSplitting HDF5 file {} into training, validation and testing samples:".format('/'.join((outdir,h5_file))))
        train_name = 'train.h5'
        val_name = 'valid.h5'
        test_name = 'test.h5'
        SplitH5(h5_file, split_ratio,cwd=outdir,copts=compression_opts, train_name=train_name,val_name=val_name,test_name=test_name,verbose=True,seed=configurator.GetSplitSeed())

    # Optionally delete the full HDF5 file.
    if(delete_full):
        comm = ['rm','{}/{}'.format(outdir,h5_file)]
        sub.check_call(comm)

    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_time_readable = str(datetime.timedelta(seconds=elapsed_time))
    print('\n#############################')
    print('Done. Time elapsed = {:.1f} seconds.'.format(elapsed_time))
    print('({})'.format(elapsed_time_readable))
    print('#############################\n')

if __name__ == '__main__':
    main(sys.argv)