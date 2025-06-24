import sys,os,pathlib,time,datetime,re,shlex,uuid,atexit
import argparse as ap
import subprocess as sub
from util.generation import PythiaGenerator
from util.simulation import DelphesSimulator
from util.conversion import Processor, RemoveFailedFromHDF5, SplitH5, AddEventIndices, ConcatenateH5, MergeH5, AddConstantValue, AddMetaDataWithReference
from util.hepmc.hepmc import CompressHepMC
# from util.hepmc.setup import HepMCSetup
from util.config import Configurator,GetConfigFileContent, GetConfigDictionary
from util.args import parse_mc_steps, FloatListAction, none_or_str




def trace_hepmc3_imports():
    """Add import tracing to see what imports pyHepMC3 first"""
    original_import = __builtins__.__import__

    def traced_import(name, *args, **kwargs):
        if 'pyHepMC3' in name or 'HepMC3' in name:
            import traceback
            print(f"\n=== IMPORTING {name} ===")
            traceback.print_stack()
            print("=" * 40)
        return original_import(name, *args, **kwargs)

    __builtins__.__import__ = traced_import

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

    # trace_hepmc3_imports()

    parser = ap.ArgumentParser()

    parser.add_argument('-n',            '--nevents',           type=int,          required=True,            help='Number of events per pt bin.')
    parser.add_argument('-steps',        '--steps',             type = parse_mc_steps, default = 'generation simulation reconstruction', help='Comma- or space-separated list of step. Options are [generation,simulation,reconstruction].')
    parser.add_argument('-p',            '--ptbins',            action = FloatListAction, default=[-1,-1], nargs='*',    help='Transverse momentum bin edges, for outgoing particles of the hard process. Can be a list of floats, or a string of comma- or space-separated floats. In GeV.')
    parser.add_argument('-o',            '--outfile',           type=str,          default='events.h5',      help='Output HDF5 file name.')
    parser.add_argument('-O',            '--outdir',            type=none_or_str,  default=None,             help='Output directory.')
    parser.add_argument('-v',            '--verbose',           type=int,          default=0,                help='Verbosity.')
    parser.add_argument('-f',            '--force',             type=int,          default=0,                help='Whether or not to force generation -- if true, will possibly overwrite existing HepMC files in output directory.')
    parser.add_argument('-c',            '--compress',          type=int,          default=0,                help='Whether or not to compress HepMC files.')
    parser.add_argument('-rng',          '--rng',               type=int,          default=None,             help='Pythia RNG seed. Will override the one provided in the config file.')
    parser.add_argument('-npc',          '--nentries_per_chunk',type=int,          default=int(1e4),         help='Number of entries to process per chunk, for jet clustering & conversion to HDF5.')
    parser.add_argument('-pb',           '--progress_bar',      type=int,          default=1,                help='Whether or not to print progress bar during event generation')
    parser.add_argument('-sp',           '--split',             type=int,          default=1,                help='Whether or not to split HDF5 file into training/validation/testing files.')
    parser.add_argument('-tf',           '--train_fraction',    type=float,        default=0.7,              help='Fraction of events to place in the training file.')
    parser.add_argument('-vf',           '--val_fraction',      type=float,        default=0.2,              help='Fraction of events to place in the validation file.')
    parser.add_argument('-df',           '--delete_full',       type=int,          default=0,                help='Whether or not to delete the full HDF5 file after splitting into train/validation/testing files.')
    parser.add_argument('-ds',           '--delete_stats',      type=int,          default=1,                help='Whether or not to delete the full stats file. This file\'s info is merged into HDF5 dataset, but the file may be useful in some advanced use cases.')
    parser.add_argument('-co',           '--compression_opts',  type=int,          default=7,                help='Compression option for final HDF5 file (0-9). Higher value means more compression.')
    parser.add_argument('-pc',           '--pythia_config',     type=none_or_str,  default=None,             help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-index_offset', '--index_offset',      type=int,          default=0,                help='Offset for event_idx.')
    parser.add_argument('-config',       '--config',            type=str,          default=None,             help='Path to configuration Python file. Default will use config/config.py .')

    # DELPHES-related arguments
    parser.add_argument('-delphes',      '--delphes',           type=int,          default=-1,               help='Optional Delphes flag -- will override the \'delphes\' option in the config file. 0 == False, 1 == True, otherwise ignored.')
    parser.add_argument('-del_delphes',  '--del_delphes',       type=int,          default=0,                help='Whether or not to delete DELPHES/ROOT files.')

    args = vars(parser.parse_args())

    start_time = time.time()

    # Produce a random identifier for this dataset.
    unique_id = str(uuid.uuid4())
    unique_id_short = str(uuid.uuid4())[:5] # a second, shorter random string -- probably more convenient to use, at the risk of a higher (but still tiny) collision rate

    steps = args['steps']
    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
    h5_file = args['outfile']
    outdir = args['outdir']
    verbose = args['verbose'] > 0
    compress_hepmc = args['compress']
    force = args['force']
    pythia_rng = args['rng']
    nentries_per_chunk = args['nentries_per_chunk']
    progress_bar = args['progress_bar']
    compression_opts = args['compression_opts']
    pythia_config = args['pythia_config']
    index_offset = args['index_offset']
    delete_full_stats = args['delete_stats'] > 0
    config_file = args['config']

    delphes_override = args['delphes']
    delete_delphes = args['del_delphes']

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
    # We import this from a user-supplied file, by default it is config/config.py.
    this_dir = os.path.dirname(os.path.abspath(__file__))
    if(config_file is None):
        config_file = '{}/config/config.py'.format(this_dir)

    print('Using configuration file: {} .'.format(config_file))
    config_dictionary = GetConfigDictionary(config_file)
    configurator = Configurator(config_dictionary=config_dictionary)

    # Try to correct the configuration.
    if(not configurator.GetStatus()):
        print('Attempting to correct filepaths. This may not work.')
        configurator.CorrectFilepaths(this_dir)

    if(not configurator.GetStatus()):
        print('Error: Configuration has bad status. Exiting.')
        assert(False)

    # # Set up FastJet -- we will need this later on (except for the special use case of no jet clustering!).
    # # To keep our printouts clean, we are initializing FastJet here instead of later on in a loop, so that
    # # we can get the FastJet banner printout out of the way. We remove the banner with some ANSI printing
    # # hackery, since it's really not useful and clutters up our printout (we acknowledge the use in the
    # # documentation, having this unavoidable printout for a single package's use is quite gratuitous).
    # print(13 * '\n')
    # dummy_processor = Processor(configurator)
    # line_up = '\033[1A'
    # line_clear = '\x1b[2K'
    # for i in range(13):
    #     print(line_up, end=line_clear)

    # Setting the verbosity for the HDF5 conversion.
    # If there are many events it might take a bit, so some printout
    # is helpful to monitor the progress.
    h5_conversion_verbosity = 0
    if(nevents_per_bin >= 100): h5_conversion_verbosity = 1
    elif(nevents_per_bin >= 10000): h5_conversion_verbosity = 2

    if(delphes_override == 0):
        configurator.SetDelphesConfig(False)
    elif(delphes_override == 1):
        configurator.SetDelphesConfig(True)
    simulation_type = configurator.GetSimulationType()

    # Keep track of some files we create.
    hepmc_files = []
    delphes_files = []
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

    #=========================
    # STEP 1: Generation
    #=========================

    if('generation' in steps):
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

    else:
        pythia_rng = -1 # TODO: Make this better -- currently this means seed is unknown!

    print()
    for i in range(nbins):
        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.
        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]

        #TODO: Toggle ROOT vs. ASCII
        hepmc_extension = 'hepmc'
        if(configurator.GetHepMCFormat().lower() == 'root'):
            hepmc_extension = 'root'

        hep_file = 'events_{}-{}.{}'.format(float_to_str(pt_min),float_to_str(pt_max),hepmc_extension)

        generator = PythiaGenerator(pt_min,pt_max, configurator, pythia_rng,pythia_config_file=pythia_config)
        # generator.SetEventFilter(configurator.GetEventFilter())
        # generator.SetEventFilterFlag(configurator.GetEventFilterFlag())

        generator.SetOutputDirectory(outdir)
        generator.SetFilename(hep_file, rename_extra_files=False)
        generator.SetProgressBar(progress_bar)

        hepfile_exists = pathlib.Path('{}/{}'.format(outdir,hep_file)).exists()
        generate = True
        if(not 'generation' in steps):
            generate = False
        elif(hepfile_exists and not force):
            print('\tHepMC3 file {}/{} already found, skipping its generation.'.format(outdir,hep_file))
            generate = False

        if(generate): generator.Generate(nevents_per_bin)

        filter_flag_files.append(generator.GetEventFilterFlagFilename())
        hepmc_files.append(hep_file)

    #===============================
    # STEP 2: Simulation (optional)
    #===============================

    if('simulation' in steps):
        simulator = None
        if(simulation_type == 'delphes'):
            sim_logfile = '{}/delphes.log'.format(outdir)
            simulator = DelphesSimulator(configurator,outdir,logfile=sim_logfile)

        if(simulator is not None):
            simulator.SetInputs(hepmc_files)
            simulator.Process()
            delphes_files = simulator.GetOutputFiles()

    #======================================================
    # STEP 3: Conversion + Reconstruction + Post-processing
    #======================================================
    if('reconstruction' in steps):
        # Now put everything into an HDF5 file.
        # We can optionally package things into separate HDF5 files, one per pT-hat bin, and then concatenate those later.
        # This could be useful for certain use cases.
        if(verbose): print('\nRunning jet clustering and producing final HDF5 output.\n')
        processor = Processor(configurator)
        processor.SetNentriesPerChunk(100) # the larger this is, the larger the chunks in memory (and higher the memory usage)
        processor.SetDelphesFiles(delphes_files)
        processor.SetOutputDirectory(outdir)

        h5_files = []
        print('\nProducing separate HDF5 files for each pT bin, and then concatenating these.')
        delete_individual_h5 = True
        nentries_per_chunk = int(nentries_per_chunk/nbins)
        for i in range(nbins):
            # TODO: Rework this a little. Should just generically loop over HepMC files, since they might have an external source and not be pt-binned.
            pt_min = pt_bin_edges[i]
            pt_max = pt_bin_edges[i+1]
            hepmc_files = [hepmc_files[i]]
            h5_file_individual = h5_file.replace('.h5','_{}-{}.h5'.format(float_to_str(pt_min),float_to_str(pt_max)))
            processor.SetProgressBarPrefix('\tConverting HepMC3 -> HDF5 for pT bin [{},{}]:'.format(pt_min,pt_max))

            processor.Process(hepmc_files,h5_file_individual,verbosity=h5_conversion_verbosity)

            # To each HDF5 event file, we will add event indices. These may be useful/necessary for the post-processing step.
            # We will remove these indices when concatenating files, since as one of our last steps we'll add indices again
            # but with respect to the full event listing (not just w.r.t. events in each generation pT bin).
            AddEventIndices(h5_file_individual,cwd=outdir,copts=compression_opts)

            # Optional post-processing. Any post-processing steps have been configured in the config file, config/config.py.
            processor.PostProcess(hepmc_files,[h5_file_individual])

            # # Optionally add any "event_filter_flags" that were set in the configuration. If none were set, this doesn't do anything.
            # processor.MergeEventFilterFlag(h5_file_individual,filter_flag_files[i],copts=compression_opts)

            h5_file_individual = '/'.join((outdir,h5_file_individual))
            h5_files.append(h5_file_individual)

        print('\n\tConcatenating HDF5 files. Will drop the "event_idx" key, \n\tthis was used internally for any post-processing steps.\n\tIndices will be recomputed and added at the end.')
        ConcatenateH5(h5_files,'/'.join((outdir,h5_file)),copts=compression_opts,delete_inputs=delete_individual_h5,ignore_keys=['event_idx'],verbose=True,silent_drop=True)

        if(simulation_type == 'delphes'):
            if(delete_delphes):
                # Cleanup: Delete the jet files -- which are Delphes/ROOT files, since they can always be recreated from the (compressed) HepMC files.
                delphes_files = ['{}/{}'.format(outdir,x) for x in delphes_files]
                comm = ['rm'] + delphes_files
                sub.check_call(comm)

        else:
            #Cleanup: Compress the HepMC files.
            if(compress_hepmc): CompressHepMC(hepmc_files,True,cwd=outdir)

        # Add some event indices to our dataset.
        if(index_offset < 0): index_offset = 0
        print('\tAdding event indices to file {}.'.format('/'.join((outdir,h5_file))))
        print('\t\tStarting at event_idx = {}.'.format(index_offset))
        AddEventIndices(h5_file,cwd=outdir,copts=compression_opts,offset=index_offset)

        # Remove any failed events (e.g. detector-level events with no jets passing cuts).
        # print('\tRemoving any failed events from file {}.\n\tThese may be events where there weren\'t any jets passing the requested cuts.'.format('/'.join((outdir,h5_file))))
        print('\tRemoving any failed events from file {}.'.format('/'.join((outdir,h5_file))))
        RemoveFailedFromHDF5(h5_file,cwd=outdir)

        # Now, add some metadata to the file -- we use the HDF5 file attributes to store lists of metadata, and create columns that reference these lists.
        # This is handled correctly by metadata.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=pythia_rng,                                                 key='pythia_random_seed'    ) # Add the Pythia8 RNG seed. Storing this way doesn't really save space -- it's just an int -- but we'll do this for consistency with how metadata is handled.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=" ".join(map(shlex.quote, sys.argv[1:])),                   key='command_line_arguments') # Add the command line arguments.
        AddMetaDataWithReference(h5_file,cwd=outdir,value='\n'.join(GetConfigFileContent(config_file)),               key='config_file'           ) # Add the full config file a string.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=start_time,                                                 key='timestamp'             ) # Add the epoch time for the start of generation.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(start_time)), key='timestamp_string_utc'  ) # Add the epoch time for the start of generation, as a string.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=unique_id,                                                  key='unique_id'             ) # A unique random string to identify this dataset.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=unique_id_short,                                            key='unique_id_short'       ) # A unique random string to identify this dataset. (A shorter one, at the risk of increased likelihood of collisions)
        AddMetaDataWithReference(h5_file,cwd=outdir,value=get_git_revision_short_hash(),                              key='git_hash'              ) # Git hash for the data generation code.
        AddMetaDataWithReference(h5_file,cwd=outdir,value=configurator.GetPythiaConfigFileContents(pythia_config),    key='pythia_config'         ) # Pythia configuration (except for "\hat{p_T}", which is handled externally)

        # TODO: Might want to think about offering the ability to split the HepMC3 and Delphes files too?
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