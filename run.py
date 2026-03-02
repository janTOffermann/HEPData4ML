import sys,os,pathlib,re,shlex,glob
import argparse as ap
import subprocess as sub
from contextlib import nullcontext
from util.generation.generation import PythiaGenerator
from util.simulation.simulation import DelphesSimulator
from util.reconstruction.conversion import Processor
from util.hepmc.hepmc import CompressHepMC
from util.config.config import Configurator,GetConfigFileContent, GetConfigDictionary
from util.config.args import parse_mc_steps, FloatListAction, none_or_str
from util.metadata.meta import MetaDataHandler

from util.misc.timing import BasicTimer, profiling_context

# Convenience function for file naming
def float_to_str(value):
    value_str = str(value)
    value_str = value_str.replace('.',',')
    return re.sub(',0$','',value_str)

def main(args):
    parser = ap.ArgumentParser()

    parser.add_argument('-n',            '--nevents',           type=int,          required=True,            help='Number of events per pt bin.')
    parser.add_argument('-steps',        '--steps',             type = parse_mc_steps, default = 'generation pileup simulation reconstruction', help='Comma- or space-separated list of step. Options are [generation,pileup, simulation,reconstruction].')
    parser.add_argument('-p',            '--ptbins',            action = FloatListAction, default=[-1,-1], nargs='*',    help='Transverse momentum bin edges, for outgoing particles of the hard process. Can be a list of floats, or a string of comma- or space-separated floats. In GeV.')
    parser.add_argument('-o',            '--outfile',           type=str,          default='events',      help='Output ntuple file name, excluding extension.')
    parser.add_argument('-O',            '--outdir',            type=none_or_str,  default=None,             help='Output directory.')
    parser.add_argument('-v',            '--verbose',           action='store_true',                         help='Verbosity.')
    parser.add_argument('-f',            '--force',             type=int,          default=0,                help='Whether or not to force generation -- if true, will possibly overwrite existing HepMC files in output directory.')
    parser.add_argument('-c',            '--compress',          type=int,          default=0,                help='Whether or not to compress HepMC files from the generation step.')
    parser.add_argument('-npc',          '--nentries_per_chunk',type=int,          default=int(1e4),         help='Number of entries to process per chunk, for jet clustering & conversion to HDF5.')
    parser.add_argument('-pb',           '--progress_bar',      type=int,          default=1,                help='Whether or not to print progress bar during event generation')
    parser.add_argument('-co',           '--compression_opts',  type=int,          default=7,                help='Compression option for final HDF5 file (0-9). Higher value means more compression.')
    parser.add_argument('-pc',           '--pythia_config',     type=none_or_str,  default=None,             help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-index_offset', '--index_offset',      type=int,          default=0,                help='Offset for Event.Index.')
    parser.add_argument('-config',       '--config',            type=str,          default=None,             help='Path to configuration Python file. Default will use config/config.py .')

    # Overrides for things set in the config file.
    # TODO: There are a lot of settings in the config. Is there some nice programmatic way
    #       to make lots of corresponding arguments?
    parser.add_argument('-rng',          '--rng',               type=int,          default=None,             help='Pythia RNG seed. Overides the config file.')
    parser.add_argument('-pileup',       '--pileupFiles',       type=str,          default=None,             help='Glob-compatible string for input pileup files, for the pileup step. Overides the config file.')

    # Flag for when running as a parallelized job. Shouldn't need to be touched by user.
    parser.add_argument('-condor','--condor', action='store_true',help='Flag to be set when this is an HTCondor job; for advanced usage (you typically should *not* set this, for internal use.).')
    parser.add_argument('-condor_job_number','--condor_job_number', type=int, default= None,help='Number of this HTCondor job; for advanced usage (you typically should *not* set this, for internal use.).')
    parser.add_argument('-n_condor_jobs','--n_condor_jobs', type=int, default= None,help='Number of HTCondor jobs; for advanced usage (you typically should *not* set this, for internal use.).')

    # Time profiling
    parser.add_argument('-profile','--profile',action='store_true')

    args = vars(parser.parse_args())

    timer = BasicTimer()
    timer.start_main()

    steps = args['steps']
    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
    ntuple_file = args['outfile']
    outdir = args['outdir']
    verbose = args['verbose']
    compress_hepmc = args['compress']
    force = args['force']
    nentries_per_chunk = args['nentries_per_chunk']
    progress_bar = args['progress_bar']
    pythia_config = args['pythia_config']
    index_offset = args['index_offset']
    config_file = args['config']

    nbins = len(pt_bin_edges) - 1

    pythia_rng = args['rng']
    pileup_input_files = args['pileupFiles']

    # Arguments to be used by HTCondor jobs.
    # TODO: Maybe find another way to handle these? A wrapper script for condor?
    is_condor_job = args['condor']
    condor_job_number = args['condor_job_number']
    n_condor_jobs = args['n_condor_jobs']

    do_profile = args['profile']

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

    metadata_handler = MetaDataHandler(configurator)

    # # Set up FastJet -- we will need this later on (except for the special use case of no jet clustering!).
    # # To keep our printouts clean, we are initializing FastJet here instead of later on in a loop, so that
    # # we can get the FastJet banner printout out of the way. We remove the banner with some ANSI printing
    # # hackery, since it's really not useful and clutters up our printout (we acknowledge the use in the
    # # documentation, having this unavoidable printout for a single package's use is quite gratuitous).
    # if('reconstruction' in steps):
    #     print(13 * '\n')
    #     dummy_processor = Processor(configurator)
    #     line_up = '\033[1A'
    #     line_clear = '\x1b[2K'
    #     for i in range(13):
    #         print(line_up, end=line_clear)

    # Setting the verbosity for the HDF5 conversion.
    # If there are many events it might take a bit, so some printout
    # is helpful to monitor the progress.
    ntuple_verbosity = 0
    if(nevents_per_bin >= 100): ntuple_verbosity = 1
    elif(nevents_per_bin >= 10000): ntuple_verbosity = 2

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

    # Optional time profiling.
    profile_context = profiling_context() if do_profile else nullcontext()

    with profile_context as profiler:
        #=========================
        # STEP 0: Metadata
        #=========================
        # We can now stash some things away in the metadata handler.
        metadata_handler.AddElement('Metadata.CommandLineArguments'," ".join(map(shlex.quote, sys.argv[1:])))
        metadata_handler.AddElement('Metadata.ConfigurationFile','\n'.join(GetConfigFileContent(config_file)))

        #=========================
        # STEP 1: Generation
        #=========================
        hepmc_files = []
        if('generation' in steps):
            timer.start_timestamp('generation')

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

        # Stash some stuff in metadata.
        # TODO: Would be nice to do this within the PythiaGenerator()?
        metadata_handler.AddElement('Metadata.Generation.PythiaRandomSeed',pythia_rng)
        metadata_handler.AddElement('Metadata.Generation.PythiaConfiguration',configurator.GetPythiaConfigFileContents(pythia_config))

        print()
        for i in range(nbins):
            # Generate a HepMC file containing our events.
            # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
            # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.
            pt_min = pt_bin_edges[i]
            pt_max = pt_bin_edges[i+1]

            hepmc_extension = 'hepmc'
            if(configurator.GetHepMCFormat().lower() == 'root'):
                hepmc_extension = 'hepmc.root'

            hep_file = 'events_{}.{}'.format(i,hepmc_extension)

            generator = PythiaGenerator(pt_min,pt_max, configurator, pythia_rng,pythia_config_file=pythia_config)
            generator.SetMetadataHandler(metadata_handler)

            generator.SetOutputDirectory(outdir)
            generator.SetFilename(hep_file)
            generator.SetProgressBar(progress_bar)

            hepfile_exists = pathlib.Path('{}/{}'.format(outdir,hep_file)).exists()
            generate = 'generation' in steps
            if(hepfile_exists and not force):
                print('\tHepMC3 file {}/{} already found, skipping its generation.'.format(outdir,hep_file))
                generate = False

            if(generate):
                generator.Generate(nevents_per_bin)

            hepmc_files.append(hep_file)

        if('generation' in steps):
            timer.end_timestamp('generation')

        #===================================
        # STEP 1.5: Metadata for HepMC3/ROOT
        #===================================
        # If HepMC3/ROOT files were generated, we stash the metadata in them too.
        # This is particularly helpful for the pileup step, as it can read this metadata
        # if using these files in another run.
        # Note that if plaintext HepMC3 was used, this feature is unavailable, and it may
        # make the final dataset harder to reproduce (since you won't directly have info on
        # how the pileup_handlers' input files were created!).
        if(hepmc_extension == 'root'):
            metadata_handler.StashMetaDataInROOTFile(hepmc_files,cwd=outdir)

        #===============================
        # STEP 2: Pileup (optional)
        #===============================
        # In this step, we produce pileup HepMC/ROOT "sidecar" files
        # (or just fetch existing ones). Note that producing these from scratch
        # is a relatively slow process (specifically, combining HepMC3 events is slow).
        # Thus it is preferrable to produce these in a dedicated run, and then later
        # simply use a pileup_handler that fetches these pre-mixed events.
        pileup_handler = None
        pileup_files = None
        if('pileup' in steps):
            timer.start_timestamp('pileup')
            pileup_handler = configurator.GetPileupHandler()

            if(pileup_handler is not None): # if it's set to None, we just skip the pileup step
                pileup_handler.SetInputDirectory(outdir)
                pileup_handler.SetOutputDirectory(outdir)
                pileup_handler.SetMetadataHandler(metadata_handler)
                pileup_handler.SetConfigurator(configurator)
                pileup_handler.SetHTCondorInfo(is_condor_job,condor_job_number,n_condor_jobs)

                if(pileup_input_files is not None): # overriding config file
                    pileup_handler.SetPileupFiles(pileup_input_files)

                # Special cases, where we use the Pythia RNG seed
                if(pileup_handler.GetRNGSeed() < 0): # Case 1: Seed in the config file is negative.
                    pileup_handler.SetRNGSeed(pythia_rng)

                elif(args['rng'] is not None): # Case 2: The Pythia RNG seed was specified at command line -- in practice we may want to then use this for pileup too (e.g. HTCondor usage).
                    pileup_handler.SetRNGSeed(pythia_rng)

                # Now, we actually run the pileup production.
                pileup_files = pileup_handler.Process(hepmc_files)

                # TODO: Now we fetch some information from the pileup handler, that will propagate into the final dataset:
                # info on the number of interactions per bunch crossing, the actual indices of pileup events used,
                # plus some other optional pieces of info that depend on what handler we used and how it was configured.

                #TODO Rework this bit -- pileup handling was changed since this snippet was written
                # #===================================
                # # STEP 2.5: Metadata for HepMC3/ROOT
                # #===================================
                # # Similar to Step 1.5 -- we again add metadata to HepMC3/ROOT files.
                # # We need to do this again since, if we're doing this pileup step,
                # # the HepMC3/ROOT files have been overwritten. Plus, there's more
                # # metadata to add to them now.
                # if(hepmc_extension == 'root'):
                #     metadata_handler.StashMetaDataInROOTFile(hepmc_files,cwd=outdir)
            timer.end_timestamp('pileup')

        #===============================
        # STEP 3: Simulation (optional)
        #===============================
        #
        # Here is where we invoke detector simulation.
        # For now, the only option is fast detector simulation
        # with Delphes.
        # Note that some of the of the code in Step 4 is specialized
        # for handling Delphes output; if other detector sims are
        # introduced this may require some add-ons for Step 4.
        #
        delphes_files = []
        simulation_type=None
        if('simulation' in steps):
            timer.start_timestamp('simulation')
            simulator = None
            simulation_type = configurator.GetSimulationType()
            if(simulation_type == 'delphes'):
                sim_logfile = '{}/delphes.log'.format(outdir)
                simulator = DelphesSimulator(configurator,outdir,logfile=sim_logfile)

            if(simulator is not None):
                simulator.SetMetadataHandler(metadata_handler)
                simulator.SetInputs(hepmc_files)
                if(pileup_files is not None):
                    simulator.SetPileupInputs(pileup_files)
                simulator.Process()
                delphes_files = simulator.GetOutputFiles()
            timer.end_timestamp('simulation')
        else: # try to pick up any available delphes files #TODO: Make this more robust
            delphes_files = glob.glob('{}/*.delphes.root'.format(outdir))
            delphes_files = [x.replace(outdir + '/','') for x in delphes_files]


        #=========================================================
        # STEP 4: HDF5 conversion + Reconstruction/Post-processing
        #=========================================================
        #
        # A lot of stuff happens here. We produce the (ROOT/HDF5) n-tuples.
        # This is also where we'll run things like jet clustering, which
        # will act upon those n-tuples as an afterburner (or "post-processor")
        # and add new branches to them.
        #

        processor = None
        if('reconstruction' in steps):
            timer.start_timestamp('reconstruction')
            # Do reco and put everything into an n-tuple file.
            print('\nRunning reconstruction and producing final N-tuple output.\n')
            processor = Processor(configurator)
            processor.SetOutputDirectory(outdir)
            processor.SetMetadataHandler(metadata_handler)

            # List of ntuple filenames, for the temporary files
            # (one per input HepMC/Delphes). We'll merge them after the loop.
            ntuple_files = []

            # Create the filename for the unified ntuple file.
            ntuple_file = '{}.{}'.format(ntuple_file,processor.GetOutputExtension())

            if(verbose): print('\nProducing separate N-tuple files for each pT bin, and then concatenating these.')
            nentries_per_chunk = int(nentries_per_chunk/nbins)

            for i, hepmc_file in enumerate(hepmc_files):
                delphes_file = None
                pileup_file = None
                if(len(delphes_files) > 0):
                    delphes_file = delphes_files[i]
                if(pileup_files is not None):
                    pileup_file = pileup_files[i]

                # TODO: Rework this a little. Should just generically loop over HepMC files, since they might have an external source and not be pt-binned.
                ntuple_file_individual = hepmc_file.split('/')[-1].replace(hepmc_extension,processor.GetOutputExtension())

                processor.SetProgressBarPrefix('\tProducing N-tuple for file {}/{}:'.format(i+1,len(hepmc_files)))
                processor.ProcessFull(hepmc_file,pileup_file, delphes_file, ntuple_file_individual,verbosity=ntuple_verbosity)

                # # Add information from the pileup handler (if any).
                # # TODO: This may need a little reworking? The handling of filenames might be a little fragile.
                # if(pileup_handler is not None):
                #     pileup_handler.AddPileupInfoToH5(ntuple_file_individual,cwd=outdir,file_key=hepmc_file)

                ntuple_file_individual = '/'.join((outdir,ntuple_file_individual))
                ntuple_files.append(ntuple_file_individual)

            # Now, concatenate the ntuple files together.
            processor.ConcatenateNtuples(ntuple_files,ntuple_file,format='root',delete_inputs=True) # will prepend outdir to the output_file argument (this is all a little messy)

            #Cleanup: Compress the HepMC files.
            if(compress_hepmc): CompressHepMC(hepmc_files,True,cwd=outdir)

            # Now, add event indices to the datset. Each of the `ntuple_files` had them (in case useful
            # for any post-processing), but they are dropped when we call processor.ConcatenateNtuples.
            # Now, we add them across all the events in the concatenated set.
            if(index_offset < 0):
                index_offset = 0
            processor.AddEventIndices(ntuple_file,offset=index_offset)

            #======================================================
            # STEP 4.5: Metadata (into N-tuple).
            #======================================================
            # Now, add some metadata to the file.
            metadata_handler.AddMetaDataWithReferenceRoot(ntuple_file,cwd=outdir, tree_name=processor.GetTreeName())

            # TODO: Might want to think about offering the ability to split upstream files too? Could be complicated...
            # if(split_files):
            #     # Now split the HDF5 file into training, testing and validation samples.
            #     split_ratio = (train_frac,val_frac,test_frac)
            #     print("\tSplitting HDF5 file {} into training, validation and testing samples:".format('/'.join((outdir,ntuple_file))))
            #     train_name = 'train.h5'
            #     val_name = 'valid.h5'
            #     test_name = 'test.h5'
            #     SplitH5(ntuple_file, split_ratio,cwd=outdir,copts=compression_opts, train_name=train_name,val_name=val_name,test_name=test_name,verbose=True,seed=configurator.GetSplitSeed())

            # # Optionally delete the full N-tuple file.
            # if(delete_full):
            #     comm = ['rm','{}/{}'.format(outdir,ntuple_file)]
            #     sub.check_call(comm)
            timer.end_timestamp('reconstruction')
        timer.end_main()
        print('\n#############################')
        timer.summarize_time()
        # Give a further breakdown of the post-processing.
        if(processor is not None):
            if(processor.post_processing is not None):
                for i,post_proc in enumerate(processor.post_processing):
                    post_proc.SummarizeRuntime(level=2)
            print('\n#############################')


        if(profiler is not None): profiler.report()

if __name__ == '__main__':
    main(sys.argv)