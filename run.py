import sys,os,pathlib,time,datetime
import argparse as ap
import subprocess as sub
from util.generation import Generator
from util.delphes import BuildDelphes, HepMC3ToDelphes
from util.conversion import Processor, RemoveFailedFromHDF5, SplitH5, AddEventIndices, ConcatenateH5, MergeStatsInfo
from util.config import GetDelphesConfig, GetDelphesCard

def none_or_str(value): # see https://stackoverflow.com/a/48295546
    if value == 'None':
        return None
    return value

# TODO: We might want to move this to somewhere in the util library, to clean up this file.
def CompressHepMC(files, delete=True, cwd=None):
    for file in files:
        compress_file = file.replace('.hepmc','.tar.bz2')
        if(cwd is not None): compress_file = '{}/{}'.format(cwd,compress_file)
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],file.split('/')[-1]]
        # if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)
        if(delete):
            sub.check_call(['rm',file.split('/')[-1]],shell=False,cwd=cwd)
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n',          '--nevents',           type=int,          required=True,            help='Number of events per pt bin.')
    parser.add_argument('-p',          '--ptbins',            type=int,          required=True, nargs='+', help='Transverse momentum bin edges.')
    parser.add_argument('-o',          '--outfile',           type=str,          default='events.h5',      help='Output HDF5 file name.')
    parser.add_argument('-O',          '--outdir',            type=none_or_str,  default=None,             help='Output directory.')
    parser.add_argument('-g',          '--generation',        type=int,          default=True,             help='Whether or not to do event generation.')
    parser.add_argument('-s',          '--sep_truth',         type=int,          default=True,             help='Whether or not to store truth-level particles in separate arrays.')
    parser.add_argument('-ns',         '--n_sep_truth',       type=int,          default=-1,               help='How many truth particles to save in separate arrays -- will save the first n as given by the truth selection.')
    parser.add_argument('-d',          '--diagnostic_plots',  type=int,          default=False,            help='Whether or not to make diagnostic plots.')
    parser.add_argument('-v',          '--verbose',           type=int,          default=0,                help='Verbosity.')
    parser.add_argument('-h5',         '--hdf5',              type=int,          default=1,                help='Whether or not to produce final HDF5 files. If false, stops after HepMC or Delphes/ROOT file production.')
    parser.add_argument('-f',          '--force',             type=int,          default=0,                help='Whether or not to force generation -- if true, will possibly overwrite existing HepMC files in output directory.')
    parser.add_argument('-c',          '--compress',          type=int,          default=0,                help='Whether or not to compress HepMC files.')
    parser.add_argument('-cd',         '--clean_delphes',     type=int,          default=0,                help='Whether or not to clean up DELPHES/ROOT files.')
    parser.add_argument('-rng',        '--rng',               type=int,          default=None,             help='Pythia RNG seed. Will override the one provided in the config file.')
    parser.add_argument('-npc',        '--nentries_per_chunk',type=int,          default=int(1e4),         help='Number of entries to process per chunk, for jet clustering & conversion to HDF5.')
    parser.add_argument('-pb',         '--progress_bar',      type=int,          default=1,                help='Whether or not to print progress bar during event generation')
    parser.add_argument('-sp',         '--split',             type=int,          default=1,                help='Whether or not to split HDF5 file into training/validation/testing files.')
    parser.add_argument('-tf',         '--train_fraction',    type=float,        default=0.7,              help='Fraction of events to place in the training file.')
    parser.add_argument('-vf',         '--val_fraction',      type=float,        default=0.2,              help='Fraction of events to place in the validation file.')
    parser.add_argument('-co',         '--compression_opts',  type=int,          default=7,                help='Compression option for final HDF5 file (0-9). Higher value means more compression.')
    parser.add_argument('-separate_h5','--separate_h5',       type=int,          default=1,                help='Whether or not to make separate HDF5 files for each pT bin.')
    parser.add_argument('-pc',         '--pythia_config',     type=none_or_str,  default=None,             help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-debug',      '--debug',             type=int,          default=0,                help='If > 0, will record the full final-state (i.e. before jet clustering/selection) in a separate key.')
    args = vars(parser.parse_args())

    start_time = time.time()

    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
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
    nbins = len(pt_bin_edges) - 1
    pythia_rng = args['rng']
    nentries_per_chunk = args['nentries_per_chunk']
    progress_bar = args['progress_bar']
    compression_opts = args['compression_opts']
    separate_h5 = args['separate_h5'] > 0
    pythia_config = args['pythia_config']
    debug = args['debug'] > 0

    split_files = args['split'] > 0
    train_frac = args['train_fraction']
    val_frac = args['val_fraction']
    test_frac = 1. - train_frac - val_frac
    if(test_frac < 0. and split_files):
        print('Error: Requested training fraction and validation fraction sum to more than 1, this leaves no events for the test file.')
        assert(False)

    # Setting the verbosity for the HDF5 conversion.
    # If there are many events it might take a bit, so some printout
    # is helpful to monitor the progress.
    h5_conversion_verbosity = 0
    if(nevents_per_bin >= 100): h5_conversion_verbosity = 1
    elif(nevents_per_bin >= 10000): h5_conversion_verbosity = 2

    use_delphes = GetDelphesConfig()

    # Keep track of some files we create.
    jet_files = []
    truth_files = []
    hist_files = []
    stat_files = []
    final_state_truth_overlap_indices_files = [] # only used in certain cases, these files record indices of particles that show up in both truth and final-state selections (indices w.r.t. final-state HepMC files)

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
    comm = ['cp','config/config.py','{}/config.py'.format(outdir)]
    sub.check_call(comm)

    if(do_generation):
        if(verbose):
            print('\n=================================')
            print('Running Pythia8 event generation.')
            print('=================================\n')

        if(pythia_rng is not None):
            print('\tSetting Pythia RNG seed to {}. (overriding config)'.format(pythia_rng))
        if(pythia_config is not None):
            print('\tSetting Pythia process configuration from {}. (overriding config)'.format(pythia_config))

        if(verbose):
            print('\tGenerating {} events per \\hat{p_T} bin, with the following bin edges (in GeV):')
            for bin_edge in pt_bin_edges:
                print('\t\t{}'.format(bin_edge))

    for i in range(nbins):
        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.
        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]
        hep_file = 'events_{}-{}.hepmc'.format(pt_min,pt_max)

        generator = Generator(pt_min,pt_max, pythia_rng,pythia_config_file=pythia_config)
        generator.SetOutputDirectory(outdir)

        hist_filename = 'hists_{}.root'.format(i)
        hist_files.append(hist_filename)
        generator.SetHistFilename(hist_filename)

        stat_filename = 'stats_{}.h5'.format(i)
        stat_files.append(stat_filename)
        generator.SetStatsFilename(stat_filename)

        generator.SetFilename(hep_file)
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

        overlap_index_file = generator.GetIndexOverlapFilename()
        final_state_truth_overlap_indices_files.append(overlap_index_file)

        if(use_delphes): # Case 1: Using Delphes
            if(verbose and i == 0): print('Running Delphes on HepMC files.')
            # Pass the HepMC file to Delphes. Will output a ROOT file.
            delphes_dir = BuildDelphes() # build Delphes if it does not yet exist
            hep_file_with_extension = '{}/{}'.format(outdir,hep_file)
            delphes_card = GetDelphesCard() # will default to the ATLAS card that is shipped with Delphes
            delphes_file = hep_file.replace('.hepmc','.root')
            delphes_file = HepMC3ToDelphes(hepmc_file=hep_file, output_file=delphes_file, delphes_dir=delphes_dir, cwd=outdir, delphes_card=delphes_card)
            jet_files.append(delphes_file)

            # We can now compress the HepMC file. Once it's compressed, we delete the original file to recover space.
            if(compress_hepmc): CompressHepMC([hep_file_with_extension],True)

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
    processor = Processor(use_delphes)
    processor.SetOutputDirectory(outdir)
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
    else:
        # Slightly more complex way of running things -- we first create pT-binned HDF5 files, corresponding with the binning of the
        # HepMC/Delphes files. These binned files are then concatenatd into a single HDF5 file with all events.
        # Making the binned files first may allow us to (more easily) run some algorithms where we may want to compare across
        # the HepMC/Delphes and HDF5 files, since they are all binned it'll be easy to match events across files.
        h5_files = []
        print('\tProducing separate HDF5 files for each pT bin, and then concatenating these.')
        delete_individual_h5 = False
        nentries_per_chunk = int(nentries_per_chunk/nbins)
        for i in range(nbins):
            pt_min = pt_bin_edges[i]
            pt_max = pt_bin_edges[i+1]
            jet_file = [jet_files[i]]
            truth_file = [truth_files[i]]
            h5_file_individual = h5_file.replace('.h5','_{}-{}.h5'.format(pt_min,pt_max))
            processor.SetProgressBarPrefix('\tClustering jets & preparing data for pT bin [{},{}]:'.format(pt_min,pt_max))

            processor.Process(jet_file,truth_file,h5_file_individual,verbosity=h5_conversion_verbosity,nentries_per_chunk=nentries_per_chunk)

            # To each HDF5 event file, we will add event indices. These may be useful/necessary for the post-processing step.
            # We will remove these indices when concatenating files, since as one of our last steps we'll add indices again
            # but with respect to the full event listing (not just w.r.t. events in each generation pT bin).
            AddEventIndices(h5_file_individual,cwd=outdir,copts=compression_opts)

            # Optional post-processing. Any post-processing steps have been configured in the config file, config/config.py.
            processor.PostProcess(jet_file,[h5_file_individual],[final_state_truth_overlap_indices_files[i]])

            h5_files.append('/'.join((outdir,h5_file_individual)))

        print('\n\tConcatenating HDF5 files. Will drop the "event_idx" key, this was used internally for any post-processing steps.\n\tIndices will be recomputed and added at the end.')
        ConcatenateH5(h5_files,'/'.join((outdir,h5_file)),copts=compression_opts,delete_inputs=delete_individual_h5,ignore_keys=['event_idx'],verbose=True)

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
    delete_individual_stats = False
    try: ConcatenateH5(stat_files,stats_filename,cwd=outdir,delete_inputs=delete_individual_stats, copts=9)
    except: pass # as long as the full stats file exists, it's okay if the individual ones were deleted already

    try:
        MergeStatsInfo(h5_file,stats_filename,cwd=outdir,delete_stats_file=False, copts=9)
    except:
        print('Warning: Stats information not found!')
        pass

    # Now also compress the truth files.
    if(compress_hepmc): CompressHepMC(truth_files,True,cwd=outdir)

    # Add some event indices to our dataset.
    print('\tAdding event indices to file {}.'.format('/'.join((outdir,h5_file))))
    AddEventIndices(h5_file,cwd=outdir,copts=compression_opts)

    # Remove any failed events (e.g. detector-level events with no jets passing cuts).
    print('\tRemoving any failed events from file {}. These may be events where there weren\'t any jets passing the requested cuts.'.format('/'.join((outdir,h5_file))))
    RemoveFailedFromHDF5(h5_file,cwd=outdir)

    # TODO: Might want to think about offering the ability to split the HepMC3 and Delphes files too?
    #       Could be useful for certain post-processing where we need to access those files, and we
    #       want to work on the split files (insead of the full events.h5).
    #       Before splitting those, we'd want to effectively join them together first, they are currently
    #       split by pT bin.
    if(split_files):
        # Now split the HDF5 file into training, testing and validation samples.
        split_ratio = (train_frac,val_frac,test_frac)
        print("\tSplitting HDF5 file {} into training, validation and testing samples:".format('/'.join((outdir,h5_file))))
        split_ratio_sum = train_frac + val_frac + test_frac
        train_name = 'train.h5'
        val_name = 'valid.h5'
        test_name = 'test.h5'
        SplitH5(h5_file, split_ratio,cwd=outdir,copts=compression_opts, train_name=train_name,val_name=val_name,test_name=test_name,verbose=True)

    # Optionally delete the full HDF5 file.
    delete_h5 = False # TODO: Make this configurable
    if(delete_h5):
        comm = ['rm',h5_file]
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