import sys,os,glob,uuid,pathlib
import argparse as ap
import subprocess as sub
import util.utils as utils

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n',                '--nevents',         type=int, help='Number of events per pt bin.', required=True)
    parser.add_argument('-p',                '--ptbins',          type=int, help='Transverse momentum bin edges.', required=True, nargs='+', )
    parser.add_argument('-O',                '--outdir',          type=str, help='Output directory for the jobs', default=None)
    parser.add_argument('-R',                '--rundir',          type=str, help='Run directory -- where the condor jobs will go.', default='run0')
    parser.add_argument('-s',                '--sep_truth',       type=int, help='Whether or not to store truth-level particles in separate arrays.', default=1)
    parser.add_argument('-ns',               '--n_sep_truth',     type=int, help='How many truth particles to save in separate arrays -- will save the first n as given by the truth selection.', default=-1)
    parser.add_argument('-h5',               '--hdf5',            type=int, help='Whether or not to produce final HDF5 files. If false, stops after HepMC or Delphes/ROOT file production.',default=1)
    parser.add_argument('-rng',              '--rng',             type=int, help='Pythia RNG seed offset -- seed will be offset + job number.',default=1)
    parser.add_argument('-sp',               '--split',           type=int, help='Whether or not to split HDF5 file into training/validation/testing files.',default=1)
    parser.add_argument('-pc',               '--pythia_config',   type=str, help='Path to Pythia configuration template (for setting the process).',default=[None], nargs='+')
    parser.add_argument('-N',                '--Njobs',           type=int, help='Number of jobs per config.',default=1)
    parser.add_argument('-sq',               '--shortqueue',      type=int, help='Whether or not to use the condor short queue (for UChicago Analysis Facility).',default=0)
    parser.add_argument('-ncpu',             '--n_cpu',           type=int, help='Number of (logical) CPU cores per job.',default=1)
    parser.add_argument('-nblas',            '--n_openblas',      type=int, help='Sets the $OPENBLAS_NUM_THREADS variable (used for numpy multithreading). Advanced usage.',default=16)
    parser.add_argument('-event_idx_offset', '--event_idx_offset',type=int, help='Initial offset for event_idx. Advanced usage.',default=0)
    parser.add_argument('-batch_name',       '--batch_name',      type=str, help='Name for the condor batch. If left empty, unused.',default=None)
    parser.add_argument('-git',              '--git_option',      type=int, help='If < 0, ships code directly from this repo. If == 0, runs \'git clone\' within each job. If > 0, runs \'git clone\' locally and ships result to jobs.',default=0)
    parser.add_argument('-branch',           '--branch',          type=str, help='If using git for jobs, which branch to clone. Defaults to current.', default=None)
    args = vars(parser.parse_args())

    nevents = args['nevents']
    ptbins = args['ptbins']
    rundir = args['rundir']
    outdir = args['outdir']
    if(outdir is None):
        outdir = rundir + '/output'
    sep_truth = args['sep_truth']
    n_sep_truth = args['n_sep_truth']
    do_h5 = args['hdf5']
    rng_seed = args['rng']
    split = args['split']
    pythia_configs = args['pythia_config']
    njobs = args['Njobs']
    short_queue = args['shortqueue'] > 0
    ncpu = args['n_cpu']
    nblas = args['n_openblas']
    event_idx_offset_initial = args['event_idx_offset']
    batch_name = args['batch_name']
    git_option = args['git_option']
    git_branch = args['branch']

    rundir = str(pathlib.Path(rundir).absolute())
    outdir = str(pathlib.Path(outdir).absolute())

    # Prepare the job and output directories.
    os.makedirs(rundir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    # Prepare a plaintext file with all the different sets of arguments, for the various jobs.
    jobs_filename = '{}/arguments.txt'.format(rundir)
    ptbins_str = ','.join([str(x) for x in ptbins])
    template = '{} {} {} {} {} {} {} {} {}'.format(
        nevents, # $1 Number of events per pT bin.
        ptbins_str, #$2 pT bins (list of bin edges)
        sep_truth, #$3 whether or not to save separate truth particle keys (0 or 1)
        n_sep_truth, #$4 number of separate truth particles to save (i.e. save the first n from the list)
        do_h5, #$5 whether or not to do jet clustering and make HDF5 files (HepMC -> HDF5). If not, stops after making the HepMC files.
        '{}', #$6 RNG seed for generation. (can be used to overwrite the builtin config file)
        split, #$7 whether or not to split final HDF5 file into train/validation/test files. Only relevant if making the HDF5 file.
        '{}', #$8 Pythia config (can be used to overwrite the builtin config file)
        '{}' #$9 Event index offset.
    ) # double-quotes for silly condor argument syntax
    with open(jobs_filename,'w') as f:
        rng_counter = 0
        for i,pythia_config in enumerate(pythia_configs):
            for j in range(njobs):
                rng = rng_seed + rng_counter
                event_idx_offset = event_idx_offset_initial + j * (len(ptbins) - 1) * nevents # len(ptbins) gives # of bin edges, which is number of bins + 1
                command_arguments = template.format(rng,pythia_config,event_idx_offset) + '\n'
                f.write(command_arguments)
                rng_counter += 1

    this_dir = os.path.dirname(os.path.abspath(__file__))

    config_file = '{}/../config/config.py'.format(this_dir)
    comm = ['cp',config_file,rundir]
    sub.check_call(comm)

    # Determine how the package will be handled -- do we clone it as part of the job, clone it here and ship it,
    # or just ship this existing repo?
    payload_string = ''
    if(git_branch is None): git_branch = utils.GetGitBranch()

    if(git_option < 0): git_option = -1
    elif(git_option > 0): git_option = 1

    if(git_option == -1):
        # We have to ship the local repo.
        payload = 'payload.tar.bz2'
        gitdir = 'HEPData4ML'
        utils.PreparePayload(rundir,payload,gitdir)
        payload_string = ', ../{}'.format(payload)

    elif(git_option == 1):
        # We have to do a local git clone
        payload = 'payload.tar.bz2'
        gitdir = 'HEPData4ML'
        utils.PreparePayloadFromClone(rundir,payload,gitdir,branch=git_branch)
        payload_string = ', ../{}'.format(payload)

    else:
        print(63 * '-')
        print('HEPData4ML repository will be cloned internally by condor jobs.')
        print(63 * '-')

    # Now we need to make a condor submission file for this set of jobs. It will be based off a template file.
    condor_template = '{}/util/condor_templates/condor_template.sub'.format(this_dir)
    with open(condor_template,'r') as f:
        condor_submit_lines = f.readlines()

    short_queue_line = ' +queue="short"'
    if(not short_queue): short_queue_line = '#' + short_queue_line

    batch_line = 'batch_name = {}'
    if(batch_name is not None): batch_line = batch_line.format(batch_name)
    else: batch_line = '#' + batch_line

    for i,line in enumerate(condor_submit_lines):
        condor_submit_lines[i] = condor_submit_lines[i].replace("$BATCH_NAME",batch_line + '\n')
        condor_submit_lines[i] = condor_submit_lines[i].replace("$OUTDIR",outdir)
        condor_submit_lines[i] = condor_submit_lines[i].replace("$N_THREAD",str(nblas))
        condor_submit_lines[i] = condor_submit_lines[i].replace("$N_CPU",str(ncpu))
        condor_submit_lines[i] = condor_submit_lines[i].replace('$SHORT_QUEUE',short_queue_line + '\n')
        condor_submit_lines[i] = condor_submit_lines[i].replace("$GIT_OPTION",str(git_option))
        condor_submit_lines[i] = condor_submit_lines[i].replace("$PAYLOAD_STRING",payload_string)
        condor_submit_lines[i] = condor_submit_lines[i].replace("$GIT_BRANCH",git_branch)

    condor_submit_file = '{}/condor.sub'.format(rundir)
    with open(condor_submit_file,'w') as f:
        for line in condor_submit_lines:
            f.write(line)

    # Copy the condor executable to the job folder.
    executable = '{}/util/condor_job.sh'.format(this_dir)
    comm = ['cp',executable,rundir]
    sub.check_call(comm)

    print('Prepared jobs. Using RNG seeds:')
    print('\tstart = {}'.format(rng_seed))
    print('\t  end = {}'.format(rng))

if(__name__=='__main__'):
    main(sys.argv)