import sys,os,glob,uuid,pathlib
import argparse as ap
import subprocess as sub

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n', '--nevents', type=int, help='Number of events per pt bin.', required=True)
    parser.add_argument('-p', '--ptbins', type=int, nargs='+', help='Transverse momentum bin edges.', required=True)
    parser.add_argument('-O', '--outdir', type=str, help='Output directory for the jobs', default=None)
    parser.add_argument('-R', '--rundir', type=str, help='Run directory -- where the condor jobs will go.', default='run0')
    parser.add_argument('-s', '--sep_truth',type=int, help='Whether or not to store truth-level particles in separate arrays.', default=1)
    parser.add_argument('-ns', '--n_sep_truth',type=int, help='How many truth particles to save in separate arrays -- will save the first n as given by the truth selection.', default=-1)
    parser.add_argument('-h5', '--hdf5',type=int,help='Whether or not to produce final HDF5 files. If false, stops after HepMC or Delphes/ROOT file production.',default=1)
    parser.add_argument('-rng','--rng',type=int,help='Pythia RNG seed offset -- seed will be offset + job number.',default=1)
    parser.add_argument('-sp','--split',type=int,default=1,help='Whether or not to split HDF5 file into training/validation/testing files.')
    parser.add_argument('-pc','--pythia_config',type=str, nargs='+', default=[None],help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-N','--Njobs',type=int,default=1,help='Number of jobs per config.')
    parser.add_argument('-sq', '--shortqueue', type=int,default=0,help='Whether or not to use the condor short queue (for UChicago Analysis Facility).')
    parser.add_argument('-nblas','--n_openblas',type=int,default=16,help='Sets the $OPENBLAS_NUM_THREADS variable (used for numpy multithreading). Advanced usage.')
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
    nblas = args['n_openblas']

    rundir = str(pathlib.Path(rundir).absolute())
    outdir = str(pathlib.Path(outdir).absolute())

    # Prepare the job and output directories.
    os.makedirs(rundir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    # Prepare a plaintext file with all the different sets of arguments, for the various jobs.
    jobs_filename = '{}/job_arguments.txt'.format(rundir)
    ptbins_str = ','.join([str(x) for x in ptbins])
    template = '{} {} {} {} {} {} {} {}'.format(nevents,ptbins_str,sep_truth, n_sep_truth, do_h5, '{}',split,'{}') # double-quotes for silly condor argument syntax
    with open(jobs_filename,'w') as f:
        rng_counter = 0
        for i,pythia_config in enumerate(pythia_configs):
            for j in range(njobs):
                rng = rng_seed + rng_counter
                command_arguments = template.format(rng,pythia_config) + '\n'
                f.write(command_arguments)
                rng_counter += 1

    this_dir = os.path.dirname(os.path.abspath(__file__))

    # We also ship the config.py file separately so that it can be easily modified outside the payload.
    # There is one inside the payload too but it will be overwritten by this one within the job.
    config_file = '{}/../config/config.py'.format(this_dir)
    comm = ['cp',config_file,rundir]
    sub.check_call(comm)

    # Now we need to make a condor submission file for this set of jobs. It will be based off a template file.
    condor_template = '{}/util/condor_templates/condor_template.sub'.format(this_dir)
    with open(condor_template,'r') as f:
        condor_submit_lines = f.readlines()

    short_queue_line = ' +queue="short"'
    if(not short_queue): short_queue_line = '#' + short_queue_line

    for i,line in enumerate(condor_submit_lines):
        condor_submit_lines[i] = condor_submit_lines[i].replace("$OUTDIR",outdir)
        condor_submit_lines[i] = condor_submit_lines[i].replace("$N_THREAD",nblas)
        condor_submit_lines[i] = condor_submit_lines[i].replace('$SHORT_QUEUE',short_queue_line + '\n')

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