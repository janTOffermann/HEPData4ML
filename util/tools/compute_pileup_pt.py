
import sys,os,pathlib,glob
import argparse as ap
from typing import Union
import subprocess as sub

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../..' )
from util.pileup.pileup import PileupOverlayPtFilter
from util.config.config import Configurator
from util.condor.condor_utils import FetchRequirements, GetGitBranch, PreparePayload, PreparePayloadFromClone

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFiles',type=str,help='Glob-compatible string for input HepMC3/ROOT pileup files.',required=True)
    parser.add_argument('-o','--outputDirectory',type=str,help='Output directory. If provided, will copy input files there first, and act on those.',default=None)
    parser.add_argument('-fastjet','--fastjetDir',type=str,help='FastJet installation directory. Defaults to `None`, which will use local install.',default=None)

    parser.add_argument('-condor','--condor',action='store_true',help='Whether or not to prep condor jobs to run this in parallel.')
    parser.add_argument('-r','--runDirectory',type=str,default=None,help='[condor only] Run directory.')
    parser.add_argument('-m','--mode',type=int,default=-1,help='[condor only] Running mode.')
    parser.add_argument('-n','--nFiles',type=int,default=-1,help='Number of files to run on.')


    args = vars(parser.parse_args())

    input_files = args['inputFiles']
    output_directory = args['outputDirectory']
    fastjet_dir = args['fastjetDir']
    condor = args['condor']
    run_directory = args['runDirectory']
    run_mode = args['mode']
    n_files = args['nFiles']

    # Create a configurator -- used under-the-hood by pileup_handler for some FastJet configuration.
    config_dictionary = {'reconstruction':{'fastjet_dir':fastjet_dir}}
    configurator = Configurator(config_dictionary=config_dictionary)

    if(output_directory is not None):

        input_files_new = []
        os.makedirs(output_directory,exist_ok=True)
        input_file_list = glob.glob(input_files,recursive=True)

        for i,file in enumerate(input_file_list):
            if(n_files > 0 and i >= n_files):
                break

            new_file = '{}/{}'.format(output_directory,file.split('/')[-1])
            if(pathlib.Path(new_file).exists()):
                print('Warning: Output file {} already exists.'.format(new_file))
            else:
                print('Copying input file {}/{} to output (where it will be modified by job).'.format(i+1,len(input_file_list)))
                command = ['cp',file,new_file]
                sub.check_call(command)
            input_files_new.append(new_file)
        input_files = input_files_new

    if(not condor):
        pileup_handler = PileupOverlayPtFilter(input_files,precompute=True)
        pileup_handler.SetConfigurator(configurator)
        pileup_handler.Initialize() # this will launch the computation of the leading jet pt
        return

    assert(run_directory is not None)

    condor_runner = CondorRunner()
    condor_runner.SetRunDirectory(run_directory)
    condor_runner.SetInputFiles(input_files)
    condor_runner.SetPayloadMode(run_mode)
    condor_runner.Prepare()

    return

class CondorRunner:

    def __init__(self):
        self.short_queue = False
        self.batch_name = 'HEPData4ML::compute_pileup_pt'
        self.requirements = None
        self.payload_mode = 0
        self.payload_string = ''
        self.this_dir = os.path.dirname(os.path.abspath(__file__))
        self.njobs = 0

        self.git_branch = None
        self.condor_template = None

    def SetPayloadMode(self,val:int):
        self.payload_mode = val

    def SetShortQueueFlag(self,val:bool):
        self.short_queue = val

    def SetInputFiles(self,val:Union[str,list]):
        if(isinstance(val,list)):
            self.input_files = val
        else:
            self.input_files = glob.glob(val,recursive=True)
        return

    def SetRunDirectory(self,val:str):
        self.rundir = val

    def SetCondorTemplate(self,val:str):
        path = pathlib.Path(val)

        if(not path.exists()):
            print('Warning: Condor template {} does not exist.'.format(val))
            assert False
        self.condor_template = str(pathlib.Path(val).absolute())

    def Prepare(self):

        # check that the template submission file exists
        if(self.condor_template is None):
            self.SetCondorTemplate('{}/../../util/condor/condor_templates/condor_template.sub'.format(self.this_dir))
        if(not pathlib.Path(self.condor_template).exists()):
            print('Error: Condor submission file template not found: {}'.format(self.condor_template))
            return

        # Prepare the job and output directories.
        os.makedirs(self.rundir,exist_ok=True)

        # Prepare a plaintext file with all the different sets of arguments, for the various jobs.
        arguments_filename = '{}/arguments.txt'.format(self.rundir)
        with open(arguments_filename,'w') as f:
            for i,input_file in enumerate(self.input_files):
                f.write('{}\n'.format(str(pathlib.Path(input_file).absolute())))
                self.njobs += 1

        # Prepare the payload & its options.
        self._payload()

        # Prepare the condor submission file.
        self._write_condor_submission_file()

        # Copy executable to the run directory
        self._copy_executable()

        # Make subdirs.
        self._prepare_subdirs()


    def _payload(self):
        ######################
        # Handling the payload
        ######################

        # Determine how the package will be handled. Do we:
        # 1) Clone it from git as part of the job itself?
        # 2) Clone it here in the interactive session, and ship it as an archive to the job?
        # 3) Ship this existing repository as an archive to the job?
        # 4) Point the jobs to this existing repository?
        #
        # Note that option 4 requires that the condor workers have access to this filesystem.
        # For systems such as the condor queue on BRUX (Brown University), where jobs aren't
        # shipped to some temporary directory but rather run directly out of "initialdir",
        # this might be preferrable or necessary.

        if(self.git_branch is None): self.git_branch = GetGitBranch()

        if(self.payload_mode == 1):
            print(63 * '-')
            print('HEPData4ML repository will be cloned internally by condor jobs.')
            print(63 * '-')

        elif(self.payload_mode == 2):
            # We have to do a local git clone
            payload = 'payload.tar.gz'
            gitdir = 'HEPData4ML'
            PreparePayloadFromClone(self.rundir,payload,gitdir,branch=self.git_branch)
            self.payload_string = ', ../{}'.format(payload)

        elif(self.payload_mode == 3):
            # The "local" run mode -- a bit special, originally designed for Brown University BRUX system.
            # We will point the workers at this repo.
            self.this_dir = os.path.dirname(os.path.abspath(__file__))
            gitdir = str(pathlib.Path('{}/../../'.format(self.this_dir)).absolute())
            self.payload_mode = gitdir
            self.git_branch = ''

        else:
            # We have to ship the local repo.
            payload = 'payload.tar.gz'
            gitdir = 'HEPData4ML'
            PreparePayload(self.rundir,payload,gitdir)
            self.payload_string = ', ../{}'.format(payload)

    def _write_condor_submission_file(self):

        ########################
        # Condor submission file
        ########################

        # Now we need to make a condor submission file for this set of jobs. It will be based off a template file.
        with open(self.condor_template,'r') as f:
            condor_submit_lines = f.readlines()

        # We are recycling a template file for the main script, so we change the arguments section
        for i,line in enumerate(condor_submit_lines):
            if('arguments' in line and 'queue' not in line):
                condor_submit_lines[i] = 'arguments               = "$(job_arguments) $(Process) $PAYLOAD_MODE $GIT_BRANCH"\n'
            elif('transfer_input_files') in line:
                condor_submit_lines[i] = 'transfer_input_files    = $PAYLOAD_STRING\n'

        # The short queue is something specific to the UChicago Analysis Facility condor queue.
        short_queue_line = ' +queue="short"'
        if(not self.short_queue): short_queue_line = '#' + short_queue_line

        batch_line = 'batch_name = {}'
        if(self.batch_name is not None): batch_line = batch_line.format(self.batch_name)
        else: batch_line = '#' + batch_line

        requirements_string = FetchRequirements(self.requirements)
        if(requirements_string is not None):
            if(len(requirements_string) > 0):
                requirements_string = 'requirements            = {}'.format(requirements_string)

        for i,line in enumerate(condor_submit_lines):
            condor_submit_lines[i] = condor_submit_lines[i].replace("$BATCH_NAME",batch_line + '\n')
            condor_submit_lines[i] = condor_submit_lines[i].replace('$ADDITIONS',short_queue_line + '\n')
            condor_submit_lines[i] = condor_submit_lines[i].replace("$PAYLOAD_MODE",str(self.payload_mode))
            condor_submit_lines[i] = condor_submit_lines[i].replace("$PAYLOAD_STRING",self.payload_string)
            condor_submit_lines[i] = condor_submit_lines[i].replace("$GIT_BRANCH",self.git_branch)
            condor_submit_lines[i] = condor_submit_lines[i].replace('$REQUIREMENTS',requirements_string)
            condor_submit_lines[i] = condor_submit_lines[i].replace('$N_CPU',str(1))


        condor_submit_file = '{}/condor.sub'.format(self.rundir)
        with open(condor_submit_file,'w') as f:
            for line in condor_submit_lines:
                f.write(line)


    def _copy_executable(self):
        # Copy the condor executable to the job folder.
        executable = '{}/../../util/condor/executables/compute_pileup_pt.sh'.format(self.this_dir)
        comm = ['cp',executable,self.rundir + '/condor_job.sh']
        sub.check_call(comm)

    def _prepare_subdirs(self):
        for i in range(self.njobs):
            subdir = '{}/job{}'.format(self.rundir,i)
            os.makedirs(subdir,exist_ok=True)
        return

if(__name__=='__main__'):
    main(sys.argv)