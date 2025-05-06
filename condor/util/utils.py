import os,glob,pathlib,uuid
import subprocess as sub
# Use this file to place functions that will otherwise clutter up prep_condor

def PreparePayload(rundir,payload,gitdir='HEPData4ML'):
    """
    This function packages up the HEPData4ML repository into
    a tar.gz archive. This is to be used when one needs to
    ship the repository to an HTCondor job.
    """

    this_dir = os.path.dirname(os.path.abspath(__file__))

    tmp_dir = 'dir_{}'.format(uuid.uuid4())
    os.makedirs(tmp_dir)

    # Only fetch the things we really need to copy. For example, we don't copy this condor subdir!
    paths = [
        str(pathlib.Path('{}/../../run.py'.format(this_dir)).absolute()),
        str(pathlib.Path('{}/../../config'.format(this_dir)).absolute()),
        str(pathlib.Path('{}/../../util'.format(this_dir)).absolute()),
        str(pathlib.Path('{}/../../setup'.format(this_dir)).absolute())
    ]

    gitdir_full = '{}/{}'.format(tmp_dir,gitdir)
    os.makedirs(gitdir_full)
    for path in paths:
        command = ['cp','-r',path,'{}/'.format(gitdir_full)]
        sub.check_call(command)

    # cleanup
    cleanup =  glob.glob('{}/**/__pycache__'.format(gitdir_full),recursive=True)
    cleanup += glob.glob('{}/util/**/build'.format(gitdir_full),recursive=True)
    cleanup += glob.glob('{}/fastjet'.format(gitdir_full),recursive=True)
    cleanup += glob.glob('{}/delphes'.format(gitdir_full),recursive=True)
    for obj in cleanup:
        command = ['rm','-r',obj]
        sub.check_call(command)

    payload_full = '{}/{}'.format(tmp_dir,payload)
    command = ['tar','-czf',payload,gitdir]
    sub.check_call(command,cwd=tmp_dir)
    command = ['rm','-r',gitdir_full]
    sub.check_call(command)
    command = ['mv',payload_full,rundir]
    sub.check_call(command)
    command = ['rm','-r',tmp_dir]
    sub.check_call(command)
    return

def PreparePayloadFromClone(rundir,payload,gitdir='HEPData4ML',branch='main'):
    """
    This function packages up the HEPData4ML repository into
    a tar.gz archive. This is to be used when one needs to
    ship the repository to an HTCondor job.
    It first fetches the repository via git.
    """
    repo = 'git@github.com:janTOffermann/HEPData4ML.git'
    tmp_dir = 'dir_{}'.format(uuid.uuid4())
    os.makedirs(tmp_dir)
    command = ['git','clone','-b',branch,repo,gitdir]
    sub.check_call(command,cwd=tmp_dir)
    gitdir_full = '{}/{}'.format(tmp_dir,gitdir)
    command = ['tar','-czf',payload,gitdir]
    sub.check_call(command,cwd=tmp_dir)
    payload_full = '{}/{}'.format(tmp_dir,payload)
    command = ['rm','-rf',gitdir_full]
    sub.check_call(command)
    command = ['mv',payload_full,rundir]
    sub.check_call(command)
    command = ['rm','-r',tmp_dir]
    sub.check_call(command)
    return

def GetGitBranch():
    command = ['git','rev-parse','--abbrev-ref', 'HEAD']
    result = sub.check_output(command).decode('utf-8')
    result = result.replace('\n','')
    return result