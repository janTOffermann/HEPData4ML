# For loading our custom VectorCalcs library.
import ROOT as rt
import sys,os,glob,pathlib
import subprocess as sub
import util.qol_utils.qol_util as qu

def BuildVectorCalcs(force=False):

    # Check if the VectorCalcs library is already built, by trying to load it.
    # We will disable printouts for this check, since it's okay if it fails.
    if(not force):
        with qu.stdout_redirected():
            try:
                LoadVectorCalcs()
                return
            except: pass

    script_dir = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root')
    executable = '/bin/bash'
    script = 'build.sh'
    comm = ['.',script]
    env = os.environ.copy()
    print('Building VectorCalcs library.')
    sub.check_call(comm,shell=False,cwd=script_dir,env=os.environ.copy(),executable=executable,stderr=sub.DEVNULL,stdout=sub.DEVNULL)

def LoadVectorCalcs():
    # Load our custom ROOT library.
    try:
        a = rt.VectorCalcs
        return
    except: pass

    # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
    # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
    # TODO: Is there a more elegant workaround for this?
    lib_extensions = ['so','dylib']
    custom_lib_paths = [os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root/vectorcalcs/build/lib/libVectorCalcs.{}'.format(x)) for x in lib_extensions]
    custom_inc_paths = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/root/vectorcalcs/build/include/vectorcalcs')
    custom_inc_paths = glob.glob(custom_inc_paths + '/**/*.h',recursive=True)

    # Check for any of the libraries.
    found_libary = False
    for libpath in custom_lib_paths:
        if(pathlib.Path(libpath).exists()):
            custom_lib_path = libpath
            found_libary = True
            break

    if(not found_libary):
        print('Error: The VectorCalcs lib has not been built!')
        assert False

    for inc_path in custom_inc_paths:
        command = '#include "{}"'.format(inc_path)
        status = rt.gInterpreter.Declare(command)
        assert status, 'The following header file did not load properly: {}'.format(inc_path)

    status = rt.gSystem.Load(custom_lib_path)
    assert status == 0, 'The VectorCalcs lib did not load properly.'
    return
