import os, glob, pathlib
import subprocess as sub
import util.qol_utils.qol_util as qu

def BuildDelphes(delphes_dir=None, j=4, force=False, verbose=False):
    if(delphes_dir is None):
        delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../delphes'

    # Check if Delphes is already built at destination.
    if(not force):
        ex = glob.glob('{}/**/DelphesHepMC3'.format(delphes_dir),recursive=True)
        if(len(ex) > 0):
            if(verbose): print('Found DelphesHepMC3 @ {}'.format(ex[0]))
            return

    # Make the Delphes dir if it does not exist.
    try: os.makedirs(delphes_dir)
    except: pass

    # Put output into log files.
    logfile = '{}/log.stdout'.format(delphes_dir)
    errfile = '{}/log.stderr'.format(delphes_dir)

    with open(logfile,'w') as f, open(errfile,'w') as g:

        # Fetch Delphes source.
        print('Downloading Delphes.')
        delphes_download = 'https://github.com/delphes/delphes/archive/refs/tags/3.5.1pre01.tar.gz'
        delphes_file = 'delphes-{}'.format(delphes_download.split('/')[-1]) # TODO: A bit hacky
        # Depending on Linux/macOS, we use wget or curl.
        has_wget = True
        with qu.stdout_redirected():
            try: sub.check_call('which wget'.split(' '))
            except:
                has_wget = False
                pass

        if(has_wget):
            sub.check_call('wget --no-check-certificate --content-disposition {} -O {}'.format(delphes_download,delphes_file).split(' '), shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        else:
            sub.check_call('curl -LJ {} -o {}'.format(delphes_download,delphes_file).split(' '), shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        sub.check_call(['tar', '-zxf', delphes_file], shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        sub.check_call(['rm', delphes_file], shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        # Now make.
        print('Making Delphes.')
        sub.check_call(['make', '-j{}'.format(j)],
                       shell=False, cwd = '{}/{}'.format(delphes_dir,delphes_file.replace('.tar.gz','')), stdout=f, stderr=g)
    return delphes_dir

def HepMC3ToDelphes(hepmc_file, output_file=None, delphes_card=None, delphes_dir=None, logfile=None, cwd=None, force=False):
    hepmc_file_nodir = hepmc_file
    if(cwd is not None):
        hepmc_file = '{}/{}'.format(cwd,hepmc_file)

    # default to using HepMC filename for ROOT output
    if(output_file is None):
        output_file = hepmc_file.replace('.hepmc','.root')
        if(cwd is not None): output_file_nodir = hepmc_file_nodir.replace('.hepmc','.root')

    # check if the output file already exists (i.e. look for file with same name)
    if(pathlib.Path(output_file).exists() and not force):
        print('\t\tDelphes ROOT file {} already found, skipping its generation.'.format(output_file))
        if(cwd is not None): return output_file_nodir
        return output_file

    if(delphes_dir is None):
        delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../delphes'

    delphes_ex = glob.glob('{}/**/DelphesHepMC3'.format(delphes_dir),recursive=True)
    if(len(delphes_ex) == 0):
        print('Error: Delphes executable not found at {}/**?DelphesHepMC3 .'.format(delphes_dir))
        assert(False)
    delphes_ex = delphes_ex[0]

    # default to the ATLAS Delphes card
    if(delphes_card is None): delphes_card = glob.glob('{}/**/delphes_card_ATLAS.tcl'.format(delphes_dir),recursive=True)[0]

    # Delphes will crash if output file already exists, so we need to remove it.
    try: os.remove(output_file)
    except: pass

    if(logfile is not None):
        with open(logfile,'w') as f:
            sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                           shell=False, stdout=f, stderr=f)

    else:
        sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                       shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)

    if(cwd is not None): return output_file_nodir
    return output_file # return the name (esp. useful if none was provided)