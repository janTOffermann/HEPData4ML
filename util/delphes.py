import os, glob, pathlib
import subprocess as sub
from util.qol_utils.misc import stdout_redirected

class DelphesWrapper:
    """
    A class for running (and downloading/building, if necessary) the Delphes
    fast detector simulation framework.
    """
    def __init__(self,delphes_dir=None):
        self.SetDirectory(delphes_dir)
        self.delphes_download = 'https://github.com/delphes/delphes/archive/refs/tags/3.5.1pre12.tar.gz' # LCG_105 - LCG_107 have 3.5.1pre09
        self.delphes_file = 'delphes-{}'.format(self.delphes_download.split('/')[-1]) # TODO: A bit hacky
        self.executable = None
        self.default_card = "delphes_card_ATLAS.tcl" # Default to the ATLAS card. Will search for this within the Delphes directory.

    def SetDirectory(self,delphes_dir=None):
        self.delphes_dir = delphes_dir
        if(self.delphes_dir is None):
            self.delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../delphes'
        self.delphes_dir = os.path.normpath(self.delphes_dir)

        # Make the Delphes dir if it does not exist.
        os.makedirs(self.delphes_dir,exist_ok=True)

        # Files that are used for logging progress with downloading/building Delphes.
        self.logfile = '{}/log.stdout'.format(self.delphes_dir)
        self.errfile = '{}/log.stderr'.format(self.delphes_dir)

    def GetDirectory(self):
        return self.delphes_dir

    def GetExecutable(self):
        return self.executable

    def PrepDelphes(self, j=4, force=False, verbose=False):
        """
        Checks if Delphes is built. If it's not, it will download
        and build Delphes.
        """
        # Check if Delphes is already built at destination.
        if(not force):
            ex = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)
            if(len(ex) > 0):
                self.executable = ex[0]
                if(verbose): print('Found DelphesHepMC3 @ {}'.format(self.executable))
                return

        self.DownloadDelphes() # TODO: Could check to see if source code is downloaded, though it's hard to make a perfect check.
        self.BuildDelphes(j=j)
        self.executable = ex = glob.glob('{}/**/DelphesHepMC3'.format(self.delphes_dir),recursive=True)[0]
        return

    def DownloadDelphes(self):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            # Fetch Delphes source.
            print('Downloading Delphes.')
            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget):
                sub.check_call('wget --no-check-certificate --content-disposition {} -O {}'.format(self.delphes_download,self.delphes_file).split(' '), shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            else:
                sub.check_call('curl -LJ {} -o {}'.format(self.delphes_download,self.delphes_file).split(' '), shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['tar', '-zxf', self.delphes_file], shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
            sub.check_call(['rm', self.delphes_file], shell=False,cwd=self.delphes_dir, stdout=f, stderr=g)
        return

    def BuildDelphes(self, j=4):
        with open(self.logfile,'w') as f, open(self.errfile,'w') as g:
            print('Making Delphes.')
            sub.check_call(['make', '-j{}'.format(j)],
                        shell=False, cwd = '{}/{}'.format(self.delphes_dir,self.delphes_file.replace('.tar.gz','')), stdout=f, stderr=g)
        return

    def HepMC3ToDelphes(self,hepmc_file, output_file=None, delphes_card=None, logfile=None, cwd=None, force=False):
        hepmc_file_nodir = hepmc_file
        if(cwd is not None):
            hepmc_file = '{}/{}'.format(cwd,hepmc_file)

        if(output_file is not None and cwd is not None):
            output_file = '{}/{}'.format(cwd,output_file)

        # default to using HepMC filename for ROOT output
        if(output_file is None):
            output_file = hepmc_file.replace('.hepmc','.root')

        if(cwd is not None): output_file_nodir = hepmc_file_nodir.replace('.hepmc','.root')
        else: output_file_nodir = output_file.split('/')[-1]

        # check if the output file already exists (i.e. look for file with same name)
        if(pathlib.Path(output_file).exists() and not force):
            print('\t\tDelphes ROOT file {} already found, skipping its generation.'.format(output_file))
            if(cwd is not None): return output_file_nodir
            return output_file

        if(delphes_card is None): delphes_card = glob.glob('{}/**/{}'.format(self.delphes_dir,self.default_card),recursive=True)[0]

        # Delphes will crash if output file already exists, so we need to remove it.
        try: os.remove(output_file)
        except: pass
        command = [self.executable, delphes_card, output_file, hepmc_file]
        if(logfile is not None):
            with open(logfile,'w') as f:
                sub.check_call(command, stdout=f, stderr=f)

        else:
            sub.check_call(command, stdout=sub.DEVNULL, stderr=sub.DEVNULL)

        if(cwd is not None): return output_file_nodir
        return output_file # return the name (esp. useful if none was provided)








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
        with stdout_redirected():
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

    if(output_file is not None and cwd is not None):
        output_file = '{}/{}'.format(cwd,output_file)

    # default to using HepMC filename for ROOT output
    if(output_file is None):
        output_file = hepmc_file.replace('.hepmc','.root')

    if(cwd is not None): output_file_nodir = hepmc_file_nodir.replace('.hepmc','.root')
    else: output_file_nodir = output_file.split('/')[-1]

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