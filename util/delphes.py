import sys, os, glob
import subprocess as sub

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
        sub.check_call(['wget', 'http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz'],
                       shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        sub.check_call(['tar', '-zxf', 'Delphes-3.5.0.tar.gz'],
                       shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        sub.check_call(['rm', 'Delphes-3.5.0.tar.gz'],
                       shell=False,cwd=delphes_dir, stdout=f, stderr=g)
        # Now make.
        print('Making Delphes.')
        sub.check_call(['make', '-j{}'.format(j)],
                       shell=False, cwd = '{}/Delphes-3.5.0'.format(delphes_dir), stdout=f, stderr=g)
    return delphes_dir
        
def HepMC3ToDelphes(hepmc_file, output_file=None, delphes_card=None, delphes_dir=None, logfile=None):
    if(delphes_dir is None):
        delphes_dir = os.path.dirname(os.path.abspath(__file__)) + '/../delphes'
        
    delphes_ex = glob.glob('{}/**/DelphesHepMC3'.format(delphes_dir),recursive=True)
    if(len(delphes_ex) == 0):
        print('Error: Delphes executable not found at {}/**?DelphesHepMC3 .'.format(delphes_dir))
        assert(False)
    delphes_ex = delphes_ex[0]
    
    # default to the ATLAS Delphes card
    if(delphes_card is None): delphes_card = glob.glob('{}/**/delphes_card_ATLAS.tcl'.format(delphes_dir),recursive=True)[0]
        
    # default to using HepMC filename for ROOT output
    if(output_file is None): output_file = hepmc_file.replace('.hepmc','.root')
        
    # Delphes will crash if output file already exists, so we need to remove it (TODO: careful about overwriting things).
    try: os.remove(output_file)
    except: pass
        
    if(logfile is not None):
        with open(logfile,'w') as f:
            sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                           shell=False, stdout=f, stderr=f)            
        
    else:
        sub.check_call([delphes_ex, delphes_card, output_file, hepmc_file],
                       shell=False, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
        
    #comm = '{ex} {card} {outfile} {infile}'.format(ex=delphes_ex, card=delphes_card, outfile=output_file, infile=hepmc_file)
    #sub.check_call(comm,shell=True)
    
    return output_file # return the name (esp. useful if none was provided)
        
        
    
    
    
    
    