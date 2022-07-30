import sys,os
import argparse as ap
import subprocess as sub
from util.generation import Generator
from util.delphes import BuildDelphes, HepMC3ToDelphes
from util.conversion import Processor, RemoveFailedFromHDF5, SplitHDF5
from util.config import GetDelphesConfig

def CompressHepMC(files, delete=True):
    for file in files:
        compress_file = file.replace('.hepmc','.tar.bz2')
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],file.split('/')[-1]]
        # if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)
        if(delete):
            sub.check_call(['rm',file.split('/')[-1]],shell=False,cwd=cwd)
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n', '--nevents', type=int, help='Number of events per pt bin.', required=True)
    parser.add_argument('-p', '--ptbins', type=int, nargs='+', help='Transverse momentum bin edges.', required=True)
    parser.add_argument('-o', '--outfile', type=str, help='Output HDF5 file name.', default='events.h5')
    parser.add_argument('-O', '--outdir', type=str, help='Output directory.', default=None)
    parser.add_argument('-g', '--generation',type=int, help='Whether or not to do event generation.', default=True)
    parser.add_argument('-s', '--sep_truth',type=int, help='Whether or not to store truth-level particles in separate arrays.', default=True)

    args = vars(parser.parse_args())

    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
    h5_file = args['outfile']
    outdir = args['outdir']
    do_generation =args['generation'] > 0 # TODO: Find nicer way to handle Boolean -- argparse is weird.
    separate_truth_particles = args['sep_truth'] > 0
    nbins = len(pt_bin_edges) - 1

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

    # Prepare the output directory.
    if(outdir is None): outdir = os.getcwd()
    else: os.makedirs(outdir,exist_ok=True)
    h5_file = '{}/{}'.format(outdir,h5_file)

    # # Prepare our custom ROOT library, that is used to do coordinate conversions and other 4-vector calcs.
    # BuildVectorCalcs(force=True) # make sure we re-build the library
    # LoadVectorCalcs()

    for i in range(nbins):
        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.

        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]
        hep_file = '{}/events_{}-{}.hepmc'.format(outdir,pt_min,pt_max)
        if(do_generation):
            generator = Generator(pt_min,pt_max)
            generator.SetOutputDirectory(outdir) # TODO: Currently only used for histograms, should implement so it is tacked onto hep_file etc.
            generator.Generate(nevents_per_bin,hep_file)

            # Generate(nevents_per_bin, pt_min, pt_max, hep_file)

        # Extract the truth-level particles from the full HepMC file.
        truthfile = hep_file.replace('.hepmc','_truth.hepmc') # TODO: Should be returned by Generate()
        truth_files.append(truthfile)

        if(use_delphes): # Case 1: Using Delphes
            # Pass the HepMC file to Delphes. Will output a ROOT file.
            delphes_dir = BuildDelphes() # build Delphes if it does not yet exist
            delphesfile = HepMC3ToDelphes(hepmc_file = hep_file, delphes_dir = delphes_dir)
            jet_files.append(delphesfile)

            # We can now compress the HepMC file. Once it's compressed, we delete the original file to recover space.
            CompressHepMC([hep_file],True)

        else: # Case 2: No Delphes
            jet_files.append(hep_file)

    # Now put everything into an HDF5 file.
    processor = Processor(use_delphes)
    processor.Process(jet_files,truth_files,h5_file,verbosity=h5_conversion_verbosity,separate_truth_particles=separate_truth_particles)

    if(use_delphes):
        # Cleanup: Delete the jet files, since they can always be recreated from the compressed HepMC files.
        comm = ['rm'] + jet_files
        sub.check_call(comm)

    else:
        #Cleanup: Compress the HepMC files.
        CompressHepMC(jet_files,True)

    # Cleanup.
    # Now also compress the truth files.
    for truth_file in truth_files:
        compress_file = truth_file.replace('.hepmc','.tar.bz2')
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],truth_file.split('/')[-1]]
        sub.check_call(comm,shell=False,cwd=cwd)
        sub.check_call(['rm',truth_file.split('/')[-1]],shell=False,cwd=cwd)

    # Remove any failed events (e.g. detector-level events with no jets passing cuts).
    RemoveFailedFromHDF5(h5_file)

    # Now split the HDF5 file into training, testing and validation samples.
    split_ratio = (7,2,1)
    SplitHDF5(h5_file, split_ratio)

    # Optionally delete the full HDF5 file.
    delete_h5 = False # TODO: Make this configurable
    if(delete_h5):
        comm = ['rm',h5_file]
        sub.check_call(comm)

if __name__ == '__main__':
    main(sys.argv)