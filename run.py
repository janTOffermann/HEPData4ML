import sys,os
import argparse as ap
import subprocess as sub
from util.generation import Generate, GenerateSimple, CopyTruth
from util.delphes import BuildDelphes, HepMC3ToDelphes
from util.conversion import DelphesWithTruthToHDF5, RemoveFailedFromHDF5, SplitHDF5

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-n', '--nevents', type=int, help='Number of events per pt bin.', required=True)
    parser.add_argument('-p', '--ptbins', type=int, nargs='+', help='Transverse momentum bin edges.', required=True)
    parser.add_argument('-o', '--outfile', type=str, help='Output HDF5 file name.', default='events.h5')
    parser.add_argument('-O', '--outdir', type=str, help='Output directory.', default=None)
    parser.add_argument('-g', '--generation',type=int, help='Whether or not to do event generation.', default=True)
    args = vars(parser.parse_args())

    nevents_per_bin = args['nevents']
    pt_bin_edges = args['ptbins']
    h5_file = args['outfile']
    outdir = args['outdir']
    do_generation =args['generation'] > 0 # TODO: Find nicer way to handle Boolean -- argparse is weird.
    nbins = len(pt_bin_edges) - 1

    # TODO: Make Delphes optional.
    use_delphes = True

    # Keep track of some files we create.
    jet_files = []
    truth_files = []

    # Prepare the output directory.
    if(outdir is None): outdir = os.getcwd()
    else: os.makedirs(outdir,exist_ok=True)
    h5_file = '{}/{}'.format(outdir,h5_file)

    for i in range(nbins):

        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.

        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]
        hepfile = '{}/events_{}-{}.hepmc'.format(outdir,pt_min,pt_max)
        if(do_generation): Generate(nevents_per_bin, pt_min, pt_max, hepfile)

        if(use_delphes):
            # Pass the HepMC file to Delphes. Will output a ROOT file.
            delphes_dir = BuildDelphes() # build Delphes if it does not yet exist
            delphesfile = HepMC3ToDelphes(hepmc_file = hepfile, delphes_dir = delphes_dir)
            jet_files.append(delphesfile)

        # Extract the truth-level particles from the full HepMC file.
        truthfile = hepfile.replace('.hepmc','_truth.hepmc')
        truthfile = CopyTruth(hepfile, truthfile)
        truth_files.append(truthfile)

        # Compress the full HepMC file. Will use tar (shipped with the custom conda env).
        delete_hepmc = True
        compress_file = hepfile.replace('.hepmc','.tar.bz2')
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],hepfile.split('/')[-1]]
        if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)

    # Now put everything into an HDF5 file.
    DelphesWithTruthToHDF5(
        delphes_files=jet_files,
        truth_files=truth_files,
        h5_file=h5_file
    )

    # Cleanup: Delete the jet files, they can in principle be recreated from the compressed HepMC files.
    comm = ['rm'] + jet_files
    sub.check_call(comm)

    # Now also compress the truth files.
    for truth_file in truth_files:
        compress_file = truth_file.replace('.hepmc','.tar.bz2')
        cwd = '/'.join(compress_file.split('/')[:-1])
        comm = ['tar','-cjf',compress_file.split('/')[-1],truth_file.split('/')[-1]]
        if(delete_hepmc): comm.append('--remove-files')
        sub.check_call(comm,shell=False,cwd=cwd)

    # Remove any failed events (e.g. detector-level events with no jets passing cuts).
    RemoveFailedFromHDF5(h5_file)

    # Now split the HDF5 file into training, testing and validation samples.
    split_ratio = (7,2,1)
    SplitHDF5(h5_file, split_ratio)

    # Optionally delete the full HDF5 file.
    delete_h5 = False
    if(delete_h5):
        comm = ['rm',h5_file]
        sub.check_call(comm)

if __name__ == '__main__':
    main(sys.argv)

