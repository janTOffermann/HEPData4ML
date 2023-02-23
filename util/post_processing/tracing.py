# The purpose of this code is to help with tracing particles through Delphes.
# Specifically, this is to be used for stable particles that are selected by
# the (user-specified) truth particle selection algorithm.
# One can configure the main script (run.py) to save files containing the indices
# of particles that are in both the "truth" particle list and the final state list,
# and these indices can be referenced when looking at Delphes Tower objects to determine
# if these particles hit a particular Delphes calorimeter tower. For example, this is useful
# if the truth selector is selecting stable decay products of a W boson from top decay, and we
# want to determine which tower(s) that W boson's stable daughters interacted with.

import sys
import numpy as np
import h5py as h5
import uproot as ur
import pyhepmc as hep
import subprocess as sub
import argparse as ap

class Tracer:
    def __init__(self):
        self.status = False
        self.SetVerbosity(False)
        self.energy_ratio_truth = None
        return

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetDelphesFile(self,file):
        self.delphes_file = file

    def SetIndicesFile(self,file):
        self.indices_file = file

    def SetTruthHepMCFile(self,file):
        self.truth_file = file

    def SetVerbosity(self,flag):
        self.verbose = flag

    def GetEnergyRatioTruth(self):
        return self.energy_ratio_truth

    # Reads in files, places relevant information in memory.
    def Initialize(self):
        if(self.h5_file is None):
            self.status = False
            return

        if(self.delphes_file is None):
            self.status = False
            return

        if(self.indices_file is None):
            self.status = False
            return

        if(self.truth_file is None):
            self.status = False
            return

        # Fetch the (HepMC) indices of particles in the final state that are also in the truth selection.
        f_indices = h5.File(self.indices_file,"r")
        self.daughter_indices = f_indices['indices'][:]
        f_indices.close()

        # Fetch the Delphes file.
        f_delphes = ur.open(self.delphes_file)
        t_delphes = f_delphes['Delphes']
        self.indices_delphes = t_delphes['Tower/Tower.Particles'].array()
        self.energies_delphes = t_delphes['Tower/Tower.E'].array()
        f_delphes.close()

        # Fetch the final HDF5 file, that contains selected jets' constituent four-momenta.
        f_event = h5.File(self.h5_file,'r')
        self.nobj = f_event['Nobj'][:]
        self.nobj_max = f_event['Pmu'].shape[1]
        self.event_indices = f_event['event_idx'][:]
        self.selected_tower_indices = f_event['final_state_idx'][:]
        f_event.close()
        self.nevents = self.nobj.shape[0]

        # Fetch the truth HepMC file.
        self.truth_events = []
        with hep.io.ReaderAscii(self.truth_file) as f_truth:
            for ev in f_truth:
                self.truth_events.append(ev)

        self.status = True
        return

    def Process(self):
        self.Initialize()
        if(not self.status):
            print('Error: Tracer not properly configured.')
            return

        if(self.verbose): print('Process: Number of events = {} .'.format(self.nevents))
        energy_ratio_truth = np.zeros((self.nevents,self.nobj_max))

        for i in range(self.nevents):
            unfiltered_event_idx = self.event_indices[i] # there may be more events in Delphes file than the final h5 file, since some Delphes outputs may have failed jet cuts and those events got discarded.
            daughter_idxs = np.trim_zeros(self.daughter_indices[unfiltered_event_idx,:,0]) # Note: These indices are using 1-indexing! This is how indexing is done within HepMC3.
            daughter_idxs_truth = np.trim_zeros(self.daughter_indices[unfiltered_event_idx,:,1]) # Same as above, but w.r.t. the truth HepMC3 file.

            # Get the list of selected towers and their energies.
            towers = self.indices_delphes[unfiltered_event_idx] # all towers for event i, not just those we selected via jet clustering!
            selected_towers = towers[ self.selected_tower_indices[i,:self.nobj[i]]]
            tower_energies = self.energies_delphes[unfiltered_event_idx]
            selected_tower_energies = tower_energies[ self.selected_tower_indices[i,:self.nobj[i]]]

            # Fetch the daughter particles' energies from the truth HepMC3 file.
            truth_event = self.truth_events[unfiltered_event_idx]
            truth_particles = truth_event.particles
            truth_energies = np.array([x.momentum.e for x in truth_particles],dtype=float)
            selected_truth_energies = truth_energies[daughter_idxs_truth - 1] # daughter_idxs_truth is 1-indexed, but here we need to use 0-indexing based on how this array has been made

            for j,tower in enumerate(selected_towers): # Loop over towers that are in the selected jet
                daughter_energy_sum = 0.
                contained_particles = np.array(tower['refs'],dtype=int) - 1 # TODO: Should triple-check the indexing style (0 vs. 1).
                for k,idx in enumerate(contained_particles): # Loop over particles contained by this tower
                    if(idx in daughter_idxs):
                        idx2 = np.where(daughter_idxs == idx)[0][0] # sort of inelegant, but this should always work by construction
                        daughter_energy_sum += selected_truth_energies[idx2]
                energy_ratio_truth[i,j] = daughter_energy_sum / selected_tower_energies[j]
        self.energy_ratio_truth = energy_ratio_truth
        return

def Process(h5_file,delphes_file,indices_file,truth_file,output_file=None,verbose=False):
    tracer = Tracer()
    tracer.SetH5EventFile(h5_file)
    tracer.SetDelphesFile(delphes_file)
    tracer.SetIndicesFile(indices_file)
    tracer.SetTruthHepMCFile(truth_file)
    tracer.SetVerbosity(verbose)
    tracer.Process()

    # Get the ratio of true daughter particle energy, versus the total energy deposited in the Delphes Tower.
    # Since the numerator is total true energy, and not deposited/measured energy, this ratio can be larger than 1!
    energy_ratio_truth = tracer.GetEnergyRatioTruth()

    replace_h5 = False
    if(output_file is None):
        output_file = h5_file.replace('.h5','_proc.h5')

    if(output_file == h5_file):
        replace = True

    if(not replace):
        sub.check_call(['cp',h5_file,output_file])

    f = h5.File(output_file,'a')
    key = 'energy_ratio_truth'
    copts = 7
    d = f.create_dataset(key,data=energy_ratio_truth,compression='gzip',compression_opts=copts)
    f.close()

    return output_file

# Main function lets this be run as an executable, which may be a useful feature to have.
# Ideally this is set up so that it can be automatically run by the main generation script (run.py).
def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-h5','--h5_file',type=str,required=True,help='Input HDF5 file containing *all* events from generation.')
    parser.add_argument('-d','--delphes_file',type=str,required=True,help='Input Delphes file (run.py can be used with "-cd 0" option to turn off deletion of Delphes files).')
    parser.add_argument('-i','--index_file',type=str,required=True,help='Input HDF5 file containing indices of particles in both final-state and truth selections (creation of this file set in config).')
    parser.add_argument('-t','--truth_file',type=str,required=True,help='Input HepMC file corresponding to truth-level particle selection.')
    parser.add_argument('-v','--verbose',type=int,default=0)
    args = vars(parser.parse_args())

    h5_file = args['h5_file']
    delphes_file = args['delphes_file']
    indices_file = args['index_file']
    truth_file = args['truth_file']
    verbose = args['verbose'] > 0

    Process(h5_file,delphes_file,indices_file,truth_file,verbose=verbose)
    return

if(__name__ == '__main__'):
    main(sys.argv)