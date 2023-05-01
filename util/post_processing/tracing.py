# Credit to Tim Hoffman (https://github.com/hoffmantj) for implementing
# the tracing algorithm ( Tracer.Process() ) in this code.

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
import subprocess as sub
import argparse as ap
import util.qol_utils.qol_util as qu

class Tracer:
    def __init__(self,verbose=False):
        self.status = False
        self.SetVerbosity(verbose)
        self.energy_ratio_truth = None
        self.energy_ratio_smeared = None
        self.configurator = None

        self.print_prefix = '\n\tTracer'
        self.progress_bar_length = 50
        self.progress_bar_prefix = '\tTracing daughter particles:'
        self.progress_bar_suffix = 'Complete'
        self.SetKeyPrefix('energy_ratio')

    def SetKeyPrefix(self,val):
        self.key_prefix = val

    def SetH5EventFile(self,file):
        self.h5_file = file

    def SetDelphesFile(self,file):
        self.delphes_file = file

    def SetIndicesFile(self,file):
        self.indices_file = file

    def SetVerbosity(self,flag):
        self.verbose = flag

    def GetEnergyRatioTruth(self):
        return self.energy_ratio_truth

    def GetEnergyRatioSmeared(self):
        return self.energy_ratio_smeared

    def RequiresIndexing(self):
        return True

    def SetConfigurator(self,configurator):
        self.configurator = configurator

    def Initialize(self): # TODO: May want to consider chunking things and using a buffer? Memory usage will scale better for larger files.
        """
        Reads in files, and places the relevant information in memory.
        """

        # Check if Delphes generation is on. If not, we will presumably run into some issues because Delphes files are not produced,
        # so we should disable this post-processing.
        if(not self.configurator.GetDelphesConfig()):
            print('Error: Delphes option is turned off in config. Skipping Tracer.')
            self.status = False
            return

        if(self.h5_file is None):
            self.status = False
            return

        if(self.delphes_file is None):
            self.status = False
            return

        if(self.indices_file is None):
            self.status = False
            return

        # Fetch the (HepMC) indices of particles in the final state that are also in the truth selection.
        # These use 1-indexing -- as do the Delphes towers' particle references in Tower/Tower.Particles.
        f_indices = h5.File(self.indices_file,"r")
        self.daughter_indices_fs = f_indices['indices'][:,:,0] # w.r.t. the final-state HepMC file & Delphes file
        self.daughter_indices_truth = f_indices['indices'][:,:,1] # w.r.t. the truth HepMC file
        f_indices.close()

        # Fetch the tower energy and tower particle indices from Delphes.
        # The tower particle indices use 1-indexing.
        f_delphes = ur.open(self.delphes_file)
        t_delphes = f_delphes['Delphes']
        self.indices_delphes = t_delphes['Tower/Tower.Particles'].array()['refs']
        self.tower_energies_delphes = t_delphes['Tower/Tower.E'].array()
        self.particle_energies_delphes = t_delphes['Particle/Particle.E'].array()
        f_delphes.close()

        # Fetch the final HDF5 file, that contains selected jets' constituent four-momenta.
        f_event = h5.File(self.h5_file,'r')
        self.nobj = f_event['Nobj'][:] # number of constituents in each jet
        self.nobj_max = f_event['Pmu'].shape[1]
        self.event_indices = f_event['event_idx'][:] # event index, since some events may be dropped in jet clustering. Uses 0-indexing.
        self.selected_tower_indices = f_event['final_state_idx'][:] # Index of Delphes towers in the jet. Uses 0-indexing.
        self.truth_Pmu = f_event['truth_Pmu'][:]
        self.nobj_truth = f_event['truth_Nobj'][:]
        f_event.close()
        self.nevents = self.nobj.shape[0]

        self.status = True
        return

    def Process(self):
        self.Initialize()
        if(not self.status):
            print('Error: Tracer not properly configured.')
            return

        if(self.verbose):
            print('{}.Process: Number of events = {}'.format(self.print_prefix,self.nevents))
            qu.printProgressBarColor(0,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)

        daughter_energy    = np.zeros((self.nevents, self.nobj_max))
        tower_energy       = np.full( (self.nevents, self.nobj_max), np.nan)
        truth_tower_energy = np.full( (self.nevents, self.nobj_max), np.nan)

        #true daughter to smeared total
        for i, event_idx in enumerate(self.event_indices): # Looping over events
            nobj = self.nobj[i]
            for j, tower_idx in enumerate(self.selected_tower_indices[i,:nobj]): # Looping over towers in event i
                daughter_energy[i,j] = 0 #set tower energy =0 when you touch said tower
                tower_energy[i,j] = self.tower_energies_delphes[event_idx,tower_idx]
                residual = [i for n, i in enumerate(self.indices_delphes[event_idx,tower_idx]) if i not in self.indices_delphes[event_idx,tower_idx][:n]]# NOTE: eliminate Delphes's apparent double-counting of anti-muons
                for part_idx in residual:
                    if part_idx in self.daughter_indices_fs[event_idx]:
                        # part_idx is the daughter particle's index w.r.t. the listing of final-state particles that went into Delphes.
                        # We want its index w.r.t. the listing of truth-level particles we saved. These are not only in the truth HepMC file,
                        # but are also already saved to the final HDF5 file in the "truth_Pmu" array.
                        part_idx_truth = self.daughter_indices_truth[event_idx, np.where(self.daughter_indices_fs[event_idx] == part_idx)[0][0]]
                        part_energy = self.truth_Pmu[i,part_idx_truth-1,0] # conversion from 1-indexing to 0-indexing on the 2nd axis
                        daughter_energy[i,j] += part_energy

                    # Get this particle's truth-level energy, to add to the "truth tower energy".
                    part_truth_energy = self.particle_energies_delphes[event_idx, part_idx-1] # part_idx uses 1-indexing
                    if(np.isnan(truth_tower_energy[i,j])): truth_tower_energy[i,j] = part_truth_energy
                    else: truth_tower_energy[i,j] += part_truth_energy

            if(self.verbose): qu.printProgressBarColor(i+1,self.nevents,prefix=self.progress_bar_prefix,suffix=self.progress_bar_suffix,length=self.progress_bar_length)
        self.energy_ratio_smeared = daughter_energy/tower_energy
        self.energy_ratio_truth = daughter_energy/truth_tower_energy
        return

    def __call__(self, delphes_file,h5_file,indices_file,output_file=None,verbose=None, copts=9, key=None):
        self.SetH5EventFile(h5_file)
        self.SetDelphesFile(delphes_file)
        self.SetIndicesFile(indices_file)
        if(verbose is not None): self.SetVerbosity(verbose)
        if(key is not None):
            self.key_prefix = key # TODO: Is this clear? Maybe should rename argument, but I want it to be a generic function signature
        self.Process()

        key_truth = '{}_truth'.format(self.key_prefix)
        key_smeared = '{}_smeared'.format(self.key_prefix)

        if(self.status):
            if(output_file is None): output_file = h5_file
            else: sub.check_call(['cp',h5_file,output_file])
            f = h5.File(output_file,'a')
            if(self.verbose): print('\tWriting {} to {}.'.format(key_truth,output_file))
            d = f.create_dataset(key_truth,data=self.GetEnergyRatioTruth(),compression='gzip',compression_opts=copts)
            if(self.verbose): print('\tWriting {} to {}.'.format(key_smeared,output_file))
            d = f.create_dataset(key_smeared,data=self.GetEnergyRatioSmeared(),compression='gzip',compression_opts=copts)
            f.close()
            return output_file
        return h5_file

def Process(h5_file,delphes_file,indices_file,output_file=None,verbose=False):
    tracer = Tracer(verbose)
    return tracer(delphes_file,h5_file,indices_file,output_file,verbose)

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