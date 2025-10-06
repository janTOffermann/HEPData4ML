import os, uuid
import subprocess as sub
import ROOT as rt
import numpy as np
from typing import Optional, List, Set

def ConcatenateRootTreeFiles(input_files: List[str], output_file: str, tree_name: str,
    branches_to_drop: Optional[List[str]] = None):

    if not input_files:
        raise ValueError("No input files provided")

    branches_to_drop = set(branches_to_drop) if branches_to_drop else set()

    # Determine the branch structure from the first file.
    first_file = rt.TFile.Open(input_files[0], "READ")
    if not first_file or first_file.IsZombie():
        raise IOError(f"Cannot open file: {input_files[0]}")

    first_tree = first_file.Get(tree_name)
    if not first_tree:
        raise ValueError(f"Tree '{tree_name}' not found in {input_files[0]}")

    # Get all branches and filter out the ones to drop
    all_branches = [b.GetName() for b in first_tree.GetListOfBranches()]
    branches_to_keep = FilterBranches(all_branches, branches_to_drop)
    first_file.Close()

    # Create output file and tree
    out_file = rt.TFile.Open(output_file, "RECREATE")
    # out_file.SetCompressionAlgorithm(compression)
    # out_file.SetCompressionLevel(compression_level)

    # Create TChain to read all input files
    chain = rt.TChain(tree_name)
    for f in input_files:
        chain.Add(f)

    # Disable all branches, then enable only the ones we want
    chain.SetBranchStatus("*", 0)
    for branch in branches_to_keep:
        chain.SetBranchStatus(branch, 1)

    # Use CloneTree to copy only active branches
    out_tree = chain.CloneTree(-1, "fast")

    out_file.cd()
    out_tree.Write()
    out_file.Close()

def FilterBranches(all_branches: List[str], patterns_to_drop: Set[str]) -> List[str]:
    import fnmatch

    branches_to_keep = []
    for branch in all_branches:
        should_drop = False
        for pattern in patterns_to_drop:
            if fnmatch.fnmatch(branch, pattern):
                should_drop = True
                break
        if not should_drop:
            branches_to_keep.append(branch)

    return branches_to_keep

def AddEventIndices(input_file:str,  tree_name:str, cwd=None, output_file: Optional[str]=None,index_branch_name:str="Event.Index", offset:int=0):

    if(cwd is not None):
        input_file = '{}/{}'.format(cwd,input_file)

    modify_in_place = output_file is None
    if modify_in_place:
        output_file = input_file.replace('.root','_tmp_{}.root'.format(str(uuid.uuid4())))

    # Open input file
    in_file = rt.TFile.Open(input_file, "READ")
    if not in_file or in_file.IsZombie():
        raise IOError("Cannot open file: {}".format(input_file))

    in_tree = in_file.Get(tree_name)
    if not in_tree:
        raise ValueError("Tree '{}' not found in {}".format(tree_name,input_file))

    n_entries = in_tree.GetEntries()

    # Check if branch already exists
    if in_tree.GetBranch(index_branch_name):
        print("Warning: Branch '{}' already exists in file {}.".format(index_branch_name,input_file))
        in_file.Close()
        return

    # Create output file with same compression settings
    out_file = rt.TFile.Open(output_file, "RECREATE")
    out_file.SetCompressionAlgorithm(in_file.GetCompressionAlgorithm())
    out_file.SetCompressionLevel(in_file.GetCompressionLevel())

    # Clone the tree structure (no entries yet)
    out_tree = in_tree.CloneTree(0)

    # Create array for the new branch (Python array for ROOT compatibility)
    index_value = np.zeros(1,dtype=np.dtype('uint64'))

    # Add the new branch
    new_branch = out_tree.Branch(
        index_branch_name,
        index_value,
        "{}/l".format(index_branch_name) # ULong64_t
    )

    for i in range(n_entries):
        in_tree.GetEntry(i)
        index_value[0] = i + offset
        out_tree.Fill()

    # Write and close
    out_file.cd()
    out_tree.Write()
    out_file.Close()
    in_file.Close()

    # If modifying in place, replace original file
    if modify_in_place:
        sub.check_call(['mv',output_file,input_file])
    return