import os, uuid, time, pathlib, json
import subprocess as sub
import numpy as np
import h5py as h5
import ROOT as rt
from typing import Union, List
from util.config.config import Configurator

class MetaDataHandler:

    def __init__(self,configurator:'Configurator'):
        self.metadata = {}
        self.print_prefix = 'MetaDataHandler: '
        self.configurator = configurator
        self.Initialize()

    def _init_citations(self):
        """
        This is where we initialize the citation for HEPData4ML itself!
        Using the @software format, which works nicely together with the
        custom BibLaTeX style here: https://github.com/janTOffermann/CustomNumericComp
        """
        citation = {'HEPData4ML':
            """
@software{Offermann:HEPData4ML,
    author = {Offermann, Jan Tuzli\\'c and Liu, Xiaoyang and Hoffman, Timothy},
    title = {\\texttt{HEPData4ML}},
    url = {https://github.com/janTOffermann/HEPData4ML},
    year = {2023}
}
        """
        }

        # Also adding the HepMC3 citation, as the package is used throughout.
        # TODO: Is there a better place to put this? Doesn't have a dedicated module like Pythia8/Delphes/Fastjet.
        citation['HepMC3'] = """
@article{Buckley:2019xhk,
    author = {Buckley, Andy and Ilten, Philip and Konstantinov, Dmitri and L{\\"o}nnblad, Leif and Monk, James and Pokorski, Witold and Przedzinski, Tomasz and Verbytskyi, Andrii},
    title = "{The HepMC3 event record library for Monte Carlo event generators}",
    eprint = "1912.08005",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "MPP-2019-258, MCNET-19-27, LU-TP 19-58",
    doi = "10.1016/j.cpc.2020.107310",
    journal = "Comput. Phys. Commun.",
    volume = "260",
    pages = "107310",
    year = "2021"
}
        """

        self.AddCitations(citation)


    def AddElement(self,key,val,combine_dictionaries=False):
        if(key in self.metadata.keys()):
            old_val = self.metadata[key]
            if(isinstance(old_val,dict) and isinstance(val,dict) and combine_dictionaries):
                # special case: combine dictionaries
                for k,v in old_val.items():
                    val[k] = v
            else:
                self._print('Warning: Overwriting metadata associated with key={} .'.format(key))
        self.metadata[key] = val

    def AddCitations(self,val):
        """
        Special function for appending citation information to the metadata.
        It is expected that multiple objects will contribute to this key throughout
        the workflow, each supplying a dictionary; we combine these together.
        """
        key = 'Metadata.Citations'
        self.AddElement(key,val,combine_dictionaries=True)

    def GetMetaData(self):
        return self.metadata

    def Initialize(self):
        start_time = time.time()
        self.AddElement('Metadata.Timestamp',start_time)
        self.AddElement('Metadata.Timestamp.StringUTC',time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(start_time)))
        self.AddElement('Metadata.GitHash',self._get_git_revision_short_hash())
        self.AddElement('Metadata.HostName',self._get_hostname())
        self.AddElement('Metadata.UniqueID',str(uuid.uuid4()))
        self.AddElement('Metadata.UniqueIDShort',str(uuid.uuid4())[:5]) # a second, shorter random string -- probably more convenient to use, at the risk of a higher (but still tiny) collision rate
        self._init_citations()

    def _get_git_revision_short_hash(self): # see https://stackoverflow.com/a/21901260
        cwd = os.path.dirname(os.path.abspath(__file__))
        try:
            result = sub.check_output(['git', 'rev-parse', '--short', 'HEAD'],cwd=cwd).decode('ascii').strip()
        except:
            result = 'NO_GIT_HASH'
        return result

    def _get_hostname(self):
        try:
            result = sub.check_output(['hostname']).decode('ascii').strip()
        except:
            result='NO_HOSTNAME'
        return result

    def StashMetaDataInROOTFile(self,root_file:Union[List[str],str], cwd=None, tree_name:str='hepmc3_tree'):
        """
        This utility function adds the existing metadata (in self.metadata)
        to the "UserInfo" of a TTree in the input ROOT file.
        Note that this simply stashes the metadata, it does *not* produce
        branches marking which metadata entries correspond with which events.
        This is to be used, for example, for stashing metadata into HepMC3/ROOT
        files during the "generation" step.
        """
        if(isinstance(root_file,list)):
            for entry in root_file:
                self.StashMetaDataInROOTFile(entry,cwd,tree_name)
                return

        if(cwd is not None):
            root_file = '{}/{}'.format(cwd,root_file)

        if(not pathlib.Path(root_file).exists()):
            self._print('Input file {} does not exist.'.format(root_file))
            return

        # acting on a single file
        f = rt.TFile(root_file,'UPDATE')
        keys = [x.GetName() for x in f.GetListOfKeys()]
        if(tree_name not in keys):
            self._print('Input file {} does not contain tree {}.'.format(root_file,tree_name))
            self._print('Available keys in file:')
            for key in keys:
                self._print('\t{}'.format(key))
            f.Close()
            return

        t = f.Get(tree_name)
        self._add_to_ttree(t)
        t.Write("", rt.TObject.kOverwrite) # make sure the userinfo is saved
        f.Close()
        return

    def _add_to_ttree(self,tree:rt.TTree):
        """
        We store the information in a TMap, within the TTree's UserInfo
        (which is a TList). For dictionary-type information, we serialize
        using the json package.
        """
        user_info = tree.GetUserInfo()

        for key, value in self.metadata.items():
            if isinstance(value, int):
                # TParameter<int> for ints
                param = rt.TParameter(int)(key, value)
                user_info.Add(param)
            elif isinstance(value, float):
                # TParameter<double> for floats
                param = rt.TParameter(float)(key, value)
                user_info.Add(param)
            elif isinstance(value, str):
                # TNamed for strings (name=key, title=value)
                param = rt.TNamed(key, value)
                user_info.Add(param)
            elif isinstance(value,np.ndarray): # TODO: Add support for reading
                # TList of TParameter for numpy array
                param = rt.TList()
                param.SetName(key)

                for i,entry in enumerate(value):
                    if(value.dtype==int):
                        param.Add(rt.TParameter(int)('key[{}]'.format(i),entry))
                    else: # assume float
                        param.Add(rt.TParameter(float)('key[{}]'.format(i),entry))
                user_info.Add(param)

            elif isinstance(value, dict):
                # Convert dict to JSON string and store as TNamed
                json_str = json.dumps(value)
                param = rt.TNamed(key, json_str)
                param.SetUniqueID(999)  # Custom marker for JSON data, to tell it apart from the basic string
                user_info.Add(param)
            else:
                self._print('Warning: _add_to_ttree() unable to add metadata associated with key={} to ROOT file. It is of type {}.'.format(key,type(value)))
        return

    def _read_from_ttree(self,tree:rt.TTree):
        user_info = tree.GetUserInfo()
        metadata = {}

        for i in range(user_info.GetEntries()):
            obj = user_info.At(i)
            key = obj.GetName()

            if obj.InheritsFrom("TParameter<int>"):
                metadata[key] = obj.GetVal()
            elif obj.InheritsFrom("TParameter<float>"):
                metadata[key] = obj.GetVal()
            elif obj.InheritsFrom("TParameter<double>"):
                metadata[key] = obj.GetVal()
            elif obj.InheritsFrom("TNamed"):
                value = obj.GetTitle()
                # Check if this was originally a dictionary
                if obj.GetUniqueID() == 999:
                    metadata[key] = json.loads(value)
                else:
                    metadata[key] = value
        return metadata

    def ReadMetaDataFromROOTFile(self,root_file:str, cwd=None, tree_name:str='hepmc3_tree'):
        if(cwd is not None):
            root_file = '{}/{}'.format(cwd,root_file)

        if(not pathlib.Path(root_file).exists()):
            self._print('Input file {} does not exist.'.format(root_file))
            return None

        f = rt.TFile(root_file,'READ')
        keys = [x.GetName() for x in f.GetListOfKeys()]
        if(tree_name not in keys):
            self._print('Input file {} does not contain tree {}.'.format(root_file,tree_name))
            self._print('Available keys in file:')
            for key in keys:
                print('\t{}'.format(key))
            f.Close()
            return

        t = f.Get(tree_name)
        metadata = self._read_from_ttree(t)
        f.Close()
        return metadata

    def _print(self,val):
        print('{}: {}'.format(self.print_prefix,val))
        return

    def AddMetaDataWithReference(self,ntuple_file,cwd=None,overwrite=False, copts=9):
        """
        Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
        Also creates a column in the dataset that will point to this metadata's index.
        Somewhat redundant for file generation but this type of logic will be useful when concatenating files
        with different entries in the metadata fields.
        """
        if(self.configurator.GetReconstructionOutputFormat().lower() == 'hdf5'): # TODO: Rethink this? Reconstruction will always make ROOT -- conversion to HDF5 is a post-processing step, but nice to keep ability to handle metadata there?
            self.AddMetaDataWithReferenceH5(ntuple_file,cwd,overwrite,copts)
        elif(self.configurator.GetReconstructionOutputFormat().lower() == 'root'):
            self.AddMetaDataWithReferenceRoot(ntuple_file,self.configurator.GetReconstructionTreeName(),cwd)
        else:
            self._print('Warning: AddMetaDataWithReference() not implemented for file format {}.'.format(self.configurator.GetReconstructionOutputFormat()))
        return

    def AddMetaDataWithReferenceH5(self,ntuple_file,cwd=None,overwrite=False, copts=9):
        """
        Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
        Also creates a column in the dataset that will point to this metadata's index.
        Somewhat redundant for file generation but this type of logic will be useful when concatenating files
        with different entries in the metadata fields.
        """

        if(cwd is not None): ntuple_file = '{}/{}'.format(cwd,ntuple_file)
        f = h5.File(ntuple_file,'r+')
        check_key = list(f.keys())[0]
        nevents = f[check_key].shape[0]
        metadata = f.attrs

        for key,value in self.metadata.items():

            if(key.split('.')[0] != 'Metadata'):
                key = 'Metadata.{}'.format(key)

            # Dictionaries are not supported in HDF5, but we can convert to JSON.
            if(isinstance(value,dict)):
                value = json.dumps(value)

            if((key not in metadata.keys()) or overwrite):
                f.attrs[key] = [value]
                idx = 0
            else:
                # Check if this value already exists in the metadata list.
                # NOTE: This check is a bit unnecessary based on how this class is being used,
                #       but we might as well do it for flexibility in case we leverage this
                #       function in a different way. Typically, we're never going to encounter
                #       a duplicate value.
                if(value in f.attrs[key]):
                    idx = list(f.attrs[key]).index(value)
                else:
                    f.attrs[key] = list(f.attrs[key]) + [value] # I think the list <-> array stuff should be OK here
                    idx = len(f.attrs[key]) - 1

            if(key not in f.keys()):
                f.create_dataset(key,data=np.full(nevents,idx,dtype=np.dtype('i4')),compression='gzip',compression_opts=copts)
        f.close()

    def AddMetaDataWithReferenceRoot(self,ntuple_file,tree_name, cwd=None):
        """
        Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
        Also creates a branch in the dataset that will point to this metadata's index.
        Somewhat redundant for file generation but this type of logic will be useful when concatenating files
        with different entries in the metadata fields.
        """

        if(cwd is not None): ntuple_file = '{}/{}'.format(cwd,ntuple_file)
        f = rt.TFile(ntuple_file,'UPDATE')
        t = f.Get(tree_name)

        # We are modifying the TTree "in-place", so we actually make a temporary output file,
        # which we'll copy over the input file at the end.
        output_file = ntuple_file.replace('.root','_tmp_{}.root'.format(str(uuid.uuid4())))
        out_file = rt.TFile.Open(output_file, "RECREATE")
        out_file.SetCompressionAlgorithm(f.GetCompressionAlgorithm())
        out_file.SetCompressionLevel(f.GetCompressionLevel())

        # Clone the tree structure (no entries yet)
        out_tree = t.CloneTree(0)

        index_values = self._add_list_to_ttree(out_tree)

        # Create branches for the metadata references
        branches = {}
        buffers = {key:np.full(1,index_values[key],dtype=np.dtype('uint64')) for key in self.metadata.keys()}
        for key in self.metadata.keys():
            if(key.split('.')[0] != 'Metadata'):
                key = 'Metadata.{}'.format(key)
            branches[key] = out_tree.Branch(key, buffers[key],"{}/l".format(key))
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            out_tree.Fill()

        # Write and close
        out_file.cd()
        out_tree.Write()
        out_file.Close()
        f.Close()
        sub.check_call(['mv',output_file,ntuple_file])
        return

    def _add_list_to_ttree(self,tree:rt.TTree):
        """
        We store the information in a TMap, within the TTree's UserInfo
        (which is a TList). For dictionary-type information, we serialize
        using the json package.

        Here, we store TLists of metadata -- this is the way to store
        things in the final n-tuple, so that these lists can be appended
        to when concatenating datasets.
        """
        user_info = tree.GetUserInfo()
        existing_keys = [x.GetName() for x in user_info]
        index_values = {}

        for key, value in self.metadata.items():
            is_dictionary = False
            if(key in existing_keys):
                idx = user_info.Find(key).GetEntries()
            else:
                idx = 0
            name = '{}[{}]'.format(key,idx)
            index_values[key] = idx

            if isinstance(value, int):
                # TParameter<int> for ints
                param = rt.TParameter(int)(name, value)
            elif isinstance(value, float):
                # TParameter<double> for floats
                param = rt.TParameter(float)(name, value)
            elif isinstance(value, str):
                # TNamed for strings (name=key, title=value)
                param = rt.TNamed(name, value)
            elif isinstance(value,np.ndarray): # TODO: Add support for reading
                # TList of TParameter for numpy array
                param = rt.TList()
                param.SetName(name)

                for i,entry in enumerate(value):
                    if(value.dtype==int):
                        param.Add(rt.TParameter(int)('{}[{}]'.format(name,i),entry))
                    else: # assume float
                        param.Add(rt.TParameter(float)('{}[{}]'.format(name,i),entry))

            elif isinstance(value, dict):
                # Convert dict to JSON string and store as TNamed
                is_dictionary = True
                json_str = json.dumps(value)
                param = rt.TNamed(name, json_str)
                param.SetUniqueID(999)  # Custom marker for JSON data, to tell it apart from the basic string
            else:
                self._print('Warning: _add_list_to_ttree() unable to add metadata associated with key={} to ROOT file. It is of type {}.'.format(key,type(value)))
                return

            # Now, either package this metadata value into a TList, or add it to an existing one.
            if(key in existing_keys):
                param_list = user_info.Find(key)
                param_list.Add(param)

            else:
                param_list = rt.TList()
                param_list.SetName(key) # TODO: Is it an issue if the param inside has the same name?
                param_list.Add(param)
                if(is_dictionary):
                    param_list.SetUniqueID(999)
                user_info.Add(param_list)

        return index_values

    def _read_list_from_ttree(self,tree:rt.TTree):
        user_info = tree.GetUserInfo()
        metadata = {}

        for i in range(user_info.GetEntries()):
            obj_list = user_info.At(i)
            key = obj_list.GetName()

            if obj_list.At(0).InheritsFrom("TParameter<int>"):
                metadata[key] = [obj.GetVal() for obj in obj_list]
            elif obj_list.At(0).InheritsFrom("TParameter<float>"):
                metadata[key] = [obj.GetVal() for obj in obj_list]
            elif obj_list.At(0).InheritsFrom("TParameter<double>"):
                metadata[key] = [obj.GetVal() for obj in obj_list]
            elif obj_list.At(0).InheritsFrom("TNamed"):
                # Check if this was originally a dictionary
                if obj_list.GetUniqueID() == 999: # this is set on both obj_list and the objects inside it
                    metadata[key] = [json.loads(obj.GetTitle()) for obj in obj_list]
                else:
                    metadata[key] = [obj.GetTitle() for obj in obj_list]
        return metadata