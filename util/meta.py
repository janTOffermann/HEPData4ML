import os, uuid, time, pathlib, json
import subprocess as sub
import numpy as np
import h5py as h5
import ROOT as rt
from typing import Union, List

class MetaDataHandler:

    def __init__(self):

        self.metadata = {}
        self.print_prefix = 'MetaDataHandler: '
        self.Initialize()

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

    def AddMetaDataToROOTFiles(self,root_file:Union[List[str],str], cwd=None, tree_name:str='hepmc3_tree'):
        """
        This utility function adds the existing metadata (in self.metadata)
        to the "UserInfo" of a TTree in the input ROOT file.
        This is a convenient way to stash the metadata into a HepMC3/ROOT file,
        for example -- which may be handy in use cases such as creating input
        pileup files, where we want to store the metadata but don't have a final
        HDF5 output.
        """
        if(isinstance(root_file,list)):
            for entry in root_file:
                self.AddMetaDataToROOTFiles(entry,cwd,tree_name)
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
                print('\t{}'.format(key))
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
            elif isinstance(value, dict):
                # Convert dict to JSON string and store as TNamed
                json_str = json.dumps(value)
                param = rt.TNamed(key, json_str)
                param.SetUniqueID(999)  # Custom marker for JSON data, to tell it apart from the basic string
                user_info.Add(param)
            else:
                self._print('Warning: _add_to_ttree() unable to add metadata associated with key={} to ROOT file.'.format(key))
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

def AddMetaData(h5_file,cwd=None,value='',key='metadata'):
    """
    Generic function for adding a value to the HDF5 file
    metadata container.
    """
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    f.attrs[key] = value
    f.close()

def AddMetaDataWithReference(h5_file,cwd=None,value='',key='metadata',overwrite=False, copts=9):
    """
    Adds an entry to the metadata -- if under an existing key, appends it to the list at that key.
    Also creates a column in the dataset that will point to this metadata's index.
    Somewhat redundant for file generation but this type of logic will be useful when concatenating files
    with different entries in the metadata fields.
    """

    if(key.split('.')[0] != 'Metadata'):
        key = 'Metadata.{}'.format(key)

    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    check_key = list(f.keys())[0]
    nevents = f[check_key].shape[0]
    metadata = f.attrs

    # Dictionaries are not supported in HDF5, but we can convert to JSON.
    if(isinstance(value,dict)):
        value = json.dumps(value)

    if((key not in metadata.keys()) or overwrite): f.attrs[key] = [value]
    else: f.attrs[key] = list(f.attrs[key]) + [value] # I think the list <-> array stuff should be OK here
    idx = len(f.attrs[key]) - 1

    if(key not in f.keys()):
        f.create_dataset(key,data=np.full(nevents,idx,dtype=np.dtype('i4')),compression='gzip',compression_opts=copts)

    f.close()

def GetMetadata(h5_file,cwd=None,key='metadata'):
    result = None
    if(cwd is not None): h5_file = '{}/{}'.format(cwd,h5_file)
    f = h5.File(h5_file,'r+')
    if(key in f.attrs.keys()):
        result =  f.attrs[key]
    f.close()
    return result

