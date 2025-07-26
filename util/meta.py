import os, uuid, time
import subprocess as sub
import numpy as np
import h5py as h5

class MetaDataHandler:

    def __init__(self):

        self.metadata = {}
        self.Initialize()

    def AddElement(self,key,val):
        self.metadata[key] = val

    def GetMetaData(self):
        return self.metadata

    def Initialize(self):
        start_time = time.time()
        self.AddElement('Metadata.Timestamp',start_time)
        self.AddElement('Metadata.Timestamp.StringUTC',time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(start_time)))
        self.AddElement('Metadata.GitHash',self.get_git_revision_short_hash())
        self.AddElement('Metadata.UniqueID',str(uuid.uuid4()))
        self.AddElement('Metadata.UniqueIDShort',str(uuid.uuid4())[:5]) # a second, shorter random string -- probably more convenient to use, at the risk of a higher (but still tiny) collision rate

    def get_git_revision_short_hash(self): # see https://stackoverflow.com/a/21901260
        cwd = os.path.dirname(os.path.abspath(__file__))
        try:
            result = sub.check_output(['git', 'rev-parse', '--short', 'HEAD'],cwd=cwd).decode('ascii').strip()
        except:
            result = 'NO_GIT_HASH'
        return result

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

