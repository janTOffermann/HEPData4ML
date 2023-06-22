import uuid
import numpy as np
import h5py as h5

def RN():
    return str(uuid.uuid4())

class MetaDataHandler:
    def __init__(self,filenames):
        files = [h5.File(x,'r') for x in filenames]
        attrs_list = [dict(x.attrs).copy() for x in files]

        # Now merge the dictionary is attrs_list
        self.attrs = {}
        keys = [list(x.keys()) for x in attrs_list]
        keys = list(set([j for i in keys for j in i]))

        for key in keys:
            self.attrs[key] = []
            for attrs in attrs_list:
                self.attrs[key].append(attrs[key])

            self.attrs[key] = np.unique(np.concatenate(self.attrs[key],axis=0))
            # print(key,self.attrs[key])

        for f in files:
            f.close()

    def GetMetaData(self):
        return self.attrs


