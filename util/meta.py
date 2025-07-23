import os, uuid, time
import subprocess as sub

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
        self.AddElement('timestamp',start_time)
        self.AddElement('timestamp_string_utc',time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime(start_time)))
        self.AddElement('git_hash',self.get_git_revision_short_hash())
        self.AddElement('unique_id',str(uuid.uuid4()))
        self.AddElement('unique_id_short',str(uuid.uuid4())[:5]) # a second, shorter random string -- probably more convenient to use, at the risk of a higher (but still tiny) collision rate

    def get_git_revision_short_hash(self): # see https://stackoverflow.com/a/21901260
        cwd = os.path.dirname(os.path.abspath(__file__))
        try:
            result = sub.check_output(['git', 'rev-parse', '--short', 'HEAD'],cwd=cwd).decode('ascii').strip()
        except:
            result = 'NO_GIT_HASH'
        return result