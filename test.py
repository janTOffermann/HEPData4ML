# Basic tests.
import sys,os,glob
import numpy as np

# path_prefix = os.getcwd() + '/../'
# if(path_prefix not in sys.path): sys.path.append(path_prefix)

from util.generation import Generator

print('Starting test.')

pt_min = 550
pt_max = 650

os.makedirs('test',exist_ok=True)

g = Generator(pt_min, pt_max)
g.SetOutputDirectory('test')

# The generator has automatically read in the config from config.py, in its __init__ function.
# Now let's try generating an event.
nevents = 20
g.Generate(nevents)

codes = g.GetUniqueProcessCodes()
print(codes)

xsecs = g.GetSigmaDictionary()
for key,val in xsecs.items():
    print('Process code: {}, xsec: {:.2e}, xsec error: {:.2e}'.format(key,*val))
