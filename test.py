# Basic tests.
import sys,os,glob
import numpy as np

# path_prefix = os.getcwd() + '/../'
# if(path_prefix not in sys.path): sys.path.append(path_prefix)

from util.generation import Generator

print('Starting test.')

pt_min = 550
pt_max = 650

g = Generator(pt_min, pt_max)