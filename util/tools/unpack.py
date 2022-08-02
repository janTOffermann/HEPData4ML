import subprocess as sub
import sys, os, glob

files = glob.glob("*.bz2")

for file in files:

    command = 'tar -xjf {}'.format(file).split(' ')
    sub.check_call(command)