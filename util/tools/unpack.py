import subprocess as sub
import sys, os, glob
import argparse as ap

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--input',type=str,required=True,help='Glob string for input tar.bz2 files.')
    args = vars(parser.parse_args())
    files = sorted(glob.glob(args['input']))


    for i,file in enumerate(files):
        print('{}/{}'.format(i+1,len(files)))
        command = 'tar -xjf {}'.format(file).split(' ')
        sub.check_call(command)


if(__name__ == '__main__'):
    main(sys.argv)