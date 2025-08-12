import sys,os,pathlib,glob
import argparse as ap
import subprocess as sub

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-R', '--rundir', type=str, help='Run directory -- will prep and run the jobs defined by the submit file in here.', default='run0')
    args = vars(parser.parse_args())

    rundirs = glob.glob(args['rundir'])

    for rundir in rundirs:
        argument_file = rundir + '/arguments.txt'
        submit_file = 'condor.sub'

        with open(argument_file,'r') as f:
            lines = f.readlines()

        nlines = len(lines)
        for i in range(nlines):
            d = '{}/job{}'.format(rundir,i)
            if(pathlib.Path(d).exists()):
                sub.check_call(['rm','-r',d])
            os.makedirs(d,exist_ok=True)

        comm = ['condor_submit',submit_file]
        sub.check_call(comm,cwd=rundir)
    return

if(__name__=='__main__'):
    main(sys.argv)