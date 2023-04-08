import subprocess as sub
import sys,os
import argparse as ap

def main(args):
    parser = ap.ArgumentParser()

    parser.add_argument('-i','--input',type=str,default='param_card.dat',help='Input parameter card.')
    parser.add_argument('-o','--outdir',type=str,required=True,help='Output directory.')
    parser.add_argument('-n','--nProc',type=str,default='0',help='Process number (used for output naming).')
    parser.add_argument('-e','--extension',type=str,default=None)
    parser.add_argument('-d','--directory',type=int,default=0,help='Whether or not input is a directory.')
    args = vars(parser.parse_args())

    input = args['input']
    outdir = args['outdir']
    nproc = args['nProc']
    input_extension = args['extension']
    is_directory = args['directory'] > 0

    os.makedirs(outdir,exist_ok=True) # make the output directory if it does not already exist
    input_nopath = input.split('/')[-1]
    if(not is_directory):
        if(input_extension is None):
            input_extension = input_nopath.split('.')[-1]
        output = '.'.join(input_nopath.replace('.' + input_extension,'')) + '_{}.{}'.format(nproc,input_extension)

    else:
        output = input_nopath + '_{}'.format(nproc)
    outfile = '{}/{}'.format(outdir,output)

    comm = ['cp','-r',input,outfile] # add the -r flag in case we are copying a directory
    sub.check_call(comm)
    return

if(__name__=='__main__'):
    main(sys.argv)
