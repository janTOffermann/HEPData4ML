import subprocess as sub
import sys,os
import argparse as ap

def main(args):
    parser = ap.ArgumentParser()

    parser.add_argument('-i','--input',type=str,default='param_card.dat',help='Input file.')
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

    # TODO: Quick fix for nproc: zero pad it to help with alphanumeric sorting.
    #       This is sometimes handy, e.g. when gluing together the outputs
    #       (to keep event_idx monotonically increasing). This padding is currently
    #       not adjustable, so we hope we do not have >10k jobs.
    nproc_string = str(nproc).zfill(4)

    os.makedirs(outdir,exist_ok=True) # make the output directory if it does not already exist
    input_nopath = input.split('/')[-1]
    if(not is_directory):
        if(input_extension is None):
            input_extension = input_nopath.split('.')[-1]
        output = input_nopath.replace('.' + input_extension,'') + '_{}.{}'.format(nproc_string,input_extension)

    else:
        output = input_nopath + '_{}'.format(nproc_string)
    outfile = '{}/{}'.format(outdir,output)

    comm = ['cp','-r',input,outfile] # add the -r flag in case we are copying a directory
    print('Copying:')
    print('\t Input: {}'.format(input))
    print('\tOutput: {}'.format(outfile))
    sub.check_call(comm)
    return

if(__name__=='__main__'):
    main(sys.argv)
