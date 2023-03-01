import sys,os,glob,uuid
import numpy as np
import h5py as h5
import argparse as ap
import subprocess as sub

def IsBlankLine(line):
    if(line.strip('\n').strip('\t').strip(' ') == ''): return True
    return False

def concatenate(input_patterns,output):
    # Determine what are the input files.
    input_files = []
    for pattern in input_patterns:
        if(',' in pattern):
            input_files += [x.strip() for x in pattern.split(',')]
        elif('*' in pattern or '[' in pattern):
            input_files += glob.glob(os.path.expanduser(pattern),recursive=True)
        else:
            input_files += [pattern]

    input_files.sort()
    print('Concatenating files:')
    for i,file in enumerate(input_files):
        print('\t{}: {}'.format(i,file))

    # Get the HepMC3 header from the first file.
    header_file = str(uuid.uuid4()) + '.txt'
    comm = ['head', '-n', '2',input_files[0],'>',header_file]
    sub.check_call(' '.join(comm),shell=True)

    # Get the HepMC footer from the first file.
    footer_file = str(uuid.uuid4()) + '.txt'
    comm = ['tail', '-n', '3',input_files[0], '|','head','-n','1','>',footer_file]
    sub.check_call(' '.join(comm),shell=True)

    comm = ['cat']
    for file in input_files:
        comm.append(file)
    comm.append('>')
    output_tmp = str(uuid.uuid4()) + '.hepmc'
    comm.append(output_tmp)

    header_lines = []
    footer_lines = []
    with open(header_file,'r') as f:
        header_lines += f.readlines()
    with open(footer_file,'r') as f:
        footer_lines+= f.readlines()

    sub.check_call(['rm',header_file])
    sub.check_call(['rm',footer_file])
    sub.check_call(' '.join(comm),shell=True)

    with open(output_tmp,'r') as f:
        lines = f.readlines()
    lines = [x for x in lines if x not in header_lines + footer_lines]
    lines = [x for x in lines if not IsBlankLine(x)]
    lines = header_lines + lines + footer_lines

    # Now, we have to fix the event numbers.
    event_number = 1 # uses 1-indexing
    for i,line in enumerate(lines):
        if(line[0] != 'E'): continue
        data = line.split(' ')
        data[1] = str(event_number)
        event_number += 1
        lines[i] = ' '.join(data)

    with open(output,'w') as f:
        for line in lines:
            f.write(line)

    sub.check_call(['rm',output_tmp])

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input',   nargs='+', help='Input file pattern (or a list of patterns). Each pattern can also be provided as a comma-separated string of filenames.', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file.', default=None)
    args = vars(parser.parse_args())

    input_patterns = args['input']
    output = args['output']

    concatenate(input_patterns,output)

if __name__ == '__main__':
    main(sys.argv)
