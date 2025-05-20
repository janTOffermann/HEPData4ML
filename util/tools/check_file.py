import sys,os,pathlib
import h5py as h5
import numpy as np
import argparse as ap

def _sort_key(s):
    components = s.split('.')
    return [(components[i], i, len(components[i]), components[i])
            for i in range(len(components))]

def _key_sort(keys):
    return sorted(keys, key=_sort_key)

def _check_file(filename):
    result = pathlib.Path(filename).exists()
    if(not result):
        print('Error: File {} not found.'.format(filename))
    return result

def _print_keys(filename):
    f = h5.File(filename,'r')
    # keys = sorted(list(f.keys()))
    keys = _key_sort(list(f.keys()))
    shapes = {k:f[k].shape for k in keys}
    max_key_length = np.max([len(k) for k in keys])

    # determine which keys correspond with metadata (attributes)
    metakeys = f.attrs.keys()
    metashapes = {k:f[k].shape for k in metakeys}
    for key in keys:
        if(key in metashapes.keys()):
            del shapes[key]
    f.close()

    print(2 * max_key_length * '=')
    print('Metadata:')
    for key,val in metashapes.items():
        printkey = '[' + key + ']'
        while(len(printkey) < max_key_length + 2):
            # printkey = ' ' + printkey
            printkey += '-'
        print('{}: {}'.format(printkey,val))

    print(2 * max_key_length * '-')
    print('Data:')
    for key,val in shapes.items():
        printkey = '[' + key + ']'
        while(len(printkey) < max_key_length + 2):
            # printkey = ' ' + printkey
            printkey += '-'
        print('{}: {}'.format(printkey,val))
    print(2 * max_key_length * '=')

    return

def _print2d(array):
    for i in range(array.shape[0]):
        print('[{}]: {}'.format(i,array[i]))
    return

def _print3d(array,lims=[]):
    for i in range(array.shape[0]):
        if(len(lims) > 0 and lims[0] is not None):
            if(i >= lims[0]): break
        print('[{}]:'.format(i))
        for j in range(array.shape[1]):
            if(len(lims) > 1 and lims[1] is not None):
                if(j >= lims[1]): break
            print('\t[{}]: {}'.format(j,array[i,j]))
    return

def _check_event(filename,idx=0):
    f = h5.File(filename,'r')
    metakeys = f.attrs.keys()
    keys = [x for x in _key_sort(list(f.keys())) if x not in metakeys]
    # keys = [x for x in sorted(list(f.keys())) if x not in metakeys]
    max_key_length = np.max([len(k) for k in keys])

    for key in keys:
        print(2 * max_key_length * '-')
        printkey = '[' + key + ']'
        while(len(printkey) < max_key_length + 2):
            printkey = ' ' + printkey
        print('{}:'.format(printkey))
        data = f[key][idx]

        if('Pmu' in key):
            if('Constituents' in key):
                try:
                    lkey1 = key.replace('Constituents.Pmu_cyl','N')
                    lkey1 = lkey1.replace('Constituents.Pmu','N')

                    lkey2 = key.replace('Pmu_cyl','N')
                    lkey2 = lkey2.replace('Pmu','N')

                    l1 = f[lkey1][idx] + 1 # add one, to check the zero-padding
                    l2 = np.max(f[lkey2][idx] + 1) # add one, to check the zero-padding

                    data = f[key][idx,:l1,:l2]
                except:
                    pass

            else:
                try:
                    lkey = key.replace('Pmu_cyl','N')
                    lkey = lkey.replace('Pmu','N')
                    l = f[lkey][idx] + 1 # add one, to check the zero-padding
                    data = f[key][idx,:l]
                except:
                    pass

        if(data.ndim == 1):
            print(data)
        elif(data.ndim == 2):
            _print2d(data)
        elif(data.ndim == 3):
            _print3d(data)
        else:
            print(data)

    print(2 * max_key_length * '-')

    f.close()
    return

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFile',type=str,required=True)
    args = vars(parser.parse_args())
    input_file = args['inputFile']

    # First check that the file exists
    if(not _check_file(input_file)):
        return

    # Now open the file, and gather the keys.
    _print_keys(input_file)

    # Check a single event
    print('Checking first event.')
    _check_event(input_file)
    return

if(__name__=='__main__'):
    main(sys.argv)