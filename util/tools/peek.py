import sys,os,time,datetime,re,pathlib
sys.path.append(str(pathlib.Path(os.path.dirname(os.path.abspath(__file__)) + '/../../').resolve()))
import argparse as ap
import subprocess as sub
from util.generation import Generator
from util.conversion import Processor
from util.config import Configurator, GetConfigDictionary

def none_or_str(value): # see https://stackoverflow.com/a/48295546
    if value == 'None':
        return None
    return value

def get_git_revision_short_hash(): # see https://stackoverflow.com/a/21901260
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        result = sub.check_output(['git', 'rev-parse', '--short', 'HEAD'],cwd=cwd).decode('ascii').strip()
    except:
        result = 'NO_GIT_HASH'
    return result

# Convenience function for file naming
def float_to_str(value):
    value_str = str(value)
    value_str = value_str.replace('.',',')
    return re.sub(',0$','',value_str)

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-pc',           '--pythia_config',     type=none_or_str,  default=None,             help='Path to Pythia configuration template (for setting the process).')
    parser.add_argument('-config',       '--config',            type=str,          default=None,             help='Path to configuration Python file. Default will use config/config.py .')

    args = vars(parser.parse_args())

    start_time = time.time()

    # nevents_per_bin = 1
    pt_bin_edges = [550, 560]
    verbose = True
    pythia_rng = 1
    pythia_config = args['pythia_config']
    # debug = True
    config_file = args['config']

    nbins = len(pt_bin_edges) - 1


    # Configurator class, used for fetching information from our config file.
    # We import this from a user-supplied file, by default it is config/config.py.
    this_dir = os.path.dirname(os.path.abspath(__file__))
    if(config_file is None):
        config_file = str(pathlib.Path('{}/../../config/config.py'.format(this_dir)).resolve())

    # config_dictionary = config.config
    print('Using configuration file: {} .'.format(config_file))
    config_dictionary = GetConfigDictionary(config_file)
    configurator = Configurator(config_dictionary=config_dictionary)

    print()
    dummy_processor = Processor(configurator)
    line_up = '\033[1A'
    line_clear = '\x1b[2K'
    for i in range(13):
        print(line_up, end=line_clear)

    print('\n=================================')
    print('Running Pythia8 event generation.')
    print('=================================\n')

    if(pythia_rng is not None):
        print('\tSetting Pythia RNG seed to {}. (overriding config)'.format(pythia_rng))
    else:
        pythia_rng = configurator.GetPythiaRNGSeed()
    if(pythia_config is not None):
        print('\tSetting Pythia process configuration from {}. (overriding config)'.format(pythia_config))

    print()
    for i in range(nbins):
        # Generate a HepMC file containing our events. The generation already performs filtering
        # before writing, so that we just save final-state particles & some selected truth particles.

        # If the user has opted not to do generation, the HepMC3 files must already exist (and have the right names).
        # TODO: Make the no-generation option more flexible, to pick up any existing HepMC3 files in the cwd.
        pt_min = pt_bin_edges[i]
        pt_max = pt_bin_edges[i+1]
        # hep_file = 'events_{}-{}.hepmc'.format(float_to_str(pt_min),float_to_str(pt_max))

        generator = Generator(pt_min,pt_max, configurator, pythia_rng,pythia_config_file=pythia_config,verbose=verbose)
        generator.SetEventSelection(configurator.GetEventSelection())
        generator.SetTruthSelection(configurator.GetParticleSelection())
        generator.SetFinalStateSelection(configurator.GetFinalStateSelection())
        generator.SetEventFilter(configurator.GetEventFilter())
        generator.SetEventFilterFlag(configurator.GetEventFilterFlag())
        generator.SetJetConfig(configurator.GetJetConfig())
        generator.SetNTruth(configurator.GetNPars()['n_truth'])

        # get some actual printouts
        generator.pythia_config['Print:quiet'] = 'off'
        generator.pythia_config['Stat:showProcessLevel'] = 'on'
        generator.pythia_config['Next:numberShowProcess'] = '1'
        generator.pythia_config['Next:numberShowEvent'] = '1'
        # for generator to read pythia_config again, must reinitialize
        generator.pythia.AddToConfigDict(generator.pythia_config)
        generator.pythia.InitializePythia()

        generator.GenerateSingle()

    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_time_readable = str(datetime.timedelta(seconds=elapsed_time))
    print('\n#############################')
    print('Done. Time elapsed = {:.1f} seconds.'.format(elapsed_time))
    print('({})'.format(elapsed_time_readable))
    print('#############################\n')

if __name__ == '__main__':
    main(sys.argv)