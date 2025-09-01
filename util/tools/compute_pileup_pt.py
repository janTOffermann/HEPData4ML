
import sys,os,glob
import argparse as ap

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../..' )
from util.pileup.pileup import PileupOverlayPtFilter
from util.config import GetConfigDictionary, Configurator

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFiles',type=str,help='Glob-compatible string for input HepMC3/ROOT pileup files.',required=True)
    args = vars(parser.parse_args())

    input_files = args['inputFiles']

    # Create a configurator -- used under-the-hood by pileup_handler for some FastJet configuration.
    # config_dictionary = GetConfigDictionary(config_file)
    # configurator = Configurator(config_dictionary=config_dictionary)
    config_dictionary = {'reconstruction':{'fastjet_dir':None}}
    configurator = Configurator(config_dictionary=config_dictionary)

    pileup_handler = PileupOverlayPtFilter(input_files,precompute=True)
    pileup_handler.SetConfigurator(configurator)
    pileup_handler.Initialize()

    return

if(__name__=='__main__'):
    main(sys.argv)