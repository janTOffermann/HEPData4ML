
import sys,os,glob
import argparse as ap

sys.path.append( os.path.dirname(os.path.abspath(__file__)) + '/../..' )
from util.pileup.pileup import PileupOverlayPtFilter
from util.config import GetConfigDictionary, Configurator

def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFiles',type=str,help='Glob-compatible string for input HepMC3/ROOT pileup files.',required=True)
    parser.add_argument('-fastjet','--fastjetDir',type=str,help='FastJet installation directory. Defaults to `None`, which will use local install.',default=None)
    args = vars(parser.parse_args())

    input_files = args['inputFiles']
    fastjet_dir = args['fastjetDir']

    # Create a configurator -- used under-the-hood by pileup_handler for some FastJet configuration.
    config_dictionary = {'reconstruction':{'fastjet_dir':fastjet_dir}}
    configurator = Configurator(config_dictionary=config_dictionary)

    pileup_handler = PileupOverlayPtFilter(input_files,precompute=True)
    pileup_handler.SetConfigurator(configurator)
    pileup_handler.Initialize() # this will launch the computation of the leading jet pt

    return

if(__name__=='__main__'):
    main(sys.argv)