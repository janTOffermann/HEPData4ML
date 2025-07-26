import sys,os
import ROOT as rt
import argparse as ap
from util.display.viewer import EventDisplay

def main(args):
    # NOTE: Need to run with python -i, otherwise this will crash
    # # See: https://root-forum.cern.ch/t/th3-tf3-and-gl-crashes/37863/11
    # rt.PyConfig.StartGuiThread = 'inputhook' # <- doesn't seem to work?

    # app = rt.TApplication("MyApp", 0, [])

    # # Create a simple control panel with quit button
    # main_frame = rt.TGMainFrame(rt.gClient.GetRoot(), 200, 100)
    # quit_button = rt.TGTextButton(main_frame, "Quit")
    # quit_button.Connect("Clicked()", "TPython", rt.TPython.Exec, "quit_application()")

    # main_frame.AddFrame(quit_button)
    # main_frame.SetWindowName("Control Panel")
    # main_frame.MapSubwindows()
    # main_frame.Resize(main_frame.GetDefaultSize())
    # main_frame.MapWindow()

    this_dir = os.path.dirname(os.path.abspath(__file__))

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFile',type=str,required=True)
    parser.add_argument('-c','--cardFile',type=str,default=None)
    parser.add_argument('-ei','--eventIndex',type=int,default=0)
    parser.add_argument('-mode','--mode',type=int,default=0)

    parser.add_argument('-jetPtMin','--jetPtMin',type=float,default=15.)
    parser.add_argument('-trackPtMin','--trackPtMin',type=float,default=1.)
    parser.add_argument('-truthParticlePtMin','--truthParticlePtMin',type=float,default=1.)

    args = vars(parser.parse_args())

    input_file = args['inputFile']
    card_file = args['cardFile']
    event_index = args['eventIndex']
    jet_pt_min = args['jetPtMin']
    track_pt_min = args['trackPtMin']
    truth_particle_pt_min = args['truthParticlePtMin']

    mode = args['mode']

    if(card_file is None):
        card_file = "{}/util/delphes/cards/delphes_card_CMS.tcl".format(this_dir)

    ev_disp = EventDisplay()
    ev_disp.InitializeDisplay(delphes_card = card_file, mode=mode)
    ev_disp.SetJetPtMin(jet_pt_min)
    ev_disp.SetTrackPtMin(track_pt_min)
    ev_disp.SetTruthParticlePtMin(truth_particle_pt_min)

    ev_disp.Display(input_file,event_index)


if(__name__=='__main__'):
    main(sys.argv)