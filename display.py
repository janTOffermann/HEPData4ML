import sys,os
import ROOT as rt
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


    ev_disp = EventDisplay()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    card_file = "{}/util/delphes/cards/delphes_card_CMS.tcl".format(this_dir)

    ev_disp.InitializeDisplay(delphes_card = card_file)

    # app.Run(True)

    # print("Display started.")
    # input("Press Enter to exit...")
    # app.Terminate()


if(__name__=='__main__'):
    main(sys.argv)