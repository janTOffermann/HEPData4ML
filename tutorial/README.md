# Tutorial

## Prerequisites

Before running the tutorials below, make sure that your environment is set up correctly (see the top-level README.md). The scripts will assume that this is set up.

## Tutorial 0

In this first tutorial, we will just run Pythia8 to produce some events corresponding with the processes in Pythia's `SoftQCD` category, without any special phase space cuts. To run this, simply invoke
```
./tutorial_0.sh
```

The corresponding config file is `config/config_0.py`. As you will see, this only contains a "generation" section, as we're not progressing farther than matrix element generation and showering/hadronization.

The outputs of the script will be in the `output/tutorial_0` directory; the output HepMC3 file is `events_0.root`.

## Tutorial 1

In this tutorial, we'll do something considerably more complex: we'll generate top-antitop pair production events, with a few special configurations:
- We'll run the full MC pipeline this time: generation, pileup, simulation and reconstruction steps.
    - We'll include pileup, taken from the output of Tutorial 0.
- We'll cluster large-radius (R=0.8) jets.
    - These jets will be ghost-associated with the truth-level top quark, and we'll also require that the top and bottom quarks, and the W boson (the b and W coming from the top quark decay) are within ∆R<0.8 of the jet centroid.
    - We'll also tag these jets with the Johns Hopkins top tagger, using its default setup (from Fastjet).
- We'll cluster small-radius (R=0.4) jets.
    - These will be ghost-associated with the large-radius jets -- and thus should basically serve as a way to capture some features of the large-radius jets' substructure.

To run this tutorial, invoke
```
./tutorial_1.sh
```

If this is your first time running the package, this may take a moment as a few dependencies may need to be installed; this will happen automatically and the printouts should give you a sense of the progress. For reference, on an M3 Max MacBook Pro this takes roughly a minute or so (the building of dependencies may use a handful of logical cores; it's not maxed out).

The output will be in `output/tutorial_1`, and should consist of a couple files -- the final n-tuple is `events.h5`. You can check the contents of this file by passing it to the script `../util/tools/check_file.py`, or visualize it via
```
python -i ../display.py -i output/tutorial_1/events.h5 -ei 0 -mode 0
# note the use of python -i otherwise it closes immediately; need to keep the viewer active
```
