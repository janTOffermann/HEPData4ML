# Tutorial

## Prerequisites

Before running the tutorials below, make sure that your environment is set up correctly (see the top-level README.md). The scripts will assume that this is set up.

## Tutorial 0

**Focus:** Basic event generation.

In this first tutorial, we will just run Pythia8 to produce some events corresponding with the processes in Pythia's `SoftQCD` category, without any special phase space cuts. To run this, simply invoke
```
./tutorial_0.sh
```

The corresponding config file is `config/config_0.py`. As you will see, this only contains a "generation" section, as we're not progressing farther than matrix element generation and showering/hadronization.

The outputs of the script will be in the `output/tutorial_0` directory; the output HepMC3 file is `events_0.root`.

## Tutorial 1

**Focus:** Basic event generation, with ASCII output.

This is a duplicate of **Tutorial 0**, except that it produces an ASCII-format HepMC3 file. These are a little less useful in practice because they lack some of the nice features of ROOT (specifically, random access -- which is super-useful for things like pileup handling), but they do have the benefit of being human-readable!

## Tutorial 2

**Focus:** Basic full pipeline -- event generation, fast detector simulation, object reconstruction.

In this tutorial, we'll generate top-antitop pair production events. Here's a breakdown of our configuration:
- We'll run nearly the full MC pipeline this time: generation, simulation and reconstruction steps. We're just skipping the addition of any pileup.
- We'll cluster large-radius (R=0.8) jets.
    - We'll enforce some jet-level pT and eta cuts.

To run this tutorial, invoke
```
./tutorial_1.sh
```

If this is your first time running the package, this may take a moment as a few dependencies may need to be installed; this will happen automatically and the printouts should give you a sense of the progress. For reference, on an M3 Max MacBook Pro this takes roughly a minute or so (the building of dependencies may use a handful of logical cores; it's not maxed out).

The output will be in `output/tutorial_1`, and should consist of a couple files -- the final n-tuple is `events.h5`. You can check the contents of this file by passing it to the script `../util/tools/check_file.py`.

## Tutorial 3

**Focus:** More complex full pipeline -- including basic pileup from Tutorial 0, and multiple jet definitions for reconstruction.

In this tutorial, we'll do something considerably more complex: we'll generate top-antitop pair production events again, now with a few special configurations:
- We'll run the full MC pipeline this time: generation, pileup, simulation and reconstruction steps.
    - We'll include pileup, taken from the output of Tutorial 0.
- We'll cluster large-radius (R=0.8) jets.
    - These jets will be ghost-associated with the truth-level top quark, and we'll also require that the top and bottom quarks, and the W boson (the b and W coming from the top quark decay) are within ∆R<0.8 of the jet centroid.
    - We'll also tag these jets with the Johns Hopkins top tagger, using its default setup (from Fastjet).
- We'll cluster small-radius (R=0.4) jets.
    - These will be ghost-associated with the large-radius jets -- and thus should basically serve as a way to capture some features of the large-radius jets' substructure.

To run this tutorial, invoke
```
./tutorial_2.sh
```

If this is your first time running the package, this may take a moment as a few dependencies may need to be installed; this will happen automatically and the printouts should give you a sense of the progress. For reference, on an M3 Max MacBook Pro this takes roughly a minute or so (the building of dependencies may use a handful of logical cores; it's not maxed out).

The output will be in `output/tutorial_2`, and should consist of a couple files -- the final n-tuple is `events.h5`. You can check the contents of this file by passing it to the script `../util/tools/check_file.py`, or visualize it via
```
python -i ../display.py -i output/tutorial_2/events.h5 -ei 0 -mode 0
# note the use of python -i otherwise it closes immediately; need to keep the viewer active
```

You can also use the `check_file.py` utility to print out the citations (in BibTex format) associated with the various algorithms we've used -- this includes packages like Pythia8, FastJet and Delphes, as well as algorithms like the Johns Hopkins top tagger and ghost association. This can be accomplished via
```
python ../util/tools/check_file.py -i events.h5 --citations
```

## Tutorial 4

**Focus:** Exploring pileup handling options -- on-the-fly generation.

This tutorial is similar to Tutorials 2 and 3:
- We perform the same object reconstructions as in Tutorial 2,
- We use pileup as in Tutorial 3.

However, the pileup is in this case generated on-the-fly; rather than reading in pileup events from some input HepMC3 file, we generate them with Pythia8 as needed (harnessing the same instance of the generator we used in the `generation` step). While we might sometimes be OK with recycling pileup events -- as may happen when randomly sampling them from a file -- this is a sure-fire way to generate unique pileup events. This may be especially handy when running parallelized jobs: We have some pileup handling classes like `PileupOverlayPtFilter` that are meant to avoid recycling pileup events that have a jet above some $p_{T}$ threshold (following the [typical ATLAS prescription](https://iopscience.iop.org/article/10.1088/1742-6596/513/2/022024)), these methods currently break down under parallelization where each worker/job isn't aware of whether another has already used a particular event.

## Tutorial 5

**Focus:** Exploring pileup handling options -- configuring more advanced options.

This tutorial is similar to Tutorial 5. However, we will use the `PileupOverlayPtFilter`. This pileup handler works with input HepMC3/ROOT files, and it is designed so that events with a leading jet $p_{T}$ below a user-defined cut can be reused, whereas those above cannot. Furthermore, those events above this cut will be divided up among jobs when running parallel jobs with HTCondor via `prep_condor.py` -- we don't explore that functionality in this tutorial, but it's useful to have an example of this nonetheless since it's a configuration you may wish to use in practice.

The `PileupOverlayPtFilter` is configured with `precompute=True`, which will cause it to compute the leading jet $p_{T}$ for each input pileup event (using FastJet). These will be saved in a new `TTree` within the input files; if this tree is present (from previous running), the `PileupOverlayPtFilter` will simply read this and skip the computation. Two caveats/notes:
- As a result of this approach, this method only works with HepMC3/ROOT format, and not HepMC3/ASCII. The `PileupOverlayPtFilter` will reject any input files that are not identified as being ROOT format.
- When running HTCondor jobs prepared by `prep_condor.py`, in order to prevent I/O clashes only the 1st job will attempt to write the `TTree` if it isn't found; all the others will try to read this tree if it exists, and otherwise will just compute the leading jet $p_{T}$ in memory but not save it to disk. This can potentially add quite a significant extra computational load to the jobs -- which will each be computing the same jets --  so it is strongly advised to make sure these `TTree`s are already written (either by launching a single test job first, or using a dedicated tool).