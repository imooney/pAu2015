# EEC analysis

Minimal starting point for a STAR PicoDst analysis at RCF. As-is it reads
PicoDst events and prints one line per event -- nothing else. Build your real
analysis (track loops, cuts, histograms, output trees) inside
`StPicoAnalyzer::Make()`.

## Layout

```
StRoot/StPicoAnalyzer/   custom StMaker: StPicoAnalyzer.h / .cxx
macros/readPicoDst.C     steering macro: builds the chain, runs the event loop
submit/submit.xml        SUMS job description for batch submission
```

## Build (interactively, on an rcas node)

```bash
singularity shell --shell /usr/bin/csh -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif csh                      # Makes star libraries work with Alma9
setup 64b              # After the recent switch from 32 to 64 bit this can be necessary
starver SL23c          # pick the library version matching your PicoDsts
cons                    # compiles everything under StRoot/, from the dir above StRoot
```

`cons` needs to be run from the directory that *contains* `StRoot/` (i.e. the
top of this skeleton). It will create `.sl*` object/library files and drop
compiled libraries where ROOT can `gSystem->Load()` them -- that's what
`readPicoDst.C` does for `StPicoAnalyzer`.

## Run interactively

```bash
# single file
root4star -b -q 'macros/readPicoDst.C("/path/to/some.picoDst.root", -1)'

# a list of files (one path per line)
root4star -b -q 'macros/readPicoDst.C("picodst.list", -1)'

# just the first 100 events, useful while testing
root4star -b -q 'macros/readPicoDst.C("picodst.list", 100)'
```

`StPicoDstMaker` builds a `TChain` under the hood whether you give it one file
or a list, so the same macro handles both.

## Submit via SUMS

1. Build a file list of the PicoDsts you want to process (full paths, one per
   line), e.g.:
   ```bash
   find /star/data100/reco/.../PicoDst -name "*picoDst.root" > picodst.list
   ```
2. Edit `submit/submit.xml`:
   - set the `starver` tag to match your PicoDst production
   - point `USERNAME`/`PROJECT` paths at your own RCF space
   - point the `<input>` URL at your `picodst.list`
3. From the directory containing `submit.xml`, `readPicoDst.C`, and `StRoot/`:
   ```bash
   star-submit submit.xml
   ```
   Add `-s` first to do a dry run (writes per-job scripts without submitting),
   which is worth doing the first time to sanity check paths.

Per-job logs land wherever `<stdout>`/`<stderr>` point; per-job condor/sh
scripts land in `submit/sched/` (via the `<Generator>` block), which is handy
for debugging a job that failed.

## Where to go from here

- Add histograms/output trees in `StPicoAnalyzer::Init()` / `Make()` /
  `Finish()`, and write them out (e.g. `TFile`) so batch jobs produce
  something the `<output>` block in `submit.xml` can collect.
- Add track- and event-level cuts as member cut values on `StPicoAnalyzer`,
  set from the constructor or setters called in the steering macro.
- If the analysis grows enough parameters, consider adding them as
  constructor/macro arguments rather than hardcoding, so the same compiled
  code can run different cut variations without recompiling.
