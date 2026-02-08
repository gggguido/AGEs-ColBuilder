# Collagen Crosslink Clustering Analysis

This folder contains scripts to analyze collagen microfibril extension, overlap, and gap from MD trajectories.

## What is analyzed

Each script computes, frame by frame:
- End-to-end extension from ACE and NME cap centroids
- Crosslink cluster positions from bonded marker residues
- Overlap and gap from a k=2 clustering split along the ACE->NME axis
- Uncertainty estimates from error propagation of centroid-position SEM into derived distances/ratios, combined with sliding-window bootstrap
- Replica-averaged time series and propagated replica-level errors

## Marker selection pattern

All scripts use marker arguments with a consistent interface:
- Duo crosslinks: `--m1-resname` and `--m2-resname`
- Trio crosslinks: `--m1-resname`, `--m2-resname`, and `--m3-resname`

## Scripts in this directory

- `glucosepane/glucosepane_analysis.py`: duo analysis (default markers `AGS`, `LGX`)
- `pentosidine/pentosidine_analysis.py`: duo analysis (default markers `APD`, `LPS`)
- `MOLD/mold_analysis.py`: duo analysis (default markers `LZS`, `LZD`)
- `pyd/pyd_analysis.py`: trio analysis (markers passed explicitly)
- `pyd+glucosepane/pyd+gluco_analysis.py`: trio analysis (markers passed explicitly)

## Common run format

Duo example:
```bash
python <script>.py \
  --top reduced.tpr \
  --traj run1/traj_vmd.xtc run2/traj_vmd.xtc run3/traj_vmd.xtc \
  --sel1 "resname ACE" --sel2 "resname NME" \
  --m1-resname <MARKER1> --m2-resname <MARKER2> \
  --start 0 --stop -1 --stride 1 \
  --bootstrap-window 10 --bootstrap-resamples 1000 \
  --out-prefix collagen67nm
```

Trio example:
```bash
python <script>.py \
  --top reduced.tpr \
  --traj run1/traj_vmd.xtc run2/traj_vmd.xtc run3/traj_vmd.xtc \
  --sel1 "resname ACE" --sel2 "resname NME" \
  --m1-resname <MARKER1> --m2-resname <MARKER2> --m3-resname <MARKER3> \
  --start 0 --stop -1 --stride 1 \
  --bootstrap-window 10 --bootstrap-resamples 1000 \
  --out-prefix collagen67nm
```

## Output files

For each run, scripts write four replica-mean CSV files using `--out-prefix`:
- `<prefix>_mean_all.csv`
- `<prefix>_mean_overlap.csv`
- `<prefix>_mean_gap.csv`
- `<prefix>_mean_ratio.csv`

Some scripts also write optional debug CSV files for detected crosslink components.
