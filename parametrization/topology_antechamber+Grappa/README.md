# Topology Files from Antechamber and Grappa

This folder contains topology files for each available crosslink chemistry using two parametrization schemes.

## What is included

For each chemistry, two topology files are provided:

- `topol.top`: topology obtained with the Antechamber-based parametrization.
- `grappa_topol.top`: topology obtained from the corresponding Antechamber topology using Grappa.

The currently available chemistries are:

- `glucosepane`
- `pentosidine`
- `MOLD`

In the same folder, `amber99sb-star-ildnp.ff.zip` contains the updated `amber99sb-star-ildnp.ff` force-field directory with the extra parameters required to run MD with the Antechamber-derived topologies.

## Grappa reference

Grappa is available at:

- [Grappa GitHub repository](https://github.com/graeter-group/grappa)

## Command used to generate the Grappa topologies

The Grappa topologies were generated from the Antechamber `topol.top` files with:

```bash
grappa_gmx -f topol.top -o grappa_topol.top -t grappa-1.4 -p
```
