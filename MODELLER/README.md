# MODELLER Crosslink IC Generation

This folder contains Python utilities (`create_ic_*.py`) that build CHARMM/MODELLER-style crosslink templates from a SMILES string.

## What these scripts do

Each `create_ic_*.py` script:
- reads a SMILES string for a crosslink molecule,
- generates a 3D geometry with RDKit,
- assigns atom names and atom types,
- computes bonds, angles, and dihedrals,
- writes a `.str` file that includes an **IC (Internal Coordinates) table**.

Current examples in this folder include scripts for glucosepane (AGS), pentosidine (APD), and MOLD (LZS).

## Main output

Running a script writes a `*_ic.str` file (for example `glucosepane_ic.str`, `pentosidine_ic.str`, `mold_ic.str`).
The lines starting with `IC` are the IC list that is required to define the new crosslink in ColBuilder.

## How this is used in ColBuilder

The IC list created and reported in the generated `.str` files must then be transferred into:
- `colbuilder/src/colbuilder/data/sequence/modeller/top_heav_mod.lib`

This is the library file where new crosslinks are added for ColBuilder.

For reference, an example `top_heav_mod.lib` is also reported in this repository at:
- `MODELLER/top_heav_mod.lib`

## SMILES from a drawn structure

If you do not already have a SMILES string, you can draw the molecule in a free web editor and export/copy SMILES.
A commonly used free option is:
- [Cheminfo SMILES generator/checker](https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html)

## Typical usage

```bash
python glucosepane/create_ic_glucosepane.py
python pentosidine/create_ic_pentosidine.py
python MOLD/create_ic_mold.py
```

To generate a different crosslink, edit the `smiles = "..."` value at the end of the corresponding script and run it again.
