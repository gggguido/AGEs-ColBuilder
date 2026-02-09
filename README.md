# AGEs-ColBuilder Extension Repository

This repository collects scripts, inputs, and example systems to help ColBuilder users include and validate new collagen crosslink chemistries that are not yet implemented by default.

## Author

- **Guido Giannetti**

## Publication

This work is currently in preparation.

- **Citation:** _in preparation_
- **DOI placeholder:** `DOI: TBA`

## Repository Structure

### [`AGEs_candidate_positions`](./AGEs_candidate_positions/)
Tools to generate periodic crystal-contact copies of a collagen triple helix and identify candidate crosslink sites with geometric distance filters.

- [`generate_crystalcontact_copies.py`](./AGEs_candidate_positions/generate_crystalcontact_copies.py): expands the crystallographic environment (Chimera crystal contacts).
- [`calculate_dist_CA_CA_lys_arg.py`](./AGEs_candidate_positions/calculate_dist_CA_CA_lys_arg.py): finds **LYS-ARG** candidates (glucosepane and pentosidine).
- [`calculate_dist_CA_CA_lys_lys.py`](./AGEs_candidate_positions/calculate_dist_CA_CA_lys_lys.py): finds **LYS-LYS** candidates (MOLD).
- [`README.md`](./AGEs_candidate_positions/README.md): usage details.

### [`MODELLER`](./MODELLER/)
Scripts to generate internal-coordinate (IC) tables for new crosslink residues from SMILES geometry.

- [`glucosepane/create_ic_glucosepane.py`](./MODELLER/glucosepane/create_ic_glucosepane.py)
- [`pentosidine/create_ic_pentosidine.py`](./MODELLER/pentosidine/create_ic_pentosidine.py)
- [`MOLD/create_ic_mold.py`](./MODELLER/MOLD/create_ic_mold.py)
- Example generated IC stream files:
  - [`glucosepane/glucosepane_ic.str`](./MODELLER/glucosepane/glucosepane_ic.str)
  - [`pentosidine/pentosidine_ic.str`](./MODELLER/pentosidine/pentosidine_ic.str)
  - [`MOLD/mold_ic.str`](./MODELLER/MOLD/mold_ic.str)
- ColBuilder modeller library reference file:
  - [`top_heav_mod.lib`](./MODELLER/top_heav_mod.lib)
- Folder documentation:
  - [`MODELLER/README.md`](./MODELLER/README.md)

The workflow is general and can be adapted to additional crosslink chemistries beyond the examples included here.

### [`parametrization`](./parametrization/)
Detailed example workflow for parametrizing a new crosslink molecule and obtaining an **Amber99-compatible** parameter set.

- [`parametrization/README.md`](./parametrization/README.md): complete step-by-step commands.
- [`parametrization/glucosepane/`](./parametrization/glucosepane/): generated files for the glucosepane example.

### [`all_atom_md`](./all_atom_md/)
Detailed all-atom MD protocol and inputs used after generating collagen fibrils with ColBuilder.

- [`all_atom_md/md_steps.md`](./all_atom_md/md_steps.md): complete list of GROMACS commands used in this work.
- Minimization inputs:
  - [`all_atom_md/minimization/ions.mdp`](./all_atom_md/minimization/ions.mdp)
  - [`all_atom_md/minimization/minim_vacuum.mdp`](./all_atom_md/minimization/minim_vacuum.mdp)
  - [`all_atom_md/minimization/minim_solv.mdp`](./all_atom_md/minimization/minim_solv.mdp)
- Equilibration helpers:
  - [`all_atom_md/equilibration/A.NVT/2500/update_posre_2500.py`](./all_atom_md/equilibration/A.NVT/2500/update_posre_2500.py)
  - [`all_atom_md/equilibration/B.NPT/300/update_posre_300.py`](./all_atom_md/equilibration/B.NPT/300/update_posre_300.py)
- Production inputs:
  - [`all_atom_md/production/`](./all_atom_md/production/)

### [`equilibrated_systems`](./equilibrated_systems/)
Archived equilibrated systems used in this work.

- [`equilibrated_systems/README.md`](./equilibrated_systems/README.md)
- Zipped equilibrated `.gro` systems:
  - [`equilibrated_systems/glucosepane/glucosepane_eq.gro.zip`](./equilibrated_systems/glucosepane/glucosepane_eq.gro.zip)
  - [`equilibrated_systems/pentosidine/pentosidine_eq.gro.zip`](./equilibrated_systems/pentosidine/pentosidine_eq.gro.zip)
  - [`equilibrated_systems/MOLD/MOLD_eq.gro.zip`](./equilibrated_systems/MOLD/MOLD_eq.gro.zip)
  - [`equilibrated_systems/pyd/pyd_eq.gro.zip`](./equilibrated_systems/pyd/pyd_eq.gro.zip)
  - [`equilibrated_systems/pyd+glucosepane/pyd+glucosepane_eq.gro.zip`](./equilibrated_systems/pyd+glucosepane/pyd+glucosepane_eq.gro.zip)

### [`analysis`](./analysis/)
Post-processing scripts used to analyze fibrils under load.

- [`analysis/README.md`](./analysis/README.md): analysis overview and run format.
- Scripts:
  - [`analysis/glucosepane/glucosepane_analysis.py`](./analysis/glucosepane/glucosepane_analysis.py)
  - [`analysis/pentosidine/pentosidine_analysis.py`](./analysis/pentosidine/pentosidine_analysis.py)
  - [`analysis/MOLD/mold_analysis.py`](./analysis/MOLD/mold_analysis.py)
  - [`analysis/pyd/pyd_analysis.py`](./analysis/pyd/pyd_analysis.py)
  - [`analysis/pyd+glucosepane/pyd+gluco_analysis.py`](./analysis/pyd+glucosepane/pyd+gluco_analysis.py)

The scripts compute:
- D-band length proxy (end-to-end extension)
- overlap length
- gap length
- overlap/gap ratio
- associated uncertainties

## Suggested Usage Order

1. Identify candidate residue pairs in [`AGEs_candidate_positions`](./AGEs_candidate_positions/).
2. Build crosslink IC definitions in [`MODELLER`](./MODELLER/) and update the modeller library.
3. Parametrize the new chemistry in [`parametrization`](./parametrization/) (if needed).
4. Build fibrils with ColBuilder and run MD following [`all_atom_md`](./all_atom_md/).
5. Use [`equilibrated_systems`](./equilibrated_systems/) examples as references.
6. Analyze outputs with scripts in [`analysis`](./analysis/).

## Notes

- Each subfolder contains detailed documentation:
  - [`AGEs_candidate_positions/README.md`](./AGEs_candidate_positions/README.md)
  - [`MODELLER/README.md`](./MODELLER/README.md)
  - [`parametrization/README.md`](./parametrization/README.md)
  - [`all_atom_md/md_steps.md`](./all_atom_md/md_steps.md)
  - [`equilibrated_systems/README.md`](./equilibrated_systems/README.md)
  - [`analysis/README.md`](./analysis/README.md)
