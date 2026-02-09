# AGEs-ColBuilder Extension Repository

This repository collects scripts, inputs, and example systems to help ColBuilder users include and validate new collagen crosslink chemistries that are not yet implemented by default.

## Author

- **Guido Giannetti**

## Publication

This work is currently in preparation.

- **Citation:** _in preparation_
- **DOI placeholder:** `DOI: TBA`

## Repository Structure

### `AGEs_candidate_positions`
Tools to generate periodic crystal-contact copies of a collagen triple helix and identify candidate crosslink sites with geometric distance filters.

- `generate_crystalcontact_copies.py` expands the crystallographic environment (via Chimera crystal contacts).
- `calculate_dist_CA_CA_lys_arg.py` finds **LYS-ARG** candidates (used for glucosepane and pentosidine).
- `calculate_dist_CA_CA_lys_lys.py` finds **LYS-LYS** candidates (used for MOLD).

### `MODELLER`
Scripts to generate internal-coordinate (IC) tables for new crosslink residues from SMILES geometry.

- `create_ic_*` scripts produce `.str` files containing IC entries.
- The generated IC entries are then transferred into the ColBuilder modeller library (`top_heav_mod.lib`) to enable new crosslink building.

The workflow is general and can be adapted to additional crosslink chemistries beyond the examples included here.

### `parametrization`
Detailed example workflow for parametrizing a new crosslink molecule and obtaining an **Amber99-compatible** parameter set.

- Includes a complete step-by-step README and generated files for **glucosepane**.

### `all_atom_md`
Detailed all-atom MD protocol and inputs used after generating collagen fibrils with ColBuilder.

- `md_steps.md` provides a complete list of GROMACS commands used in this work.
- Includes the `.mdp` files used for minimization, equilibration, and production pulling runs.

### `equilibrated_systems`
Archived equilibrated systems used in this work.

- Contains zipped equilibrated `.gro` systems for:
  - glucosepane
  - pentosidine
  - MOLD
  - PYD
  - PYD + glucosepane

### `analysis`
Post-processing scripts used to analyze fibrils under load.

The scripts compute:
- D-band length proxy (end-to-end extension)
- overlap length
- gap length
- overlap/gap ratio
- associated uncertainties

## Suggested Usage Order

1. Identify candidate residue pairs in `AGEs_candidate_positions`.
2. Build crosslink IC definitions in `MODELLER` and update ColBuilder modeller library.
3. Parametrize the new chemistry in `parametrization` (if needed).
4. Build fibrils with ColBuilder and run MD following `all_atom_md`.
5. Use `equilibrated_systems` examples as references/checkpoints.
6. Analyze outputs with scripts in `analysis`.

## Notes

- Each subfolder contains its own detailed documentation (`README.md` or `md_steps.md`).
- File naming and command examples are aligned across folders to support a direct end-to-end workflow.
