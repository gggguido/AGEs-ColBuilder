# AGEs Candidate Positions Workflow

This folder provides a ready-to-run workflow to:
1. generate a periodic fibril representation by expanding a full-length collagen triple helix into symmetry-related neighbors of the crystallographic unit cell;
2. identify candidate crosslink positions from inter-chain residue distances.

## Files

- `generate_crystalcontact_copies.py`
- `calculate_dist_CA_CA_lys_arg.py` (for glucosepane and pentosidine)
- `calculate_dist_CA_CA_lys_lys.py` (for MOLD)

## Step 1: Generate Crystal-Contact Copies (Periodic Representation)

Run with UCSF Chimera from this folder:

```bash
chimera --nogui --script "generate_crystalcontact_copies.py"
```

What this does:
- opens the input triple-helix PDB;
- expands symmetry-related neighbors with `crystalcontacts`;
- merges all generated copies into one output PDB.

Default input/output names are:
- Input: `Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb`
- Output: `unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb`

## Step 2A: Candidate Sites for Glucosepane and Pentosidine (LYS-ARG, CA-CA)

Run:

```bash
python3 calculate_dist_CA_CA_lys_arg.py
```

This script:
- measures CA-CA distances for LYS-ARG candidate pairs;
- reorders pairs as `ARG -> LYS`;
- removes duplicates internally (ignoring molecule index);
- writes one final deduplicated output file.

Output:
- `lys_arg_CA_distances.csv`

## Step 2B: Candidate Sites for MOLD (LYS-LYS, CA-CA)

Run:

```bash
python3 calculate_dist_CA_CA_lys_lys.py
```

This script:
- measures CA-CA distances for LYS-LYS candidate pairs;
- removes duplicates internally:
  - unique key = `("Chain 1", "LYS1 Residue", "Chain 2", "LYS2 Residue")`
  - molecule index is ignored during deduplication;
- writes one final deduplicated output file.

Output:
- `lys_lys_CA_distances.csv`

## Notes

- Distance threshold default is 15.0 Angstrom and can be changed for either script with:

```bash
python3 calculate_dist_CA_CA_lys_arg.py --distance-threshold 12.0
python3 calculate_dist_CA_CA_lys_lys.py --distance-threshold 12.0
```

- Number of workers can be changed with:

```bash
python3 calculate_dist_CA_CA_lys_arg.py --processes 8
python3 calculate_dist_CA_CA_lys_lys.py --processes 8
```
