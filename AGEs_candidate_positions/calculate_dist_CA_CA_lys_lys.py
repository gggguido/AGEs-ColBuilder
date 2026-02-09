#!/usr/bin/env python3
"""
Calculate inter-molecule LYS/LYS CA-CA distances and write deduplicated candidates.

Helper
Pipeline order:
1) Generate crystal-contact copies (Chimera):
   chimera --nogui --script "generate_crystalcontact_copies.py"

2) Calculate LYS-LYS candidate distances:
   python3 calculate_dist_CA_CA_lys_lys.py

Default input/output names are aligned with generate_crystalcontact_copies.py:
- Input PDB:
  unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb
- Output CSV (already deduplicated):
  lys_lys_CA_distances.csv

Deduplication logic
- Matches the previous LYS-LYS duplicate filter rule:
  unique key = ("Chain 1", "LYS1 Residue", "Chain 2", "LYS2 Residue")
- Molecule indices are ignored during deduplication.
"""

import argparse
import csv
import multiprocessing as mp

import MDAnalysis as mda
import numpy as np


DEFAULT_INPUT_PDB = "unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb"
DEFAULT_OUTPUT_CSV = "lys_lys_CA_distances.csv"

REGIONS = {
    "M1": "resid 1:234",
    "M2": "resid 235:468",
    "M3": "resid 469:702",
    "M4": "resid 703:936",
    "M5": "resid 937:1056",
}

REGION_PAIRS = [
    ("M1", "M5"),
    ("M5", "M4"),
    ("M4", "M3"),
    ("M3", "M2"),
    ("M2", "M1"),
]

CSV_HEADER = [
    "Region1",
    "Region2",
    "Molecule 1",
    "Chain 1",
    "LYS1 Residue",
    "Resname 1",
    "Molecule 2",
    "Chain 2",
    "LYS2 Residue",
    "Resname 2",
    "Distance (CA-CA)",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate LYS/LYS CA-CA contacts across crystal-contact copies with "
            "integrated deduplication."
        )
    )
    parser.add_argument(
        "--input-pdb",
        default=DEFAULT_INPUT_PDB,
        help="Input merged crystal-contact PDB (default: %(default)s)",
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=15.0,
        help="Maximum CA-CA distance in Angstrom (default: %(default)s)",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=20,
        help="Multiprocessing workers (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_CSV,
        help="Output CSV path (deduplicated LYS-LYS rows) (default: %(default)s)",
    )
    return parser.parse_args()


def split_molecules_correctly(universe):
    molecules = []
    current = []
    prev_chain = None

    for atom in universe.atoms:
        if atom.segid == "A" and atom.resid == 1 and prev_chain == "C":
            if current:
                molecules.append(universe.atoms[current])
            current = []

        current.append(atom.index)
        prev_chain = atom.segid

    if current:
        molecules.append(universe.atoms[current])

    return molecules


def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1.position - atom2.position)


def _pairwise_rows(
    lys_residues_a,
    lys_residues_b,
    region_a,
    region_b,
    mol_a_idx,
    chain_a,
    mol_b_idx,
    chain_b,
    distance_threshold,
):
    out = []
    for res_a in lys_residues_a:
        ca_a = res_a.atoms.select_atoms("name CA")
        if len(ca_a) != 1:
            continue

        for res_b in lys_residues_b:
            ca_b = res_b.atoms.select_atoms("name CA")
            if len(ca_b) != 1:
                continue

            dist = calculate_distance(ca_a[0], ca_b[0])
            if dist <= distance_threshold:
                out.append(
                    [
                        region_a,
                        region_b,
                        mol_a_idx + 1,
                        chain_a,
                        res_a.resid,
                        res_a.resname,
                        mol_b_idx + 1,
                        chain_b,
                        res_b.resid,
                        res_b.resname,
                        "{:.2f}".format(dist),
                    ]
                )

    return out


def calculate_distances_for_region_pair(job):
    (
        mol1,
        mol2,
        mol1_idx,
        mol2_idx,
        region1,
        region2,
        distance_threshold,
    ) = job

    chains_mol1 = {seg: mol1.select_atoms("segid {}".format(seg)) for seg in "ABC"}
    chains_mol2 = {seg: mol2.select_atoms("segid {}".format(seg)) for seg in "ABC"}

    rows = []
    for ch1_name, ch1 in chains_mol1.items():
        for ch2_name, ch2 in chains_mol2.items():
            lys_1 = ch1.select_atoms(
                "resname LYS and {}".format(REGIONS[region1])
            ).residues
            lys_2 = ch2.select_atoms(
                "resname LYS and {}".format(REGIONS[region2])
            ).residues

            rows.extend(
                _pairwise_rows(
                    lys_1,
                    lys_2,
                    region1,
                    region2,
                    mol1_idx,
                    ch1_name,
                    mol2_idx,
                    ch2_name,
                    distance_threshold,
                )
            )

    return rows


def deduplicate_lys_lys_rows(rows):
    """
    Dedup rule matching the previous LYS-LYS filter script:
    subset = ["Chain 1", "LYS1 Residue", "Chain 2", "LYS2 Residue"]
    """
    unique = []
    seen = set()

    for row in rows:
        key = (
            row[3],  # Chain 1
            int(row[4]),  # LYS1 Residue
            row[7],  # Chain 2
            int(row[8]),  # LYS2 Residue
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)

    return unique


def write_rows(path, rows):
    with open(path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(CSV_HEADER)
        writer.writerows(rows)


def main():
    args = parse_args()

    universe = mda.Universe(args.input_pdb)
    molecules = split_molecules_correctly(universe)

    molecule_pairs = [
        (molecules[i], molecules[j], i, j, r1, r2, args.distance_threshold)
        for i in range(len(molecules) - 1)
        for j in range(i + 1, len(molecules))
        for r1, r2 in REGION_PAIRS
    ]

    with mp.Pool(processes=args.processes) as pool:
        all_groups = pool.map(calculate_distances_for_region_pair, molecule_pairs)

    raw_rows = [row for group in all_groups for row in group]
    unique_rows = deduplicate_lys_lys_rows(raw_rows)
    write_rows(args.output, unique_rows)

    print("Raw rows: {}".format(len(raw_rows)))
    print("Unique LYS-LYS rows: {}".format(len(unique_rows)))
    print("Output CSV written: {}".format(args.output))


if __name__ == "__main__":
    main()
