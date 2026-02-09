#!/usr/bin/env python3
"""
Calculate inter-molecule ARG/LYS CA-CA distances and write deduplicated candidates.

Helper
Pipeline order:
1) Generate crystal-contact copies (Chimera):
   chimera --nogui --script "generate_crystalcontact_copies.py"

2) Calculate distances:
   python3 calculate_dist_CA_CA_lys_arg.py

Default input/output names are aligned with generate_crystalcontact_copies.py:
- Input PDB:
  unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb
- Output CSV (already deduplicated):
  lys_arg_CA_distances.csv
"""

import argparse
import csv
import multiprocessing as mp

import MDAnalysis as mda
import numpy as np


DEFAULT_INPUT_PDB = "unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb"
DEFAULT_OUTPUT_CSV = "lys_arg_CA_distances.csv"

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
    "Resid 1",
    "Resname 1",
    "Molecule 2",
    "Chain 2",
    "Resid 2",
    "Resname 2",
    "Distance (CA-CA)",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate ARG/LYS CA-CA contacts across crystal-contact copies with "
            "integrated reordering and deduplication."
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
        help="Output CSV path (deduplicated ARG->LYS rows) (default: %(default)s)",
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
    residues_a,
    residues_b,
    region_a,
    region_b,
    mol_a_idx,
    chain_a,
    mol_b_idx,
    chain_b,
    distance_threshold,
):
    out = []
    for res_a in residues_a:
        ca_a = res_a.atoms.select_atoms("name CA")
        if len(ca_a) != 1:
            continue

        for res_b in residues_b:
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
            arg_1 = ch1.select_atoms(
                "resname ARG and {}".format(REGIONS[region1])
            ).residues

            lys_2 = ch2.select_atoms(
                "resname LYS and {}".format(REGIONS[region2])
            ).residues
            arg_2 = ch2.select_atoms(
                "resname ARG and {}".format(REGIONS[region2])
            ).residues

            # LYS on mol1 vs ARG on mol2
            rows.extend(
                _pairwise_rows(
                    lys_1,
                    arg_2,
                    region1,
                    region2,
                    mol1_idx,
                    ch1_name,
                    mol2_idx,
                    ch2_name,
                    distance_threshold,
                )
            )

            # ARG on mol1 vs LYS on mol2
            rows.extend(
                _pairwise_rows(
                    arg_1,
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


def orient_arg_lys_rows(rows):
    """
    Orient each pair as ARG -> LYS.
    Equivalent to the swap logic in remove_crystal_duplicates.py.
    """
    oriented = []
    bad_rows = []

    for row in rows:
        # row schema follows CSV_HEADER
        if row[5] == "LYS" and row[9] == "ARG":
            row = [
                row[1],
                row[0],
                row[6],
                row[7],
                row[8],
                row[9],
                row[2],
                row[3],
                row[4],
                row[5],
                row[10],
            ]

        if row[5] == "ARG" and row[9] == "LYS":
            oriented.append(row)
        else:
            bad_rows.append(row)

    if bad_rows:
        raise RuntimeError(
            "{} row(s) are not ARG->LYS after reordering.".format(len(bad_rows))
        )

    return oriented


def deduplicate_oriented_rows(rows):
    """
    Remove oriented duplicates ignoring molecule index, matching remove_crystal_duplicates.py.
    """
    unique = []
    seen = set()

    for row in rows:
        key = (
            row[3],  # Chain 1
            int(row[4]),  # Resid 1
            row[5],  # Resname 1
            row[7],  # Chain 2
            int(row[8]),  # Resid 2
            row[9],  # Resname 2
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

    oriented = orient_arg_lys_rows(raw_rows)
    unique_rows = deduplicate_oriented_rows(oriented)
    write_rows(args.output, unique_rows)

    print("Raw rows: {}".format(len(raw_rows)))
    print("Unique ARG->LYS rows: {}".format(len(unique_rows)))
    print("Output CSV written: {}".format(args.output))


if __name__ == "__main__":
    main()
