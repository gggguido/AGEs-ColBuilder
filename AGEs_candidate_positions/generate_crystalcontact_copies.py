#!/usr/bin/env python
"""
Generate crystal-contact copies in UCSF Chimera and merge them into one PDB.

Helper
- Default output is aligned with:
  - `calculate_dist_CA_CA_lys_arg.py`
  - `calculate_dist_CA_CA_lys_lys.py`

Launch examples (from this folder):
1) Defaults (input + output already aligned to distance script):
   chimera --nogui --script "generate_crystalcontact_copies.py"

2) Custom input/output:
   chimera --nogui --script "generate_crystalcontact_copies.py --input your_model.pdb --output unit_cell_crystalcontacts_your_model.pdb"

Important
- Run this script with UCSF Chimera (not plain python).
"""

import argparse
import os

import chimera
from chimera import runCommand


DEFAULT_INPUT = "Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb"
DEFAULT_OUTPUT = "unit_cell_crystalcontacts_Rattus_norvegicus_aln_N_NOCROSS_C_NOCROSS.pdb"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate crystal-contact copies and write one merged PDB.")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input PDB (default: %(default)s)")
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=(
            "Merged output PDB. Default is aligned with both distance scripts: "
            "%(default)s"
        ),
    )
    parser.add_argument(
        "--contact-distance",
        type=float,
        default=26.0,
        help="Distance cutoff passed to Chimera crystalcontacts (default: %(default)s)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_pdb = args.input
    output_pdb = args.output

    print("Debug: Using input PDB file: {}".format(input_pdb))
    if not os.path.exists(input_pdb):
        print("Error: Specified PDB file '{}' does not exist".format(input_pdb))
        return

    try:
        chimera.openModels.open(input_pdb)
        print("Debug: Successfully opened PDB file: {}".format(input_pdb))
    except Exception as exc:
        print("Error: Failed to open PDB file '{}' with exception: {}".format(input_pdb, exc))
        return

    print("Debug: Generating crystal contacts for the unit cell")
    runCommand(
        "crystalcontacts #0 {} copies true schematic false".format(args.contact_distance)
    )

    models = sorted(
        chimera.openModels.list(modelTypes=[chimera.Molecule]),
        key=lambda model: model.id,
    )
    print("Debug: Models after crystalcontacts: {}".format(models))

    generated_pdbs = []
    for model in models:
        pdb_filename = "tmp_copy_{}.pdb".format(model.id)
        runCommand("write format pdb #{} {}".format(model.id, pdb_filename))
        generated_pdbs.append(pdb_filename)
        print("Debug: Saved model #{} as {}".format(model.id, pdb_filename))

    cryst1_line = "CRYST1   39.970   26.950  677.900  89.24  94.59 105.58 P 1           2\n"
    try:
        with open(output_pdb, "w") as combined_file:
            combined_file.write(cryst1_line)
            for pdb_filename in generated_pdbs:
                with open(pdb_filename, "r") as pdb_file:
                    for line in pdb_file:
                        if line.startswith(("ATOM", "HETATM")):
                            combined_file.write(line)
                    combined_file.write("TER\n")
            combined_file.write("END\n")
        print("Debug: Final combined PDB saved as {}".format(output_pdb))
    except Exception as exc:
        print(
            "Error: Failed to combine PDB files into '{}' with exception: {}".format(
                output_pdb, exc
            )
        )

    for pdb_filename in generated_pdbs:
        try:
            os.remove(pdb_filename)
            print("Debug: Removed temporary file '{}'".format(pdb_filename))
        except Exception as exc:
            print(
                "Warning: Could not remove temporary file '{}' with exception: {}".format(
                    pdb_filename, exc
                )
            )

    runCommand("close all")
    print("Script completed")


if __name__ == "__main__":
    main()
