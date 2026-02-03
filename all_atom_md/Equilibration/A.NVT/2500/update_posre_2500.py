#!/usr/bin/env python3
"""
Update `posre_*.itp` files by changing fx, fy, fz values from "1000" to "2500".
Only the numeric values are replaced; original spacing/tabs between columns
are preserved exactly.

Helper / usage:
- Run this script from the directory that contains your `posre_*.itp` files.
- Command: `python3 update_posre_2500.py`
- It will validate the number of matching files before proceeding.

Customize for your system:
- Update `EXPECTED_FILES` to the number of `posre_*.itp` files in your system.
  This is system-dependent (for example, 45 in this repository).

Behavior:
- Preserves comments, empty lines, and whitespace formatting.
- Prints a per-file summary of replacements.
"""
import re
import glob
import sys

# Values to replace
REPLACE_FROM = "1000"
REPLACE_TO = "2500"
EXPECTED_FILES = 45

# Data line: "atom  type   fx     fy     fz   [; optional comment]"
# Capture whitespace as-is to rebuild the line identically,
# changing only the three numeric values.
data_line = re.compile(
    r'^(\s*\d+\s+\d+\s+)'        # g1: spazi + atom + spazi + type + spazi
    r'(-?\d+(?:\.\d+)?)'         # g2: fx (numero)
    r'(\s+)'                     # g3: spazi tra fx e fy
    r'(-?\d+(?:\.\d+)?)'         # g4: fy (numero)
    r'(\s+)'                     # g5: spazi tra fy e fz
    r'(-?\d+(?:\.\d+)?)'         # g6: fz (numero)
    r'(\s*)'                     # g7: eventuali spazi prima del commento o fine riga
    r'(;.*)?$'                   # g8: commento opzionale
)

def main() -> int:
    files = sorted(glob.glob("posre_*.itp"))
    if len(files) != EXPECTED_FILES:
        print(f"[ERROR] Found {len(files)} files matching 'posre_*.itp' (expected {EXPECTED_FILES}).")
        print("Run this script in the correct directory or verify the file names, then retry.")
        return 2

    total_replacements = 0
    summary = []

    for fn in files:
        # Read while preserving original line endings (CR/LF)
        with open(fn, "r", encoding="utf-8", newline="") as f:
            lines = f.read().splitlines(True)  # keepends=True

        new_lines = []
        changed_here = 0

        for line in lines:
            # Strip only the line ending for matching; reattach it unchanged
            stripped = line.rstrip("\r\n")
            eol = line[len(stripped):]  # terminatore originale, p.es. "\n" o "\r\n"

            m = data_line.match(stripped)
            if not m:
                new_lines.append(line)
                continue

            g1, fx, g3, fy, g5, fz, g7, comment = m.groups()

            # Replace ONLY where the value is exactly "1000"
            nfx = REPLACE_TO if fx == REPLACE_FROM else fx
            nfy = REPLACE_TO if fy == REPLACE_FROM else fy
            nfz = REPLACE_TO if fz == REPLACE_FROM else fz

            # Count how many replacements happened on this line
            if nfx != fx: changed_here += 1
            if nfy != fy: changed_here += 1
            if nfz != fz: changed_here += 1

            rebuilt = f"{g1}{nfx}{g3}{nfy}{g5}{nfz}{g7}"
            if comment is not None:
                rebuilt += comment

            new_lines.append(rebuilt + eol)

        # Write the file only if at least one line changed
        if changed_here > 0:
            with open(fn, "w", encoding="utf-8", newline="") as f:
                f.writelines(new_lines)

        total_replacements += changed_here
        summary.append((fn, changed_here))

    print("Replacement complete.")
    print(f"Files processed: {len(files)} (expected {EXPECTED_FILES}).")
    print(f"Total values replaced: {total_replacements}.")
    for fn, c in summary:
        print(f" - {fn}: {c} replacements")

    return 0

if __name__ == "__main__":
    sys.exit(main())
