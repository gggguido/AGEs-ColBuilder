#!/usr/bin/env python3
"""
Update `posre_*.itp` files by changing fx, fy, fz values from "1000" to "300"
while preserving the original column alignment. When a value goes from 4 digits
(1000) to 3 digits (300), the script pads the following whitespace so the next
column starts at the exact same position as before, cascading through comments.

Helper / usage:
- Run this script from the directory that contains your `posre_*.itp` files.
- Command: `python3 update_posre_300.py`
- It will validate the number of matching files before proceeding.

Customize for your system:
- Update `EXPECTED_FILES` to the number of `posre_*.itp` files in your system.
  This is system-dependent (for example, 45 in this repository).

Behavior:
- Preserves comments, empty lines, and line endings.
- Preserves tabs in separators; only adds spaces at the end of whitespace runs.
- Prints a per-file summary of replacements.
"""
import re
import glob
import sys

FROM = "1000"
TO = "300"
EXPECTED_FILES = 45

# Match a data line from the posre block: "atom  type   fx  fy  fz  [; comment]"
data_line = re.compile(
    r'^(\s*\d+\s+\d+\s+)'        # g1: spazi + 'atom' + spazi + 'type' + spazi fino a FX
    r'(-?\d+(?:\.\d+)?)'         # g2: fx
    r'(\s+)'                     # g3: spazi tra fx e fy
    r'(-?\d+(?:\.\d+)?)'         # g4: fy
    r'(\s+)'                     # g5: spazi tra fy e fz
    r'(-?\d+(?:\.\d+)?)'         # g6: fz
    r'(\s*)'                     # g7: spazi prima di eventuale commento
    r'(;.*)?$'                   # g8: commento opzionale
)

def pad_spaces_preserving_tabs(original_ws: str, new_len: int) -> str:
    """
    Extend 'original_ws' to 'new_len' by adding ONLY SPACES at the end.
    If new_len is smaller or equal, truncate to new_len characters (rare case).
    """
    cur = len(original_ws)
    if new_len <= cur:
        return original_ws[:new_len]
    # Add only trailing spaces to reach the desired length
    return original_ws + (" " * (new_len - cur))

def main() -> int:
    files = sorted(glob.glob("posre_*.itp"))
    if len(files) != EXPECTED_FILES:
        print(f"[ERROR] Found {len(files)} files matching 'posre_*.itp' (expected {EXPECTED_FILES}).")
        print("Run this script in the correct directory or verify the file names, then retry.")
        return 2

    total_replacements = 0
    summary = []

    for fn in files:
        with open(fn, "r", encoding="utf-8", newline="") as f:
            lines = f.read().splitlines(True)  # preserva terminatori

        new_lines = []
        changed_here = 0

        for line in lines:
            body = line.rstrip("\r\n")
            eol = line[len(body):]

            m = data_line.match(body)
            if not m:
                new_lines.append(line)
                continue

            g1, fx, g3, fy, g5, fz, g7, comment = m.groups()

            # ORIGINAL column start positions (relative to the original 'body')
            orig_fy_start = len(g1) + len(fx) + len(g3)
            orig_fz_start = orig_fy_start + len(fy) + len(g5)
            orig_cmt_start = None if comment is None else (orig_fz_start + len(fz) + len(g7))

            # New values (only if they were exactly "1000")
            nfx = TO if fx == FROM else fx
            nfy = TO if fy == FROM else fy
            nfz = TO if fz == FROM else fz

            if nfx != fx: changed_here += 1
            if nfy != fy: changed_here += 1
            if nfz != fz: changed_here += 1

            # Recompute whitespace to keep the same column start positions
            # New g3 length so the new fy starts where it used to
            new_g3_len = max(1, orig_fy_start - (len(g1) + len(nfx)))
            new_g3 = pad_spaces_preserving_tabs(g3, new_g3_len)

            # New g5 length so the new fz starts where it used to
            new_fy_start = len(g1) + len(nfx) + len(new_g3)
            new_g5_len = max(1, orig_fz_start - (new_fy_start + len(nfy)))
            new_g5 = pad_spaces_preserving_tabs(g5, new_g5_len)

            # Keep the comment start (if present) identical
            if comment is not None:
                new_fz_start = new_fy_start + len(nfy) + len(new_g5)
                new_g7_len = max(0, orig_cmt_start - (new_fz_start + len(nfz)))
                new_g7 = pad_spaces_preserving_tabs(g7, new_g7_len)
            else:
                new_g7 = g7  # no comment to realign

            rebuilt = f"{g1}{nfx}{new_g3}{nfy}{new_g5}{nfz}{new_g7}"
            if comment is not None:
                rebuilt += comment

            new_lines.append(rebuilt + eol)

        # Write only if at least one replacement occurred in this file
        if changed_here > 0:
            with open(fn, "w", encoding="utf-8", newline="") as f:
                f.writelines(new_lines)

        total_replacements += changed_here
        summary.append((fn, changed_here))

    print("Replacement complete (1000 -> 300) with columns preserved.")
    print(f"Files processed: {len(files)} (expected {EXPECTED_FILES}).")
    print(f"Total values replaced: {total_replacements}.")
    for fn, c in summary:
        print(f" - {fn}: {c} replacements")

    return 0

if __name__ == "__main__":
    sys.exit(main())
