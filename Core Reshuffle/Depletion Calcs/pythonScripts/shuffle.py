"""
reshuffle_core.py
=================
Reshuffles the assembly burn-up data in step 26 (the final time step) of a
SIMULATE-style depletion file (.dep).

Usage
-----
Run the script and pass the input .dep file as the first argument:

    python reshuffle_core.py depletion.dep [output.dep]

The script will then interactively ask you for the new grid layout row by row.

Alternatively, set NEW_GRID_LAYOUT (a list-of-lists) near the bottom of this
script to skip the interactive prompt and run non-interactively.

How it works
------------
1.  Reads the core layout from the file header (assembly numbers 1-193).
2.  Reads the 2-D "Assembly Burnup Distribution" map from step 26.
3.  Reads the 3-D "EXP 3D MAP" (per-axial-layer burnup, lfa-indexed) from step 26.
4.  Accepts a NEW_GRID with the same shape as the header layout, where each
    cell contains either:
      - An assembly number (1-193)  →  take that assembly's burnup from step 26.
      - 0                           →  fresh assembly, burnup = 0.150 everywhere.
5.  Rebuilds both the 2-D map and the 3-D EXP map using reshuffled burnup values,
    preserving all original formatting/spacing exactly.
6.  Writes the modified file (original filename + "_reshuffled" suffix by default).

Core layout (from file header)
-------------------------------
Row  1:                               1   2   3   4   5   6   7               (7)
Row  2:              8   9  10  11  12  13  14  15  16  17  18                (11)
Row  3:        19  20  21  22  23  24  25  26  27  28  29  30  31             (13)
Row  4:        32  33  34  35  36  37  38  39  40  41  42  43  44             (13)
Row  5:    45  46  47  48  49  50  51  52  53  54  55  56  57  58  59         (15)
Row  6:    60  61  62  63  64  65  66  67  68  69  70  71  72  73  74         (15)
Row  7:    75  76  77  78  79  80  81  82  83  84  85  86  87  88  89         (15)
Row  8:    90  91  92  93  94  95  96  97  98  99 100 101 102 103 104         (15)
Row  9:   105 106 107 108 109 110 111 112 113 114 115 116 117 118 119         (15)
Row 10:   120 121 122 123 124 125 126 127 128 129 130 131 132 133 134         (15)
Row 11:   135 136 137 138 139 140 141 142 143 144 145 146 147 148 149         (15)
Row 12:       150 151 152 153 154 155 156 157 158 159 160 161 162             (13)
Row 13:       163 164 165 166 167 168 169 170 171 172 173 174 175             (13)
Row 14:           176 177 178 179 180 181 182 183 184 185 186                 (11)
Row 15:                   187 188 189 190 191 192 193                         (7)
"""

import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Core geometry
# ---------------------------------------------------------------------------
CORE_ROWS = [
    [1, 2, 3, 4, 5, 6, 7],
    [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
    [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
    [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44],
    [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
    [60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74],
    [75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
    [90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104],
    [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119],
    [120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134],
    [135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149],
    [150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162],
    [163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175],
    [176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186],
    [187, 188, 189, 190, 191, 192, 193],
]

TOTAL_ASSEMBLIES = 193
FRESH_BURNUP     = 0.150   # GWd/T for a fresh assembly


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def fmt_num(value, decimal_places, field_width):
    """Return value formatted to decimal_places, right-justified in field_width chars."""
    s = f"{value:.{decimal_places}f}"
    return s.rjust(field_width)


def parse_data_fields(rest_str):
    """
    Parse the value portion of a data line (after the k-row prefix).
    Returns a list of (leading_spaces: str, number_str: str) tuples that
    together reconstruct rest_str exactly.
    """
    result = []
    for m in re.finditer(r'( *)([\d.]+)', rest_str):
        result.append((m.group(1), m.group(2)))
    return result


def rebuild_data_line(prefix, fields, new_values):
    """
    Rebuild a data line using new_values, preserving the original column
    widths and decimal precision implied by the original fields.

    prefix      : the k-row prefix string, e.g. "  21    "
    fields      : list of (leading_spaces, original_number_str) from parse_data_fields
    new_values  : list of float values (same length as fields)
    """
    parts = []
    for (spaces, orig_str), v in zip(fields, new_values):
        dp = len(orig_str.split('.')[1]) if '.' in orig_str else 0
        fw = len(orig_str)
        parts.append(spaces + fmt_num(v, dp, fw))
    return prefix + ''.join(parts)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def find_nth_begin_step(lines, n):
    """Return 0-based line index of the n-th BEGIN STEP (1-based n)."""
    count = 0
    for i, line in enumerate(lines):
        if "BEGIN STEP" in line:
            count += 1
            if count == n:
                return i
    raise ValueError(f"Step {n} not found.")


def parse_2d_burnup(lines, section_start):
    """
    Parse the 'Assembly Burnup Distribution' 2-D section.

    Returns
    -------
    burnup_2d : dict  assembly_number -> float (average burnup)
    end_idx   : index of the first line AFTER the 15 data rows
    """
    burnup_2d = {}
    idx = section_start + 1  # skip header line
    for core_row in CORE_ROWS:
        raw = lines[idx].rstrip("\r\n")
        values = [float(v) for v in raw.split()]
        if len(values) != len(core_row):
            raise ValueError(
                f"2-D burnup: row length mismatch (expected {len(core_row)}, "
                f"got {len(values)}) at file line {idx+1}"
            )
        for asm, val in zip(core_row, values):
            burnup_2d[asm] = val
        idx += 1
    return burnup_2d, idx


def parse_3d_exp_map(lines, section_start):
    """
    Parse the 'EXP 3D MAP' section.

    Returns
    -------
    exp_3d      : dict  lfa -> list[float], k=21 down to k=2 (20 layers)
    block_data  : list of (line_idx, is_data, prefix, fields, lfa_columns)
                  is_data=True → data row; prefix, fields, lfa_columns populated
                  is_data=False → verbatim line (header / blank)
    end_idx     : first index after the section's data
    """
    exp_3d = {lfa: [] for lfa in range(1, TOTAL_ASSEMBLIES + 1)}
    block_data = []
    lfa_cols = []

    idx = section_start + 1  # skip " EXP 3D MAP 1.0E+00" line

    while idx < len(lines):
        line = lines[idx]
        stripped = line.strip().rstrip("\r\n")

        if stripped == "":
            block_data.append((idx, False, line, None, None))
            idx += 1
            continue

        if stripped.startswith("k lfa"):
            parts = stripped.split()
            lfa_cols = [int(x) for x in parts[2:]]
            block_data.append((idx, False, line, None, None))
            idx += 1
            continue

        m = re.match(r'^(\s+\d+\s+)(.+?)[\r\n]', line)
        if m:
            prefix = m.group(1)
            rest   = m.group(2)
            fields = parse_data_fields(rest)
            values = [float(n) for _, n in fields]

            if len(values) != len(lfa_cols):
                raise ValueError(
                    f"EXP 3D MAP: column count mismatch at line {idx+1} "
                    f"(got {len(values)}, expected {len(lfa_cols)})"
                )
            for lfa, val in zip(lfa_cols, values):
                exp_3d[lfa].append(val)

            block_data.append((idx, True, prefix, fields, list(lfa_cols)))
            idx += 1
            continue

        # Anything else (END STEP, etc.) means we're done
        break

    return exp_3d, block_data, idx


# ---------------------------------------------------------------------------
# Rebuilding
# ---------------------------------------------------------------------------

def build_2d_section(new_burnup, original_lines, section_start):
    """
    Rebuild the 2-D Assembly Burnup Distribution section.
    Preserves leading whitespace per row from the original.

    Returns (list_of_new_lines, idx_after_data).
    """
    new_lines = [original_lines[section_start]]  # keep "Assembly Burnup Distribution" header
    idx = section_start + 1

    for core_row in CORE_ROWS:
        orig = original_lines[idx].rstrip("\r\n")
        leading = len(orig) - len(orig.lstrip())
        prefix  = " " * leading
        values  = [new_burnup[asm] for asm in core_row]
        # Original format: 6-char values (X.XXX) separated by 2 spaces
        body = "  ".join(f"{v:6.3f}" for v in values)
        new_lines.append(prefix + body + "\r\n")
        idx += 1

    return new_lines, idx


def build_3d_section(new_exp_3d, block_data):
    """
    Rebuild the EXP 3D MAP data lines, substituting new burnup values while
    preserving exact column widths, spacing, and decimal precision of the original.

    For fresh assemblies (burnup 0.150) the original field structure is used as
    a template, with the value replaced by 0.150.
    """
    consumed = {lfa: 0 for lfa in range(1, TOTAL_ASSEMBLIES + 1)}
    rebuilt  = []

    for (line_idx, is_data, prefix_or_line, fields, lfa_cols) in block_data:
        if not is_data:
            rebuilt.append(prefix_or_line)
            continue

        new_vals = []
        for lfa in lfa_cols:
            c = consumed[lfa]
            new_vals.append(new_exp_3d[lfa][c])
            consumed[lfa] = c + 1

        new_line = rebuild_data_line(prefix_or_line, fields, new_vals) + "\r\n"
        rebuilt.append(new_line)

    return rebuilt


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def reshuffle(input_path, new_grid, step=26, output_path=None):
    """
    Parameters
    ----------
    input_path  : str   Path to the original .dep file.
    new_grid    : list  15 rows, each a list of ints matching CORE_ROWS lengths.
                        Values: assembly number (1-193) or 0 (fresh).
    step        : int   Which step to modify (default 26).
    output_path : str   Output path; defaults to input stem + "_reshuffled" + suffix.
    """
    # Validate grid shape
    if len(new_grid) != len(CORE_ROWS):
        raise ValueError(f"new_grid must have {len(CORE_ROWS)} rows, got {len(new_grid)}.")
    for i, (nr, cr) in enumerate(zip(new_grid, CORE_ROWS)):
        if len(nr) != len(cr):
            raise ValueError(f"Row {i+1}: expected {len(cr)} values, got {len(nr)}.")
        bad = [v for v in nr if not (v == 0 or 1 <= v <= TOTAL_ASSEMBLIES)]
        if bad:
            raise ValueError(f"Row {i+1}: invalid values {bad} (must be 0 or 1-{TOTAL_ASSEMBLIES}).")

    print(f"Reading '{input_path}' ...", flush=True)
    raw   = Path(input_path).read_bytes()
    lines = raw.decode("ascii", errors="replace").splitlines(keepends=True)
    print(f"  {len(lines)} lines read.")

    # Locate step
    step_start = find_nth_begin_step(lines, step)
    print(f"  Step {step} begins at line {step_start + 1}.")

    # Find 2-D burnup section
    burn2d_idx = None
    for i in range(step_start, len(lines)):
        if "Assembly Burnup Distribution" in lines[i]:
            burn2d_idx = i
            break
    if burn2d_idx is None:
        raise ValueError("'Assembly Burnup Distribution' not found in target step.")

    old_burnup_2d, after_2d = parse_2d_burnup(lines, burn2d_idx)
    print(f"  Parsed 2-D burnup ({len(old_burnup_2d)} assemblies).")

    # Find EXP 3D MAP section
    exp3d_idx = None
    for i in range(burn2d_idx, len(lines)):
        if lines[i].strip().startswith("EXP 3D MAP"):
            exp3d_idx = i
            break
    if exp3d_idx is None:
        raise ValueError("'EXP 3D MAP' not found in target step.")

    old_exp_3d, block_data, after_3d = parse_3d_exp_map(lines, exp3d_idx)
    n_layers = len(old_exp_3d[1])
    print(f"  Parsed 3-D EXP map ({n_layers} axial layers per assembly).")

    # Build reshuffled burnup maps
    new_burnup_2d = {}
    new_exp_3d    = {lfa: [] for lfa in range(1, TOTAL_ASSEMBLIES + 1)}

    for ng_row, cr_row in zip(new_grid, CORE_ROWS):
        for source_asm, dest_lfa in zip(ng_row, cr_row):
            if source_asm == 0:
                new_burnup_2d[dest_lfa] = FRESH_BURNUP
                new_exp_3d[dest_lfa]    = [FRESH_BURNUP] * n_layers
            else:
                new_burnup_2d[dest_lfa] = old_burnup_2d[source_asm]
                new_exp_3d[dest_lfa]    = list(old_exp_3d[source_asm])

    print("  Burnup values reshuffled.")

    # Rebuild sections
    new_2d_lines, after_new_2d = build_2d_section(new_burnup_2d, lines, burn2d_idx)
    new_3d_lines = build_3d_section(new_exp_3d, block_data)

    # Assemble output
    # Part A: everything before the 2-D burnup section
    # Part B: new 2-D lines
    # Part C: gap between 2-D section end and EXP 3D header (blank line + EXP header)
    # Part D: new 3-D lines
    # Part E: rest of file from after_3d onward
    part_c = lines[after_new_2d : exp3d_idx + 1]  # includes EXP header line
    out_lines = lines[:burn2d_idx] + new_2d_lines + part_c + new_3d_lines + lines[after_3d:]

    if output_path is None:
        p = Path(input_path)
        output_path = str(p.parent / (p.stem + "_reshuffled" + p.suffix))

    Path(output_path).write_bytes("".join(out_lines).encode("ascii"))
    print(f"  Output written to '{output_path}'  ({len(out_lines)} lines).")
    return output_path


# ---------------------------------------------------------------------------
# Interactive prompt
# ---------------------------------------------------------------------------

def print_core_layout():
    print("\n  Core layout — assembly numbers as in file header:")
    max_w = max(len(r) for r in CORE_ROWS)
    for row in CORE_ROWS:
        pad = (max_w - len(row)) // 2
        print("  " + "      " * pad + "  ".join(f"{n:3d}" for n in row))
    print()


def prompt_new_grid():
    print_core_layout()
    print("  Enter the new layout row by row.")
    print("  Each number = source assembly whose burnup goes to that position.")
    print("  Use 0 for a FRESH assembly (burnup set to 0.150).\n")

    new_grid = []
    for i, core_row in enumerate(CORE_ROWS):
        n = len(core_row)
        while True:
            raw = input(
                f"  Row {i+1:2d}  ({n} positions, "
                f"assemblies {core_row[0]}–{core_row[-1]}): "
            )
            parts = raw.strip().split()
            if len(parts) != n:
                print(f"    ✗  Need {n} values, got {len(parts)}. Try again.")
                continue
            try:
                vals = [int(x) for x in parts]
            except ValueError:
                print("    ✗  Non-integer value. Try again.")
                continue
            bad = [v for v in vals if not (v == 0 or 1 <= v <= TOTAL_ASSEMBLIES)]
            if bad:
                print(f"    ✗  Invalid numbers: {bad}. Must be 0 or 1-{TOTAL_ASSEMBLIES}.")
                continue
            new_grid.append(vals)
            break
    return new_grid


# ---------------------------------------------------------------------------
# Configuration — edit these to run non-interactively
# ---------------------------------------------------------------------------
#
# Set NEW_GRID_LAYOUT to a list-of-lists to skip the interactive prompt.
# Example (identity — no shuffling):
#
NEW_GRID_LAYOUT = [
     [10, 8, 21, 25, 29, 18, 16],
     [1, 12, 0, 0, 0, 0, 0,0, 0, 14, 7],
     [45, 0, 0, 37, 11, 3, 13, 5, 15, 39, 0, 0, 59],
     [76, 0, 23, 0, 35, 0, 83, 0, 41, 0, 27, 0, 88],
     [46, 0, 78, 0, 62, 0, 50, 0, 54, 0, 72, 0, 86, 0, 58],
     [19, 0, 61, 48, 0, 94, 9, 30, 17, 52, 0, 56, 72, 0, 31],
     [33, 0, 75, 0, 64, 32, 96, 0, 82, 44, 70, 0, 89, 0, 43],
     [92, 0, 91, 81, 0, 20, 0, 22, 0, 174, 0, 113, 103, 0, 102],
     [151, 0, 105, 0, 124, 150, 112, 0, 98, 162, 130, 0, 119, 0, 161],
     [163, 0, 121, 138, 0, 142, 177, 164, 185, 100, 0, 146, 116, 0, 175],
     [136, 0, 108, 0, 122, 0, 140, 0, 144, 0, 132, 0, 133, 0, 148],
     [106, 0, 167, 0, 153, 0, 111, 0, 159, 0, 171, 0, 118],
     [135, 0, 0, 155, 179, 189, 181, 191, 183, 157, 0, 0, 149],
     [187, 180, 0, 0, 0, 0, 0, 0, 0, 182, 193],
     [178, 176, 165, 169, 173, 186, 184],
 ]

#NEW_GRID_LAYOUT = []   # None → use interactive prompt
INPUT_FILE      = "depletion.dep"
OUTPUT_FILE     = None   # None → auto-name with "_reshuffled"
TARGET_STEP     = 26


if __name__ == "__main__":
    if len(sys.argv) > 1:
        INPUT_FILE = sys.argv[1]
    if len(sys.argv) > 2:
        OUTPUT_FILE = sys.argv[2]

    print(f"\nCore Reshuffler — modifying step {TARGET_STEP} of '{INPUT_FILE}'\n")

    if NEW_GRID_LAYOUT is not None:
        new_grid = NEW_GRID_LAYOUT
        print("Using hard-coded NEW_GRID_LAYOUT.")
    else:
        new_grid = prompt_new_grid()

    reshuffle(INPUT_FILE, new_grid, step=TARGET_STEP, output_path=OUTPUT_FILE)
    print("\nDone.")