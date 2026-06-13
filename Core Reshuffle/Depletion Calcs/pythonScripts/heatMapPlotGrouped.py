import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle

# ── Configuration ─────────────────────────────────────────────────────────────
DEP_FILE = "./reshuffle/reshuffle1.dep"
STEP     = 1          # Valid step from your file
SAVE     = './Images/BOLheatmapReshuffled.png'

# ── Core row offsets for the 193-FA PWR layout ────────────────────────────────
ROW_LAYOUT = [
    (5,  7),
    (3, 11),
    (2, 13),
    (2, 13),
    (1, 15),
    (1, 15),
    (1, 15),
    (1, 15),
    (1, 15),
    (1, 15),
    (1, 15),
    (2, 13),
    (2, 13),
    (3, 11),
    (5,  7),
]

# ── Parse ─────────────────────────────────────────────────────────────────────
def parse_step(filepath, step_target):
    nrows = len(ROW_LAYOUT)
    ncols = 19
    grid = np.full((nrows, ncols), np.nan)
    fa_grid = np.full((nrows, ncols), 0, dtype=int)
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = [line.strip() for line in f.readlines()]

    # 1. Extract Core Configuration Layout (Assembly IDs)
    config_start = None
    for i, line in enumerate(lines):
        if "core configuration" in line.lower():
            config_start = i + 1
            break
            
    if config_start is None:
        raise ValueError("Could not find 'core configuration' header in the file.")
        
    layout_lines = []
    for line in lines[config_start:]:
        if "=====" in line or "BEGIN STEP" in line:
            break
        layout_lines.append(line)
        
    fa_ids = [int(x) for x in re.findall(r'\b\d+\b', " ".join(layout_lines))]
    
    fa_idx = 0
    for r, (offset, count) in enumerate(ROW_LAYOUT):
        for c in range(offset, offset + count):
            if fa_idx < len(fa_ids):
                fa_grid[r, c] = fa_ids[fa_idx]
                fa_idx += 1

    # 2. Locate the specific Target Step block
    step_found = False
    target_idx = 0
    days, cycle_bu = 0.0, 0.0
    step_str = str(step_target)
    
    for i, line in enumerate(lines):
        if "PT RE" in line and "DAYs" in line:
            data_line = lines[i+1]
            if "---" in data_line:
                data_line = lines[i+2]
            
            tokens = data_line.split()
            if tokens and tokens[0] == step_str:
                step_found = True
                target_idx = i
                try:
                    days = float(tokens[2])
                    cycle_bu = float(tokens[7]) 
                except (IndexError, ValueError):
                    pass
                break
                
    if not step_found:
        raise ValueError(f"Could not find summary block for Step {step_target}.")

    # 3. Locate the 'Assembly Burnup Distribution' 2D map for this step
    burnup_start = None
    for i in range(target_idx, len(lines)):
        if "assembly burnup distribution" in lines[i].lower():
            burnup_start = i + 1
            break
            
    if burnup_start is None:
        raise ValueError(f"Could not find 'Assembly Burnup Distribution' map for Step {step_target}.")
        
    # 4. Read the 15 spatial lines of the core map
    burnup_lines = []
    for line in lines[burnup_start:]:
        if len(burnup_lines) == nrows:
            break
        if not line.strip(): 
            continue
        burnup_lines.append(line)
        
    # 5. Map the string rows back into the 2D grid
    for r, (offset, count) in enumerate(ROW_LAYOUT):
        row_vals = [float(x) for x in burnup_lines[r].split()]
        
        for idx, c in enumerate(range(offset, offset + count)):
            if idx < len(row_vals):
                grid[r, c] = row_vals[idx]

    return grid, fa_grid, days, cycle_bu

# ── Main Execution ────────────────────────────────────────────────────────────
grid, fa_grid, days, cycle_bu = parse_step(DEP_FILE, STEP)

nrows, ncols = grid.shape
valid_mask = ~np.isnan(grid)
valid_values = grid[valid_mask]
valid_fa_numbers = fa_grid[valid_mask]

# ── Grouping into 4 Categories via Quantiles ──────────────────────────────────
# This splits the assemblies evenly into 4 quartiles
q1 = np.percentile(valid_values, 25)
q2 = np.percentile(valid_values, 50)
q3 = np.percentile(valid_values, 75)

labels = np.zeros_like(valid_values, dtype=int)
labels[valid_values <= q1] = 0
labels[(valid_values > q1) & (valid_values <= q2)] = 1
labels[(valid_values > q2) & (valid_values <= q3)] = 2
labels[valid_values > q3] = 3

# Colors for the 4 categories (Blue, Light Blue, Orange, Red)
colors = ['#2c7bb6', '#abd9e9', '#fdae61', '#d7191c']  
category_colors = np.zeros(grid.shape, dtype=object)

group_data = {}
for cat in range(4):
    mask_in_cluster = (labels == cat)
    
    group_fas = valid_fa_numbers[mask_in_cluster]
    group_burnups = valid_values[mask_in_cluster]
    avg_bu = group_burnups.mean() if len(group_burnups) > 0 else 0.0
    
    group_data[cat] = {
        'avg': avg_bu,
        'fas': sorted(group_fas.tolist()),
        'color': colors[cat]
    }
    
    idx_coords = np.argwhere(valid_mask)[mask_in_cluster]
    for r, c in idx_coords:
        category_colors[r, c] = colors[cat]

# ── Terminal Output ───────────────────────────────────────────────────────────
print(f"\n=== Fuel Assembly Grouping for Step {STEP} ===")
names = ["Q1 (Lowest)", "Q2 (Low-Mid)", "Q3 (High-Mid)", "Q4 (Highest)"]
for cat in range(4):
    print(f"\n🔵 {names[cat]} (Avg: {group_data[cat]['avg']:.2f} GWd/T):")
    print(f"   Assemblies ({len(group_data[cat]['fas'])} total): {group_data[cat]['fas']}")
print("\n============================================\n")

# ── Plotting ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))
ax.set_facecolor("white")

for r in range(nrows):
    for c in range(ncols):
        if not valid_mask[r, c]:
            continue
            
        color_hex = category_colors[r, c]
        v = grid[r, c]
        
        # Draw discrete color block
        rect = Rectangle((c, nrows - 1 - r), 1, 1, facecolor=color_hex, edgecolor="white", lw=0.5)
        ax.add_patch(rect)
        
        # Text Contrast Calculations
        rgb = mcolors.to_rgb(color_hex)
        lum = 0.2126*rgb[0] + 0.7152*rgb[1] + 0.0722*rgb[2]
        tc = "#1a1a1a" if lum > 0.45 else "#f5f5f5"
        
        cx, cy = c + 0.5, nrows - 1 - r + 0.5
        ax.text(cx, cy + 0.13, f"{v:.1f}",
                ha="center", va="center", fontsize=5.5,
                color=tc, fontfamily="monospace")
        ax.text(cx, cy - 0.18, f"({fa_grid[r, c]})",
                ha="center", va="center", fontsize=4.5,
                color=tc, fontfamily="monospace")

ax.set_xlim(0, ncols)
ax.set_ylim(0, nrows)
ax.set_aspect("equal")
ax.axis("off")

ax.set_title(
    f"Step {STEP}  |  Day {days:.1f}  |  Cycle avg burnup: {cycle_bu:.3f} GWd/T\n"
    f"Assembly burnup — avg: {valid_values.mean():.2f}  "
    f"max: {valid_values.max():.2f}  min: {valid_values.min():.2f}  GWd/T",
    fontsize=9, pad=6
)

# ── Custom Discrete Legend ────────────────────────────────────────────────────
legend_elements = [
    Patch(facecolor=group_data[0]['color'], edgecolor='grey', label=f"Q1 (Avg: {group_data[0]['avg']:.1f} GWd/T)"),
    Patch(facecolor=group_data[1]['color'], edgecolor='grey', label=f"Q2 (Avg: {group_data[1]['avg']:.1f} GWd/T)"),
    Patch(facecolor=group_data[2]['color'], edgecolor='grey', label=f"Q3 (Avg: {group_data[2]['avg']:.1f} GWd/T)"),
    Patch(facecolor=group_data[3]['color'], edgecolor='grey', label=f"Q4 (Avg: {group_data[3]['avg']:.1f} GWd/T)")
]
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5), title="Burn-up Quartiles")

plt.tight_layout()

if SAVE:
    import os
    os.makedirs(os.path.dirname(SAVE), exist_ok=True)
    plt.savefig(SAVE, dpi=300)
    print(f"Plot saved to {SAVE}")
else:
    plt.show()