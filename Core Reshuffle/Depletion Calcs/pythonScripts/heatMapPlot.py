import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

# ── Configuration ─────────────────────────────────────────────────────────────
DEP_FILE = "resuffle1.dep"
STEP     = 26         # time step to plot (1-based)
VMIN     = 0.0        # colormap min (GWd/T)
VMAX     = 65.0       # colormap max (GWd/T)
SAVE     = './Images/test.png'       # e.g. "burnup.png", or None to display interactively

# ── Core row offsets for the 193-FA PWR layout ────────────────────────────────
# Each tuple: (column offset in the full 19-wide grid, number of assemblies)
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
def parse_step(filepath, step):
    """Return a 2-D numpy array (15 x 19) with burnup values; NaN outside core."""
    with open(filepath, "r", errors="replace") as f:
        text = f.read()

    blocks = re.findall(r"BEGIN STEP([\s\S]*?)END STEP", text)
    if step < 1 or step > len(blocks):
        raise ValueError(f"Step {step} out of range (1–{len(blocks)})")

    blk = blocks[step - 1]
    idx = blk.find("Assembly Burnup Distribution")
    if idx < 0:
        raise ValueError("'Assembly Burnup Distribution' not found in step block")

    # Collect float values until the first blank line after data starts
    values = []
    for line in blk[idx:].splitlines()[1:]:
        nums = re.findall(r"\d+\.\d+", line)
        if not nums and values:
            break
        values.extend(float(n) for n in nums)

    grid = np.full((len(ROW_LAYOUT), 19), np.nan)
    fa_grid = np.full((len(ROW_LAYOUT), 19), 0, dtype=int)
    i = 0
    fa = 1
    for r, (offset, count) in enumerate(ROW_LAYOUT):
        for c in range(count):
            grid[r, offset + c] = values[i]
            fa_grid[r, offset + c] = fa
            i += 1
            fa += 1
    return grid, fa_grid


def get_step_meta(filepath, step):
    """Return (days, cycle_burnup) for the given step."""
    with open(filepath, "r", errors="replace") as f:
        text = f.read()
    blk = re.findall(r"BEGIN STEP([\s\S]*?)END STEP", text)[step - 1]
    m = re.search(
        r"^\s*\d+\s+1\s+([\d.]+)\s+[\d.]+.*?\s+([\d.]+)\s+[\d.]+\s+\d+\s",
        blk, re.MULTILINE
    )
    return (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)


# ── Plot ──────────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    "burnup",
    [(0.00, "#E1F5EE"), (0.25, "#D0D327"), (0.55, "#F4A70F"), (1.00, "#D83030")]
)

def plot_heatmap(grid, fa_grid, step, days, cycle_bu):
    norm = mcolors.Normalize(vmin=VMIN, vmax=VMAX)
    nrows, ncols = grid.shape

    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor("#fafaf8")
    ax.set_facecolor("#fafaf8")

    for r in range(nrows):
        for c in range(ncols):
            v = grid[r, c]
            if np.isnan(v):
                continue
            color = CMAP(norm(v))
            ax.add_patch(Rectangle(
                (c, nrows - 1 - r), 1, 1,
                linewidth=0.4, edgecolor="white", facecolor=color
            ))
            lum = 0.2126*color[0] + 0.7152*color[1] + 0.0722*color[2]
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

    valid = grid[~np.isnan(grid)]
    ax.set_title(
        f"Step {step}  |  Day {days:.1f}  |  Cycle avg burnup: {cycle_bu:.3f} GWd/T\n"
        f"Assembly burnup — avg: {valid.mean():.2f}  "
        f"max: {valid.max():.2f}  min: {valid.min():.2f}  GWd/T",
        fontsize=9, pad=6
    )

    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Burnup (GWd/T)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    if SAVE:
        plt.savefig(SAVE, dpi=150, bbox_inches="tight")
        print(f"Saved → {SAVE}")
    else:
        plt.show()


# ── Run ───────────────────────────────────────────────────────────────────────
grid, fa_grid = parse_step(DEP_FILE, STEP)
days, cycle_bu = get_step_meta(DEP_FILE, STEP)
plot_heatmap(grid, fa_grid, STEP, days, cycle_bu)