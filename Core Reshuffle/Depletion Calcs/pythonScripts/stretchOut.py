import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Configuration ─────────────────────────────────────────────────────────────
# List .dep files in chronological order — they will be concatenated
DEP_FILES = [
    "stretchOut95.dep",
    "stretchOut90.dep",
    "stretchOut85.dep",
    "stretchOut80.dep",
]

SAVE = "../Depletion Calcs/Images/boronPlot.png"   # e.g. "summary_plot.png", or None to display interactively

# ── Parser ────────────────────────────────────────────────────────────────────
def parse_summary(filepath):
    """
    Extract rows from the 'summary:' table.
    Returns list of dicts with keys: days, burnup, boron.
    """
    with open(filepath, "r", errors="replace") as f:
        text = f.read()

    idx = text.rfind("summary:")          # use last occurrence if multiple
    if idx < 0:
        raise ValueError(f"No 'summary:' block found in {filepath}")

    rows = []
    for line in text[idx:].splitlines():
        # Data rows: start with an integer step number
        m = re.match(
            r"^\s*(\d+)\s+1\s+([\d.]+)\s+[\d.]+\s+[\d.(, )]+\s+([\d.]+)\s+[\d.]+\s+\d+\s+([\d.]+)",
            line
        )
        if m:
            rows.append({
                "days":   float(m.group(2)),
                "burnup": float(m.group(3)),
                "boron":  float(m.group(4)),
            })
    return rows


def concatenate_files(filepaths):
    """
    Parse each file and join into a single sorted-by-days list.
    Duplicate day values at file boundaries are deduplicated (keep first).
    """
    all_rows = []
    seen_days = set()
    for fp in filepaths:
        rows = parse_summary(fp)
        for row in rows:
            if row["days"] not in seen_days:
                all_rows.append(row)
                seen_days.add(row["days"])
    all_rows.sort(key=lambda r: r["days"])
    return all_rows


# ── Plot ──────────────────────────────────────────────────────────────────────
def plot(rows):
    days   = [r["days"]   for r in rows]
    burnup = [r["burnup"] for r in rows]
    boron  = [r["boron"]  for r in rows]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    fig.patch.set_facecolor("#fafaf8")
    for ax in (ax1, ax2):
        ax.set_facecolor("#fafaf8")
        ax.grid(True, linestyle="--", linewidth=0.5, color="#cccccc")
        ax.spines[["top", "right"]].set_visible(False)

    # ── Burnup ────────────────────────────────────────────────────────────────
    ax1.plot(days, burnup, color="#1D9E75", linewidth=1.8, marker="o",
             markersize=3, label="B (GWd/T)")
    ax1.set_ylabel("Burnup B (GWd/T)", fontsize=10)
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.legend(fontsize=9, framealpha=0.6)

    # ── Boron ─────────────────────────────────────────────────────────────────
    ax2.plot(days, boron, color="#D85A30", linewidth=1.8, marker="o",
             markersize=3, label="Boron (ppm)")
    ax2.set_ylabel("Boron concentration (ppm)", fontsize=10)
    ax2.set_xlabel("Days", fontsize=10)
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.legend(fontsize=9, framealpha=0.6)

    # Mark file boundaries with vertical lines
    if len(DEP_FILES) > 1:
        boundaries = set()
        seen = set()
        for fp in DEP_FILES[:-1]:          # boundary = last day of each file except last
            for row in parse_summary(fp):
                seen.add(row["days"])
            boundaries.add(max(seen))
        for bday in boundaries:
            for ax in (ax1, ax2):
                ax.axvline(bday, color="#888888", linewidth=0.8,
                           linestyle=":", label=f"File boundary ({bday:.0f} d)")

    fig.suptitle(
        f"PARCS depletion summary  —  "
        f"{len(days)} steps  |  "
        f"{days[0]:.1f} – {days[-1]:.1f} days",
        fontsize=11, y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if SAVE:
        plt.savefig(SAVE, dpi=150, bbox_inches="tight")
        print(f"Saved → {SAVE}")
    else:
        plt.show()


# ── Run ───────────────────────────────────────────────────────────────────────
rows = concatenate_files(DEP_FILES)
print(f"Loaded {len(rows)} steps from {len(DEP_FILES)} file(s)  "
      f"({rows[0]['days']:.1f} – {rows[-1]['days']:.1f} days)")
plot(rows)