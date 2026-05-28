"""Generate TheraEpiFlow pipeline roadmap figure (PNG + SVG).

Layout: serpentine board — steps 1-5 left→right, 6-10 right→left, 11-13 left→right.
Each station: number badge, step name, short description, tool (italic), output artifact.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# ── Palette ───────────────────────────────────────────────────────────────────
CAT_COLOR = {
    "input":    "#5C8FBF",
    "predict":  "#4F81BD",
    "filter":   "#C0504D",
    "reduce":   "#E8993A",
    "annotate": "#4BACC6",
    "curate":   "#8064A2",
    "global":   "#4CAF50",
}

# ── Station data ──────────────────────────────────────────────────────────────
# (num, label, short_desc, tool, output, category)
STATIONS = [
    ("1",  "fetch_sequences",        "Download viral protein\nsequences",
     "UniProt REST API",                       "SEQUENCES.fasta",            "input"),
    ("2",  "predict_binding",        "Predict MHC-I binding\nfor all 9-mers",
     "NetMHCpan EL + MHCFlurry",              "NET_PRED / FLURRY_PRED.csv", "predict"),
    ("3",  "consensus_filter",       "Keep dual-binders\n(≤2% rank + Calis score)",
     "IEDB API + Calis 2013",                  "CONSENSUS.csv",              "filter"),
    ("4",  "screen_toxicity",        "Remove toxic peptides\n(ToxinPred3 ≥0.38)",
     "ToxinPred3",                             "TOXICITY_SAFE.csv",          "filter"),
    ("5",  "cluster_epitopes",       "Group similar sequences\n(≥80% identity)",
     "NetworkX graph",                         "CLUSTER.csv",                "reduce"),
    ("6",  "select_representatives", "Pick best peptide\nper cluster (★)",
     "Composite score ranking",                "REPRESENTATIVES.csv",        "reduce"),
    ("7",  "search_variants",        "Find natural sequence\nvariants",
     "UniProt search + pairwise align",        "VARIANTS.fasta",             "annotate"),
    ("8",  "analyze_conservation",   "Score conservation\nacross variants",
     "Window identity + BLOSUM62",             "CONSERVATION.csv",           "annotate"),
    ("9",  "population_coverage",    "Estimate HLA population\ncoverage",
     "IEDB allele-frequency DB",               "COVERAGE.csv",               "annotate"),
    ("10", "predict_murine",         "Predict binding to\nmurine MHC-I (H-2)",
     "NetMHCpan + MHCFlurry vs H-2",          "MURINE.csv",                 "annotate"),
    ("11", "curate_murine",          "Join conservation,\ncoverage & murine data",
     "Pandas JOIN",                            "CURATE_MURINE.csv",          "curate"),
    ("12", "integrate_data",         "Stack all tracks into\nmaster table",
     "Pandas concat + view",                   "MASTER_TABLE_FULL / VIEW",   "global"),
    ("13", "generate_report",        "Build interactive\nHTML report",
     "Jinja2 + embedded JS",                   "REPORT.html",                "global"),
]

N = len(STATIONS)

# ── Serpentine grid positions ─────────────────────────────────────────────────
# Row 0 (top): 1-5  left→right
# Row 1 (mid): 6-10 right→left
# Row 2 (bot): 11-13 left→right
COLS = 5
ROW_Y  = [2.0, 1.0, 0.0]
COL_X  = [0.0, 1.0, 2.0, 3.0, 4.0]

def station_pos(idx):
    """Return (cx, cy) center for station index 0-based."""
    if idx < 5:
        row, col = 0, idx
    elif idx < 10:
        row, col = 1, 4 - (idx - 5)
    else:
        row, col = 2, idx - 10
    return COL_X[col], ROW_Y[row]

# ── Figure setup ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(18, 11))
ax.set_xlim(-0.75, 5.2)
ax.set_ylim(-0.72, 2.82)
ax.axis("off")

# ── Background bands ─────────────────────────────────────────────────────────
# PER-TRACK band (steps 1-11, rows 0-2 but only col 0-2 on row 2)
ax.add_patch(FancyBboxPatch(
    (-0.60, -0.48), 5.40, 2.78,
    boxstyle="round,pad=0.02,rounding_size=0.15",
    facecolor="#F2F4F7", edgecolor="#C8CDD6", lw=1.4, zorder=0))
ax.text(-0.55, 1.0, "PER-TRACK\n(organism × protein)",
        rotation=90, va="center", ha="center", fontsize=10.5,
        fontweight="bold", color="#5B6472")

# GLOBAL band (steps 12-13, drawn separately below)
ax.add_patch(FancyBboxPatch(
    (-0.60, -0.68), 2.42, 0.53,
    boxstyle="round,pad=0.02,rounding_size=0.15",
    facecolor="#EAF6EC", edgecolor="#BBD9BE", lw=1.4, zorder=0))
ax.text(-0.55, -0.42, "GLOBAL", rotation=90, va="center", ha="center",
        fontsize=10.5, fontweight="bold", color="#3C8741")

# ── Draw connector arrows ─────────────────────────────────────────────────────
BOX_W, BOX_H = 0.82, 0.40

def draw_arrow(ax, p1, p2, color="#7A828E"):
    x1, y1 = p1
    x2, y2 = p2
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=color,
                    lw=1.8,
                    mutation_scale=14,
                    connectionstyle="arc3,rad=0.0",
                ))

for i in range(N - 1):
    x1, y1 = station_pos(i)
    x2, y2 = station_pos(i + 1)

    if y1 == y2:
        # same row — horizontal arrow from right/left edge
        if x2 > x1:  # left→right
            draw_arrow(ax, (x1 + BOX_W/2 + 0.02, y1), (x2 - BOX_W/2 - 0.02, y2))
        else:  # right→left
            draw_arrow(ax, (x1 - BOX_W/2 - 0.02, y1), (x2 + BOX_W/2 + 0.02, y2))
    else:
        # row change — step down from bottom of box on row end, into top of box on new row
        # row 0→1: step 5 (rightmost) goes down-left to step 6 (rightmost)
        # row 1→2: step 10 (leftmost col 0) goes down to step 11
        mid_x = x1
        draw_arrow(ax,
                   (x1, y1 - BOX_H/2 - 0.02),
                   (x2, y2 + BOX_H/2 + 0.02))

# "merge all tracks" arrow from step 11 down to step 12
x11, y11 = station_pos(10)
x12, y12 = station_pos(11)
ax.annotate("",
    xy=(x12, y12 + BOX_H/2 + 0.02),
    xytext=(x11, y11 - BOX_H/2 - 0.02),
    arrowprops=dict(arrowstyle="-|>", color="#3C8741", lw=2.0,
                    mutation_scale=14,
                    connectionstyle="arc3,rad=0.0"))
ax.text((x11 + x12)/2 + 0.08, (y11 + y12)/2,
        "merge all\ntracks", fontsize=8.5, color="#3C8741",
        style="italic", va="center", ha="left")

# ── ★ annotation: step 6 → feeds 8-10 ────────────────────────────────────────
x6, y6 = station_pos(5)
x8, y8 = station_pos(7)
ax.annotate(
    "★  best representative\n    feeds steps 8-10",
    xy=(x8 + BOX_W/2 + 0.02, y8 + 0.05),
    xytext=(x6 + BOX_W/2 + 0.38, y6 - 0.18),
    fontsize=8.2, color="#8064A2", ha="left",
    arrowprops=dict(arrowstyle="-", color="#8064A2", lw=0.9, ls=":"))

# ── Station boxes ─────────────────────────────────────────────────────────────
for i, (num, name, desc, tool, out, cat) in enumerate(STATIONS):
    cx, cy = station_pos(i)
    col = CAT_COLOR[cat]

    # outer box
    ax.add_patch(FancyBboxPatch(
        (cx - BOX_W/2, cy - BOX_H/2), BOX_W, BOX_H,
        boxstyle="round,pad=0.015,rounding_size=0.08",
        facecolor=col, edgecolor="white", lw=1.6, zorder=3, alpha=0.96))

    # number badge (left side)
    ax.add_patch(plt.Circle((cx - BOX_W/2 + 0.115, cy), 0.092,
                             color="white", alpha=0.30, zorder=4))
    ax.text(cx - BOX_W/2 + 0.115, cy, num,
            va="center", ha="center", fontsize=11, fontweight="bold",
            color="white", zorder=5)

    # step name
    ax.text(cx - BOX_W/2 + 0.235, cy + 0.095, name,
            va="center", ha="left", fontsize=8.6, fontweight="bold",
            color="white", zorder=5)

    # description
    ax.text(cx - BOX_W/2 + 0.235, cy - 0.07, desc,
            va="top", ha="left", fontsize=7.0, color="white",
            style="normal", zorder=5, linespacing=1.25)

    # output artifact below box
    ax.text(cx, cy - BOX_H/2 - 0.055, out,
            va="top", ha="center", fontsize=7.2, color="#3A3F47",
            style="italic", zorder=4)

    # tool text above box
    ax.text(cx, cy + BOX_H/2 + 0.025, tool,
            va="bottom", ha="center", fontsize=6.8, color="#555A62",
            style="italic", zorder=4)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_items = [
    ("input",    "Input"),
    ("predict",  "Prediction"),
    ("filter",   "Filter"),
    ("reduce",   "Reduce / Cluster"),
    ("annotate", "Annotate (qualitative)"),
    ("curate",   "Curate (join)"),
    ("global",   "Global (all tracks)"),
]
lx, ly = 4.42, 2.58
ax.text(lx, ly, "Category", fontsize=9.5, fontweight="bold",
        color="#2A2E35", va="top")
for i, (cat, label) in enumerate(legend_items):
    yy = ly - 0.24 - i * 0.26
    ax.add_patch(FancyBboxPatch((lx, yy - 0.085), 0.22, 0.17,
                 boxstyle="round,pad=0.01,rounding_size=0.04",
                 facecolor=CAT_COLOR[cat], edgecolor="none"))
    ax.text(lx + 0.30, yy, label, va="center", ha="left",
            fontsize=8.6, color="#3A3F47")

# ── Title + subtitle ─────────────────────────────────────────────────────────
ax.set_title(
    "TheraEpiFlow — MHC-I Epitope Selection Pipeline",
    fontsize=17, fontweight="bold", color="#2A2E35", pad=14)
ax.text(0.5, 1.004,
        "Per-track processing (organism × protein)  →  global aggregation  →  interactive HTML report",
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=10.0, color="#5B6472")

plt.tight_layout()

out_base = "/home/davi/TheraEPIflow/docs/figures/pipeline_workflow"
for ext in ("png", "svg"):
    fig.savefig(f"{out_base}.{ext}", dpi=300, bbox_inches="tight")
print(f"Saved {out_base}.png / .svg")
