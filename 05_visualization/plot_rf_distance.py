"""
plot_rf_distance.py

Purpose:    MDS plot from an IQ-TREE *.rfdist matrix. Color = model, shape = matrix completeness, fill style = data type.
Usage:      python plot_rf_distance.py --rfdist rf_all.rfdist --map all_trees.map --out rf_mds.png --coords rf_mds_coords.tsv
Paper ref:  Figure 3, Suppl Fig 2 (standard RF panel)
Software:   Python 3.10+, numpy, pandas, scikit-learn, matplotlib
"""

import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------- robust *.rfdist reader ----------
def read_rfdist(path: str) -> pd.DataFrame:
    """
    Load an IQ-TREE *.rfdist file to a square numeric DataFrame.
    """
    # Detect and skip a leading "n n" size header
    skiprows = 0
    with open(path, "r") as fh:
        first = fh.readline().strip()
    if re.match(r"^\s*\d+\s+\d+\s*$", first):
        skiprows = 1

    # Read remaining lines as whitespace-delimited
    raw = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, skiprows=skiprows)
    raw = raw.dropna(how="all")
    if raw.empty:
        raise ValueError("RF file appears empty after skipping header.")

    r, c = raw.shape

    def _is_number(s):
        try:
            float(s)
            return True
        except Exception:
            return False

    # Case A: header row present -> (n+1) x (n+1)
    if c == r and r > 1 and not _is_number(str(raw.iloc[0, 1])):
        header = raw.iloc[0, 1:].astype(str).tolist()
        data = raw.iloc[1:, :]
        data.index = data.iloc[:, 0].astype(str)
        df = data.iloc[:, 1:].copy()
        df.columns = [str(x) for x in header]

    # Case B: no explicit header -> n x (n+1), first col = row labels
    elif c == r + 1 and r > 0:
        idx = raw.iloc[:, 0].astype(str)
        df = raw.iloc[:, 1:].copy()
        df.index = idx
        df.columns = [str(x) for x in idx.values]  # symmetric matrix

    else:
        # Fallback: try header with index in col 0
        df = pd.read_csv(path, sep=r"\s+", header=0, index_col=0, skiprows=skiprows)

    # numeric & align
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[df.columns.astype(str), df.columns.astype(str)]
    if df.shape[0] != df.shape[1]:
        raise ValueError("RF matrix is not square after parsing/alignment.")
    return df

# ---------- map reader with ID bridge ----------
def read_map(path: str, rfdist_index: pd.Index) -> pd.DataFrame:
    """
    Load 3-column map: TreeID, Name, Path.
    Parses Name -> (Model label, Matrix int, DataType).
    DataType options: 'Good', 'Bad', 'Raw'
    """
    meta = pd.read_csv(path, sep=r"\s+|\t+", engine="python",
                       names=["TreeID", "Name", "Path"])
    rf_ids = list(map(str, rfdist_index))
    map_ids = set(meta["TreeID"].astype(str))

    # Bridge ID schemes if necessary
    if not set(rf_ids).issubset(map_ids):
        bridged = []
        ok = True
        for rid in rf_ids:
            m = re.match(r"^Tree(\d+)$", rid)
            if not m:
                ok = False
                break
            bridged.append(f"Tree_{int(m.group(1)) + 1:03d}")
        if ok and set(bridged).issubset(map_ids):
            meta = meta.set_index("TreeID").loc[bridged].reset_index()
            meta["TreeID"] = rf_ids
        else:
            meta2 = meta.copy()
            def to_rf(x):
                m = re.match(r"^Tree_(\d+)$", str(x))
                return f"Tree{int(m.group(1)) - 1}" if m else str(x)
            meta2["RF_ID"] = meta2["TreeID"].apply(to_rf)
            if set(rf_ids).issubset(set(meta2["RF_ID"])):
                meta = (meta2.set_index("RF_ID").loc[rf_ids].reset_index().rename(columns={"RF_ID": "TreeID"}))
            else:
                raise KeyError("Tree IDs in rfdist and map do not match.")
    else:
        meta = meta.set_index("TreeID").loc[rf_ids].reset_index()

    # Parse logic
    def parse(name: str):
        # Format: [prefix_]code_matrix
        base = str(name).strip()
        parts = base.rsplit("_", 1)



        if len(parts) != 2 or not parts[1].isdigit():
            return ("unknown", "Good", None)



        full_code = parts[0].lower()
        matrix = int(parts[1])



        # 1. Determine Data Type (Good, Bad, Raw)
        if full_code.startswith("bad_"):
            data_type = "Bad"
            core_code = full_code[4:] # strip 'bad_'
        elif full_code.startswith("raw_"):
            data_type = "Raw"
            core_code = full_code[4:] # strip 'raw_'
        elif full_code.startswith("unroot_"):
             # Handle unroot prefix if it slipped in, recursively
             return parse(base[7:])
        else:
            data_type = "Good"
            core_code = full_code

        # 2. Determine Model Label
        if core_code == "lg":            label = "LG"
        elif core_code == "ghost":       label = "GHOST"
        elif core_code == "pmsf":        label = "PMSF"
        elif core_code == "dayhoff6":    label = "Dayhoff6"
        elif core_code == "astral":      label = "ASTRAL"
        elif core_code == "astral_old":  label = "ASTRAL"
        elif core_code == "wastral":     label = "wASTRAL"
        elif core_code == "astral4":     label = "ASTRAL_IV"
        else:                            label = core_code.upper()

        return (label, data_type, matrix)

    parsed = meta["Name"].apply(parse)
    meta["Model"] = [x[0] for x in parsed]
    meta["DataType"] = [x[1] for x in parsed]
    meta["Matrix"] = [x[2] for x in parsed]
    return meta

def main():
    ap = argparse.ArgumentParser(description="MDS plot from IQ-TREE RF distance matrix.")
    ap.add_argument("--rfdist", required=True, help="Path to *.rfdist file")
    ap.add_argument("--map", help="3-col map (TreeID, Name, Path)")
    ap.add_argument("--ntaxa", type=int, default=None, help="Normalize RF by 2*(n-3)")
    ap.add_argument("--out", default="rf_mds.png", help="Output file")
    ap.add_argument("--coords", default=None, help="Save coordinates to TSV")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--jitter", type=float, default=0.0)
    args = ap.parse_args()

    # Data Loading
    D = read_rfdist(args.rfdist)
    if args.ntaxa:
        D = D.astype(float) / float(2 * (args.ntaxa - 3))

    # MDS
    mds = MDS(n_components=2, dissimilarity="precomputed",
              random_state=args.seed, n_init=8, max_iter=3000)
    coords = mds.fit_transform(D.values)
    coords_df = pd.DataFrame(coords, index=D.index, columns=["MDS1", "MDS2"])

    # Metadata
    if args.map:
        meta = read_map(args.map, D.index)
    else:
        meta = pd.DataFrame({"TreeID": D.index, "Model": "unknown",
                             "DataType": "Good", "Matrix": 50})
    plotdf = meta.copy()
    plotdf["MDS1"] = coords_df.loc[plotdf["TreeID"]].values[:, 0]
    plotdf["MDS2"] = coords_df.loc[plotdf["TreeID"]].values[:, 1]



    # Optional Jitter
    if args.jitter > 0:
        rng = np.random.default_rng(args.seed)
        plotdf["MDS1"] += rng.normal(0, args.jitter, len(plotdf))
        plotdf["MDS2"] += rng.normal(0, args.jitter, len(plotdf))

    # --- COLOR PALETTE (Requested) ---
    color_map = {
        "PMSF":       "#D62728",  # Red
        "wASTRAL":    "#8C564B",  # Brown
        "ASTRAL":     "#FF7F00",  # Orange (Old Astral)
        "LG":         "#FFD700",  # Yellow/Gold
        "GHOST":      "#0000FF",  # Blue
        "Dayhoff6":   "#008000",  # Green
        "ASTRAL_IV":  "#800080",  # Purple
        "unknown":    "#D3D3D3"
    }



    # --- SHAPE PALETTE ---
    shape_map = {
        50: "o", 60: "v", 70: "D", 80: "s", 90: "^", 100: "P"
    }

    # ---- Plotting ----
    # 10 inches wide, 7 inches high
    fig, ax = plt.subplots(figsize=(10.0, 7.0))



    # --- MANUAL LAYOUT CONTROL ---
    # Reserve the left 70% for the plot, leaving 30% empty on the right.
    # right=0.70 means plot goes from left=0.08 to x=0.70.
    plt.subplots_adjust(left=0.08, right=0.70, top=0.92, bottom=0.10)

    ax.axhline(0, color="#e0e0e0", lw=1, zorder=0)
    ax.axvline(0, color="#e0e0e0", lw=1, zorder=0)

    # We iterate by group to apply styles efficiently
    # DataType Logic:
    # Good = Filled solid
    # Bad  = Hollow (face="none", edge=color)
    # Raw  = Hatched (face=color, hatch="////")



    for (model, matrix, dtype), sub in plotdf.groupby(["Model", "Matrix", "DataType"]):
        marker = shape_map.get(matrix, "o")
        c_val  = color_map.get(model, "gray")



        kwargs = {
            "s": 120,
            "marker": marker,
            "alpha": 0.9,
            "label": f"{model}_{matrix}"
        }

        if dtype == "Bad":
            # Hollow
            kwargs["facecolors"] = "none"
            kwargs["edgecolors"] = c_val
            kwargs["linewidths"] = 1.5
        elif dtype == "Raw":
            # Hatched / Patterned
            kwargs["facecolors"] = c_val
            kwargs["edgecolors"] = "black"
            kwargs["linewidths"] = 0.5
            kwargs["hatch"] = "////"
            kwargs["alpha"] = 0.7
        else:
            # Good (Standard)
            kwargs["facecolors"] = c_val
            kwargs["edgecolors"] = "black"
            kwargs["linewidths"] = 0.5

        ax.scatter(sub["MDS1"], sub["MDS2"], **kwargs)

    ax.set_xlabel("MDS Dimension 1")
    ax.set_ylabel("MDS Dimension 2")
    ax.set_title("Robinson-Foulds Distances", fontsize=14)

    # ---- LEGENDS (Stacked Vertically in the empty space) ----
    # bbox_to_anchor 1.05 places legends safely in the right-hand empty space



    # 1. Model Legend
    used_models = sorted(plotdf["Model"].unique())
    model_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map.get(m, "gray"),
               markeredgecolor='k', markersize=10, label=m)
        for m in used_models
    ]
    leg_model = ax.legend(handles=model_handles, title="Model",
                          loc="upper left", bbox_to_anchor=(1.05, 1.0))
    ax.add_artist(leg_model)

    # 2. Matrix Completeness Legend
    used_mats = sorted([m for m in plotdf["Matrix"].unique() if pd.notna(m)])
    mat_handles = [
        Line2D([0], [0], marker=shape_map.get(m, "o"), color='w',
               markerfacecolor='gray', markeredgecolor='k', markersize=10, label=f"{m}%")
        for m in used_mats
    ]
    leg_mat = ax.legend(handles=mat_handles, title="Matrix Completeness",
                        loc="upper left", bbox_to_anchor=(1.05, 0.60))
    ax.add_artist(leg_mat)

    # 3. Data Type Legend (Circle Patch with Hatch)
    class HatchedCircleHandler(object):
        def legend_artist(self, legend, orig_handle, fontsize, handlebox):
            x0, y0 = handlebox.xdescent, handlebox.ydescent
            width, height = handlebox.width, handlebox.height
            patch = plt.Circle([x0 + width/2, y0 + height/2], radius=6,
                               facecolor='gray', edgecolor='k', hatch='////',
                               lw=0.5, transform=handlebox.get_transform())
            handlebox.add_artist(patch)
            return patch

    good_handle = Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                         markeredgecolor='k', markersize=10, label="Retained (Solid)")



    bad_handle = Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
                        markeredgecolor='gray', markeredgewidth=1.5, markersize=10, label="Excluded (Hollow)")

    # "Raw" handle is now mapped to "All (Stripe)" label
    raw_handle = Line2D([0], [0], label="All (Stripe)")

    ax.legend(handles=[good_handle, raw_handle, bad_handle],
              handler_map={raw_handle: HatchedCircleHandler()},
              title="Data Type",
              loc="upper left", bbox_to_anchor=(1.05, 0.30))

    # Save
    fig.savefig(args.out, dpi=args.dpi)
    print(f"[OK] Figure saved to {args.out}")



    if args.coords:
        out = plotdf[["TreeID", "Name", "Model", "Matrix", "DataType", "MDS1", "MDS2"]]
        out.to_csv(args.coords, sep="\t", index=False)
        print(f"[OK] Coords saved to {args.coords}")

if __name__ == "__main__":

    main()