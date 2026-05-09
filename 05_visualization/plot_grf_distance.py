"""
plot_grf_distance.py

Purpose:    Plot pre-computed GRF/MDS coordinates (from grf_dis_plot_clean.R). Color = model, shape = matrix completeness, fill style = data type.
Usage:      python plot_grf_distance.py --coords_in grf_mds_coords.tsv --map all_trees.map --out grf_mds.png
Paper ref:  Suppl Fig 2 (GRF panels)
Software:   Python 3.10+, numpy, pandas, matplotlib
"""

import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------- map reader with ID bridge ----------
def read_map(path: str, tree_ids: list) -> pd.DataFrame:
    """
    Load 3-column map: TreeID, Name, Path.
    Parses Name -> (Model label, Matrix int, DataType).
    DataType options: 'Good', 'Bad', 'Raw'
    """
    meta = pd.read_csv(path, sep=r"\s+|\t+", engine="python",
                       names=["TreeID", "Name", "Path"])
    rf_ids = list(map(str, tree_ids))
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
                raise KeyError("Tree IDs in coordinates and map do not match.")
    else:
        meta = meta.set_index("TreeID").loc[rf_ids].reset_index()

    # Parse logic
    def parse(name: str):
        base = str(name).strip()
        parts = base.rsplit("_", 1)

        if len(parts) != 2 or not parts[1].isdigit():
            return ("unknown", "Good", None)

        full_code = parts[0].lower()
        matrix = int(parts[1])

        # 1. Determine Data Type
        if full_code.startswith("bad_"):
            data_type = "Bad"
            core_code = full_code[4:]
        elif full_code.startswith("raw_"):
            data_type = "Raw"
            core_code = full_code[4:]
        elif full_code.startswith("unroot_"):
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
    ap = argparse.ArgumentParser(description="Plot pre-calculated MDS coordinates.")
    ap.add_argument("--coords_in", required=True, help="Path to TSV with TreeID, MDS1, MDS2")
    ap.add_argument("--map", help="3-col map (TreeID, Name, Path)")
    ap.add_argument("--out", default="grf_mds.png", help="Output file")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--jitter", type=float, default=0.0)
    args = ap.parse_args()

    # Data Loading
    try:
        coords_df = pd.read_csv(args.coords_in, sep=r"\t|,|\s+", engine="python")
    except Exception as e:
        raise ValueError(f"Failed to read coordinates file: {e}")

    if not {"TreeID", "MDS1", "MDS2"}.issubset(coords_df.columns):
        raise ValueError("The input coordinates file must contain 'TreeID', 'MDS1', and 'MDS2' columns.")
    tree_ids = coords_df["TreeID"].tolist()

    # Metadata
    if args.map:
        meta = read_map(args.map, tree_ids)
    else:
        meta = pd.DataFrame({"TreeID": tree_ids, "Model": "unknown",
                             "DataType": "Good", "Matrix": 50})

    # Merge coords and metadata
    plotdf = pd.merge(coords_df, meta, on="TreeID")

    # Optional Jitter
    if args.jitter > 0:
        rng = np.random.default_rng(args.seed)
        plotdf["MDS1"] += rng.normal(0, args.jitter, len(plotdf))
        plotdf["MDS2"] += rng.normal(0, args.jitter, len(plotdf))

    # --- COLOR PALETTE ---
    color_map = {
        "PMSF":       "#D62728",
        "wASTRAL":    "#8C564B",
        "ASTRAL":     "#FF7F00",
        "LG":         "#FFD700",
        "GHOST":      "#0000FF",
        "Dayhoff6":   "#008000",
        "ASTRAL_IV":  "#800080",
        "unknown":    "#D3D3D3"
    }

    # --- SHAPE PALETTE ---
    shape_map = {
        50: "o", 60: "v", 70: "D", 80: "s", 90: "^", 100: "P"
    }

    # ---- Plotting ----
    fig, ax = plt.subplots(figsize=(10.0, 7.0))
    plt.subplots_adjust(left=0.08, right=0.70, top=0.92, bottom=0.10)

    ax.axhline(0, color="#e0e0e0", lw=1, zorder=0)
    ax.axvline(0, color="#e0e0e0", lw=1, zorder=0)

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
            kwargs["facecolors"] = "none"
            kwargs["edgecolors"] = c_val
            kwargs["linewidths"] = 1.5
        elif dtype == "Raw":
            kwargs["facecolors"] = c_val
            kwargs["edgecolors"] = "black"
            kwargs["linewidths"] = 0.5
            kwargs["hatch"] = "////"
            kwargs["alpha"] = 0.7
        else:
            kwargs["facecolors"] = c_val
            kwargs["edgecolors"] = "black"
            kwargs["linewidths"] = 0.5
        ax.scatter(sub["MDS1"], sub["MDS2"], **kwargs)

    ax.set_xlabel("MDS Dimension 1")
    ax.set_ylabel("MDS Dimension 2")
    ax.set_title("Generalized Robinson-Foulds Distances", fontsize=14)

    # ---- LEGENDS ----
    used_models = sorted(plotdf["Model"].unique())
    model_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map.get(m, "gray"),
               markeredgecolor='k', markersize=10, label=m)
        for m in used_models
    ]
    leg_model = ax.legend(handles=model_handles, title="Model",
                          loc="upper left", bbox_to_anchor=(1.05, 1.0))
    ax.add_artist(leg_model)

    used_mats = sorted([m for m in plotdf["Matrix"].unique() if pd.notna(m)])
    mat_handles = [
        Line2D([0], [0], marker=shape_map.get(m, "o"), color='w',
               markerfacecolor='gray', markeredgecolor='k', markersize=10, label=f"{m}%")
        for m in used_mats
    ]
    leg_mat = ax.legend(handles=mat_handles, title="Matrix Completeness",
                        loc="upper left", bbox_to_anchor=(1.05, 0.60))
    ax.add_artist(leg_mat)

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
    raw_handle = Line2D([0], [0], label="All (Stripe)")

    ax.legend(handles=[good_handle, raw_handle, bad_handle],
              handler_map={raw_handle: HatchedCircleHandler()},
              title="Data Type",
              loc="upper left", bbox_to_anchor=(1.05, 0.30))

    fig.savefig(args.out, dpi=args.dpi)
    print(f"[OK] Figure saved to {args.out}")

if __name__ == "__main__":
    main()