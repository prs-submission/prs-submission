"""
Regenerate the paper's outputs.csv-derived figures in one pass.

Reads final/outputs.csv (one row per (dataset, constant_multiplier) run)
and writes:
  fig_clique_minor.png     -- largest K_h found per dataset             (paper Fig. 4)
  fig_runtime.png          -- running time at the largest h found       (paper Fig. 5)
  fig_relative_sep.png     -- ln(|S|)/ln(n), smallest separator found   (paper Fig. 8)

Also reads final/Final results - export.csv (one row per dataset, with a
hand-computed treewidth upper bound and densest-subgraph K_h estimate --
see BASELINE_NAME_MAP below for how its "network" column maps onto
outputs.csv's "dataset" paths) and writes:
  fig_minor_comparison.png -- PRS K_h vs. densest-subgraph estimate     (paper Fig. 7)
  fig_treewidth_ratio.png  -- PRS separator vs. treewidth bound         (paper Fig. 9)

Usage:
    python make_figures.py [outputs.csv] [export_csv]

A dataset is classified as "social" iff its path contains "/social/",
matching the DATASETS layout convention in shallow_separator.py.
"""

import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

DATA_FILE = sys.argv[1] if len(sys.argv) > 1 else "outputs.csv"
EXPORT_FILE = sys.argv[2] if len(sys.argv) > 2 else "Final results - export.csv"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

ROAD_COLOR = "#e07b39"
SOCIAL_COLOR = "#5b8dd9"

# "network" column in "Final results - export.csv" -> "dataset" column in
# outputs.csv. 
BASELINE_NAME_MAP = {
    "Oldenburg": "datasets/Oldenburg Road Network's Edges.csv",
    "San Francisco": "datasets/San Francisco Road Network's Edges.csv",
    "San Joaquin": "datasets/San Joaquin Road Network's Edges.csv",
    "North America Road": "datasets/North America Road Network's Edges.csv",
    "California Road": "datasets/California Road Network's Edges.csv",
    "Luxembourg": "datasets/luxembourg.csv",
    "Minnesota": "datasets/minnesota.txt",
    "Belgium": "datasets/belgium.csv",
    "USroads": "datasets/usroads.csv",
    "road-asia-osm": "datasets/road-asia-osm.txt",
    "road-germany-osm": "datasets/road-germany-osm.txt",
    "road-GB-osm": "datasets/road-great-britain-osm.txt",
    "road-italy-osm": "datasets/road-italy-osm.csv",
    "road-netherlands-osm": "datasets/road-netherlands-osm.csv",
    "road-CA-osm": "datasets/road-roadNet-CA.csv",
    "road-PA-osm": "datasets/road-roadNet-PA.csv",
    "road-euroroad-osm": "datasets/road-euroroad.txt",
    "facebook-combined": "datasets/social/facebook_combined.csv",
    "tvshow_edges": "datasets/social/tvshow_edges.csv",
    "public_figure_edges": "datasets/social/public_figure_edges.csv",
    "musae_RU": "datasets/social/musae_RU_edges.csv",
    "musae_PTBR": "datasets/social/musae_PTBR_edges.csv",
    "musae_git": "datasets/social/musae_git_edges.csv",
    "musae_facebook": "datasets/social/musae_facebook_edges.csv",
    "lastfm_asia": "datasets/social/lastfm_asia_edges.csv",
    "deezer_europe": "datasets/social/deezer_europe_edges.csv",
    "HR_edges": "datasets/social/HR_edges.csv",
}


def is_social(dataset_path):
    return "/social/" in dataset_path.replace("\\", "/")


def fmt_n(x, _):
    if x >= 1e6:
        return f"{x/1e6:.0f}M"
    if x >= 1e3:
        return f"{int(x/1e3)}K"
    return str(int(x))


def load_rows():
    with open(os.path.join(OUT_DIR, DATA_FILE), newline="") as f:
        return list(csv.DictReader(f))


def load_baseline():
    """Return {dataset_path: {"treewidth": int|None, "alg2_est": float}}."""
    path = os.path.join(OUT_DIR, EXPORT_FILE)
    if not os.path.exists(path):
        print(f"Note: {EXPORT_FILE} not found; skipping baseline/treewidth figures.")
        return {}

    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            network = row["network"].strip()
            dataset = BASELINE_NAME_MAP.get(network)
            if dataset is None:
                print(f"Warning: no dataset mapping for export network '{network}', skipping.")
                continue
            tw = row["treewidth"].strip()
            out[dataset] = {
                "treewidth": int(tw) if tw and tw != "n/a" else None,
                "alg2_est": float(row["alg_2.1 est"]),
            }
    return out


def per_dataset_best(rows):
    """
    Group rows by dataset. For each dataset, keep:
      - n
      - social flag
      - max largest_h_minor_model across all constant_multiplier runs,
        together with that row's final_run_seconds (paper Fig. 4 & 5:
        "the largest h found" and "running time to find it")
      - min separator_size across all constant_multiplier runs,
        together with n for computing ln(|S|)/ln(n) (paper Fig. 8)
    """
    by_dataset = defaultdict(list)
    for r in rows:
        by_dataset[r["dataset"]].append(r)

    best = {}
    for dataset, group in by_dataset.items():
        n = int(group[0]["n"])
        social = is_social(dataset)

        best_minor_row = max(group, key=lambda r: int(r["largest_h_minor_model"]))
        best_sep_row = min(
            (r for r in group if r["separator_size"]),
            key=lambda r: int(r["separator_size"]),
            default=None,
        )

        best[dataset] = {
            "n": n,
            "social": social,
            "h_max": int(best_minor_row["largest_h_minor_model"]),
            "h_max_runtime": float(best_minor_row["final_run_seconds"]),
            "sep_min": int(best_sep_row["separator_size"]) if best_sep_row else None,
        }
    return best


def split_road_social(best, key):
    road, social = [], []
    for d in best.values():
        val = d[key]
        if val is None:
            continue
        pt = (d["n"], val)
        (social if d["social"] else road).append(pt)
    return sorted(road), sorted(social)


def scatter_fig(road, social, ylabel, title, out_name, yscale="log", y_int=True):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(*zip(*road), color=ROAD_COLOR, marker="o", s=55, zorder=3, label="Road networks")
    ax.scatter(*zip(*social), color=SOCIAL_COLOR, marker="o", s=55, zorder=3, label="Social networks")

    ax.set_xscale("log")
    if yscale:
        ax.set_yscale(yscale)

    ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1, 3], numticks=20))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_n))
    if y_int and yscale == "log":
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1, 2, 5], numticks=20))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(x))))

    ax.set_xlabel("Number of vertices $n$", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=10, frameon=True)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, out_name)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out} ({len(road)} road, {len(social)} social points)")
    plt.close(fig)


def minor_comparison_fig(best, baseline):
    """PRS K_h vs. densest-subgraph estimate, per dataset (paper Fig. 7)."""
    ss_road, ss_social, alg_road, alg_social = [], [], [], []
    for dataset, d in best.items():
        b = baseline.get(dataset)
        if b is None:
            continue
        pt_ss = (d["n"], d["h_max"])
        pt_alg = (d["n"], b["alg2_est"])
        if d["social"]:
            ss_social.append(pt_ss)
            alg_social.append(pt_alg)
        else:
            ss_road.append(pt_ss)
            alg_road.append(pt_alg)

    if not (ss_road or ss_social):
        print("No matching baseline rows found; skipping fig_minor_comparison.png.")
        return

    ss_road.sort(); alg_road.sort(); ss_social.sort(); alg_social.sort()

    fig, ax = plt.subplots(figsize=(9, 5))
    for pts, color, marker, linestyle, label in [
        (ss_road, ROAD_COLOR, "o", "-", "ShallowSeparator -- road"),
        (alg_road, ROAD_COLOR, "s", "--", "Alg. 2 estimate -- road"),
        (ss_social, SOCIAL_COLOR, "o", "-", "ShallowSeparator -- social"),
        (alg_social, SOCIAL_COLOR, "s", "--", "Alg. 2 estimate -- social"),
    ]:
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.plot(xs, ys, color=color, linestyle=linestyle, marker=marker, markersize=6,
                 markerfacecolor="none" if marker == "s" else color,
                 markeredgewidth=1.6, linewidth=1.8, label=label,
                 zorder=4 if marker == "o" else 3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1, 3], numticks=20))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_n))
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1, 2, 3, 5, 7], numticks=20))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(x))))

    ax.set_xlabel("Number of vertices $n$", fontsize=11)
    ax.set_ylabel("Clique minor size $h$", fontsize=11)
    ax.set_title("Largest clique minor found vs. densest-subgraph estimate", fontsize=12, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=9, frameon=True)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig_minor_comparison.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out} ({len(ss_road)} road, {len(ss_social)} social points)")
    plt.close(fig)


def treewidth_ratio_fig(best, baseline):
    """ln(|S|)/ln(n) relative to ln(treewidth)/ln(n), per dataset (paper Fig. 9)."""
    road, social = [], []
    for dataset, d in best.items():
        b = baseline.get(dataset)
        if b is None or b["treewidth"] is None or d["sep_min"] is None:
            continue
        prs_rel = math.log(d["sep_min"]) / math.log(d["n"])
        tw_rel = math.log(b["treewidth"]) / math.log(d["n"])
        ratio = prs_rel / tw_rel
        (social if d["social"] else road).append((d["n"], ratio))

    if not (road or social):
        print("No matching treewidth rows found; skipping fig_treewidth_ratio.png.")
        return

    road.sort(); social.sort()

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(*zip(*road), color=ROAD_COLOR, marker="o", s=55, zorder=3, label="Road networks")
    if social:
        ax.scatter(*zip(*social), color=SOCIAL_COLOR, marker="o", s=55, zorder=3, label="Social networks")
    ax.axhline(1.0, color="black", linewidth=0.9, linestyle="--", zorder=2, label="PRS = treewidth bound")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=10))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_n))

    ax.set_xlabel("Number of vertices $n$", fontsize=11)
    ax.set_ylabel(r"$\ln(|S|)/\ln(n)$  /  $\ln(\mathrm{tw})/\ln(n)$", fontsize=11)
    ax.set_title("PRS separator size vs. treewidth upper bound", fontsize=12, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=10, frameon=True)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "fig_treewidth_ratio.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out} ({len(road)} road, {len(social)} social points)")
    plt.close(fig)


def main():
    rows = load_rows()
    best = per_dataset_best(rows)
    baseline = load_baseline()

    road, social = split_road_social(best, "h_max")
    scatter_fig(road, social, "$h$", "Largest clique minor $K_h$ found", "fig_clique_minor.png")

    road, social = split_road_social(best, "h_max_runtime")
    scatter_fig(road, social, "Time (seconds)", "Running time to find largest clique minor", "fig_runtime.png")

    def rel_sep(d):
        return math.log(d["sep_min"]) / math.log(d["n"]) if d["sep_min"] else None

    road, social = [], []
    for d in best.values():
        r = rel_sep(d)
        if r is None:
            continue
        (social if d["social"] else road).append((d["n"], r))
    road.sort()
    social.sort()
    scatter_fig(
        road, social, r"$\ln(|S|)\,/\,\ln(n)$", "Relative separator size",
        "fig_relative_sep.png", yscale=None,
    )

    if baseline:
        minor_comparison_fig(best, baseline)
        treewidth_ratio_fig(best, baseline)


if __name__ == "__main__":
    main()
