"""
Regenerate the paper's appendix figures/tables from outputs.csv.

Reads final/outputs.csv (one row per (dataset, constant_multiplier) run)
and writes:
  fig_appendix_memory.png       -- peak RSS (MB) vs. n, at the (dataset,
                                    const) row used for the main-body
                                    minor-size result (paper Fig. 4's
                                    underlying run)
  fig_appendix_runtime_full.png -- end-to-end time_elapsed_seconds (search
                                    + final run, at the const that gave the
                                    largest h*) vs. n, contrasted with the
                                    main body's final_run_seconds-only view
                                    (paper Fig. 5)
  fig_appendix_ell_sensitivity.png -- largest h* found vs. constant
                                    multiplier c, one line per dataset,
                                    split into road/social panels

Usage:
    python make_appendix_figures.py [outputs.csv]
"""

import csv
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

DATA_FILE = sys.argv[1] if len(sys.argv) > 1 else "outputs.csv"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

ROAD_COLOR = "#e07b39"
SOCIAL_COLOR = "#5b8dd9"


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


def group_by_dataset(rows):
    by_dataset = defaultdict(list)
    for r in rows:
        by_dataset[r["dataset"]].append(r)
    return by_dataset


def best_minor_row(group):
    """Same selection rule as make_figures.py: the run with the largest h*."""
    return max(group, key=lambda r: int(r["largest_h_minor_model"]))


def scatter_fig(road, social, ylabel, title, out_name, yscale="log"):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(*zip(*road), color=ROAD_COLOR, marker="o", s=55, zorder=3, label="Road networks")
    ax.scatter(*zip(*social), color=SOCIAL_COLOR, marker="o", s=55, zorder=3, label="Social networks")

    ax.set_xscale("log")
    if yscale:
        ax.set_yscale(yscale)

    ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1, 3], numticks=20))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_n))

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


def memory_fig(by_dataset):
    road, social = [], []
    for dataset, group in by_dataset.items():
        row = best_minor_row(group)
        if not row["peak_rss_mb"]:
            continue
        pt = (int(row["n"]), float(row["peak_rss_mb"]))
        (social if is_social(dataset) else road).append(pt)
    road.sort(); social.sort()
    scatter_fig(
        road, social, "Peak RSS (MB)", "Peak memory usage",
        "fig_appendix_memory.png", yscale="log",
    )


def runtime_full_fig(by_dataset):
    road, social = [], []
    for dataset, group in by_dataset.items():
        row = best_minor_row(group)
        if not row["time_elapsed_seconds"]:
            continue
        pt = (int(row["n"]), float(row["time_elapsed_seconds"]))
        (social if is_social(dataset) else road).append(pt)
    road.sort(); social.sort()
    scatter_fig(
        road, social, "Time (seconds)",
        "End-to-end time (search + final run) at the largest $h$ found",
        "fig_appendix_runtime_full.png", yscale="log",
    )


def _short_label(dataset):
    name = dataset.rsplit("/", 1)[-1]
    for suffix in (" Road Network's Edges.csv", "_edges.csv", ".csv", ".txt"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _plot_sensitivity_panel(ax, by_dataset, social, title):
    series = []
    for dataset, group in by_dataset.items():
        if is_social(dataset) != social:
            continue
        pts = sorted(
            (int(r["constant_multiplier"]), int(r["largest_h_minor_model"]))
            for r in group
        )
        if pts:
            series.append((dataset, pts))
    # Order legend/colors by final (largest-c) h value so nearby lines get
    # visually distinct colors instead of clustering by insertion order.
    series.sort(key=lambda item: item[1][-1][1])

    # Okabe-Ito: a color-vision-deficiency-safe categorical palette. With
    # 15+ series, no palette makes every line individually distinguishable
    # by color alone, so we also cycle marker shape and line style every
    # len(OKABE_ITO) lines -- two series that land on the same color get a
    # different marker/dash combination instead.
    OKABE_ITO = [
        "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#000000",
    ]
    MARKERS = ["o", "s", "^", "D", "v", "P"]
    LINESTYLES = ["-", "--", "-.", ":"]

    for i, (dataset, pts) in enumerate(series):
        xs, ys = zip(*pts)
        color = OKABE_ITO[i % len(OKABE_ITO)]
        marker = MARKERS[(i // len(OKABE_ITO)) % len(MARKERS)]
        linestyle = LINESTYLES[(i // len(OKABE_ITO)) % len(LINESTYLES)]
        ax.plot(xs, ys, marker=marker, markersize=4, linewidth=1.3,
                 linestyle=linestyle, label=_short_label(dataset),
                 color=color, alpha=0.9)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Constant multiplier $c$", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(x))))
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=7, frameon=True, loc="center left",
              bbox_to_anchor=(1.02, 0.5), borderaxespad=0)


def ell_sensitivity_fig(by_dataset):
    """One line per dataset: h* found vs. constant multiplier c. Road and
    social networks are rendered as separate figures (rather than side by
    side) so each gets enough width for its legend and for individual
    lines to stay visually distinguishable on a log y-axis."""
    for social, name, title in [
        (False, "road", "Road networks"),
        (True, "social", "Social networks"),
    ]:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        _plot_sensitivity_panel(ax, by_dataset, social, title)
        fig.suptitle(
            r"Sensitivity of largest clique minor found to $c$"
            "\n"
            r"(where $\ell = c\sqrt{n}/(h\sqrt{\ln n})$)",
            fontsize=12, fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        out = os.path.join(OUT_DIR, f"fig_appendix_ell_sensitivity_{name}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out}")
        plt.close(fig)


def main():
    rows = load_rows()
    by_dataset = group_by_dataset(rows)
    memory_fig(by_dataset)
    runtime_full_fig(by_dataset)
    ell_sensitivity_fig(by_dataset)


if __name__ == "__main__":
    main()
