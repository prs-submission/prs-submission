# ShallowSeparator

A Python implementation of the Plotkin-Rao-Smith (PRS) algorithm for finding balanced vertex separators and clique minors in real-world graphs.

This code accompanies an anonymous submission to ALENEX 2027, *Large Clique Minors Or Balanced Separators in (Road) Networks: An Experimental Study*.

---

## What it does

Given a graph *G* and an integer parameter *h*, **ShallowSeparator** outputs one of two things:

- A **balanced vertex separator** *S* — a set of vertices whose removal leaves every connected component with at most 2*n*/3 vertices, or
- A **K_h clique minor model** — *h* vertex-disjoint connected subgraphs of *G* such that every pair has an edge between them, certifying that *G* is not K_h-minor-free.

This is the *win-win framework*: the algorithm either solves the problem (separator) or explains why the standard theory predicts it cannot (large clique minor).

To find the **largest** h for which a K_h minor exists, the implementation uses an exponential search (doubling h until a separator is found) followed by a binary search over the resulting interval.

---

## Installation

```bash
pip install -r requirements.txt
```

Python 3.10 or later is required.

---

## Usage

1. Place your edge-list CSV files in a `datasets/` subdirectory.
2. Edit the `DATASETS` config table and `STRATEGY` / `CONST_VALUES` constants at the bottom of `shallow_separator.py` (see "Dataset config" below).
3. Run:

```bash
python shallow_separator.py
```

Results are appended to `outputs.csv`. Validation evidence (see "Correctness validation" below) is appended to `validation_evidence.jsonl` whenever `DEBUG_VALIDATE` is `True` (the default).

---

## Dataset config

Every dataset source used in this study has a different CSV layout: the Li et al. [LCH+05] files have a 4-column header (`Edge ID, Start Node ID, End Node ID, L2 Distance`) with **0-indexed** node IDs in columns 1,2; the Network Repository [RA15] files have a 2-column layout (no real header, just a line to skip) with **1-indexed** node IDs in columns 0,1.

```python
DATASETS = [
    {"path": "datasets/Oldenburg Road Network's Edges.csv",
     "col_u": 1, "col_v": 2, "index_base": 0},
    {"path": "datasets/belgium.csv",
     "col_u": 0, "col_v": 1, "index_base": 1},
    ...
]
```

- `col_u` / `col_v`: 0-indexed CSV columns holding the two endpoint IDs.
- `index_base`: `1` if node IDs in the file start at 1, `0` if they start at 0.

Add new datasets by adding an entry here, not by editing `find_n` / `initialize_G_and_H` directly.

---

## Input format

Each dataset must be a CSV file with one header/skip row, followed by one edge per row. Column positions and index base vary by file and are declared per file in the `DATASETS` config table — see "Dataset config" above.

---

## Output columns

| Column | Description |
|--------|-------------|
| `dataset` | Path to the input file |
| `n` | Number of vertices |
| `date` | Date the experiment was run |
| `constant_multiplier` | The constant *c* used to set ℓ = *c* · √*n* / (*h* √ln *n*) |
| `largest_h_minor_model` | Largest *h** for which a K_h minor was found |
| `smallest_h_separator` | *h** + 1; the value of *h* used to obtain the reported separator |
| `separator_size` | Number of vertices in the separator *S* |
| `max_bfs_depth` | Maximum BFS tree depth observed across all iterations |
| `iterations` | Total number of while-loop iterations |
| `first_condition_calls` | Iterations that took the shallow-tree (Case 1) branch |
| `line11_calls` | Iterations that took the deep-tree (Case 2) branch |
| `line11_median_returns` | (Strategy C only) Times the median layer itself was the chosen separator |
| `time_elapsed_seconds` | Wall-clock time for the entire `find_largest_minor_h` call at this constant multiplier. This is *not* simply `search_seconds + final_run_seconds`: `find_largest_minor_h` also makes an untimed intermediate call, between the search and the final run, that independently re-derives and validates the K_h* minor model (see `minor_model_certificate` below) — its cost is included in `time_elapsed_seconds` but not broken out into either sub-field. Typically comparable in magnitude to `final_run_seconds` |
| `search_seconds` | Wall-clock time for the exponential + binary search over *h*, i.e. everything *before* h* is known |
| `final_run_seconds` | Wall-clock time for the single `shallowSeparator` call at *h** + 1, once *h** is already known — this is what the paper's Figure 5 (running time assuming the optimal *h* is known in advance) reports |
| `peak_rss_mb` | Peak resident set size of the whole process (MB) at the time this row was written, via `resource.getrusage().ru_maxrss`. Monotonic non-decreasing across the run — not scoped to just this row — but needs no extra dependency and reflects true (not just Python-object) memory use |

`shallow_separator.py`'s own sweep loop also writes a per-dataset `outputs_sweep_summary.csv` (`sweep_total_seconds`, `peak_rss_mb`), but every value in it is exactly derivable from `outputs.csv` (sum of `time_elapsed_seconds`, max of `peak_rss_mb`, grouped by dataset) and it does not cover datasets whose `(dataset, const)` grid was split across multiple machines, so it is not included here.

---

## Regenerating the paper's figures

- **`make_figures.py [outputs.csv] [export_csv]`** — regenerates the main-body figures (largest clique minor found, running time, relative separator size, comparison to the densest-subgraph baseline, comparison to the treewidth bound) directly from `outputs.csv`. The baseline/treewidth comparison figures also read `Final results - export.csv`, a per-dataset table of treewidth upper bounds and densest-subgraph clique-minor estimates (see "Reproducing the paper's supporting baselines" below); a dataset with no matching row in that file is simply skipped for those two figures rather than erroring.
- **`make_appendix_figures.py [outputs.csv]`** — regenerates the appendix figures (peak memory usage, full end-to-end search cost, and per-dataset sensitivity of the largest clique minor found to the constant multiplier `c`).

Both scripts require `matplotlib` in addition to `sortedcontainers`. Re-running either after `outputs.csv` is updated (e.g. once a pending dataset finishes) regenerates every figure from the current data with no other steps.

---

## Layer-selection strategies

Line 11 of Algorithm 1 selects a BFS layer to add to the separator. Three strategies are implemented (set `STRATEGY` in the `__main__` block):

| Strategy | Description |
|----------|-------------|
| `'A'` | **Earliest valid** — first layer from the root satisfying Equation (3) |
| `'B'` | **Smallest valid** — fewest-vertex layer satisfying Equation (3) |
| `'C'` | **Median-preferred** (default) — search outward from the median layer; return the first valid layer found, breaking ties by size |

Strategy C is theoretically motivated (median layers yield better balance).

---

## Correctness validation

`shallow_separator.py` includes a set of validation functions for independently checking the algorithm's inputs and outputs against the original graph *G*:

| Function | Checks |
|----------|--------|
| `validate_input_graph` | Edge endpoints are in the range implied by the file's declared `col_u`/`col_v`/`index_base`, adjacency is symmetric (undirected), self-loops/duplicate edges are reported, computed *n*/*m* match expected values (if given), and the graph is connected. **Raises if the input is disconnected** — the implementation's first BFS runs over the full vertex list without first restricting to the largest connected component, so a disconnected input does not match the algorithm's assumptions. An out-of-range endpoint here means the `DATASETS` config entry for that file (column positions or index base) is wrong. For the 3 datasets that are genuinely disconnected in raw form (`usroads`, `minnesota`, `road-euroroad`), `restrict_to_largest_component` is used instead to reduce to the largest connected component before running. |
| `validate_separator` | **Two checks folded into one record and one `valid` flag**, since both bound the same object (the returned separator *S*) even though they test different properties: (1) *construction* — every vertex in *S* belongs to *G*, and removing *S* leaves every remaining component with at most 2*n*/3 vertices (checked with exact integer arithmetic: `3 * largest <= 2 * n`); (2) *size* — if the run's `l`/`h` are supplied, `\|S\|` respects the paper's Equation (1) bound (`n/l + 2(h-1)(h-2)*l*ln(n)`), evaluated at the actual `l` and `h` the run used (not the closed-form-optimal `l` from Equation 2). A run can have a well-constructed *S* that still violates the size bound, or vice versa is not possible (the bound is a guarantee whenever the algorithm returns a separator at all) — the `within_size_bound` field in the record distinguishes the size check from the (implicit) balance check if you need to tell them apart. |
| `validate_minor_model` | A **separate, unrelated** check of the returned *K<sub>h</sub>* minor model (not the separator): the branch sets number exactly *h*, are nonempty, contain only vertices of *G*, are pairwise disjoint, each induces a connected subgraph, and every pair of branch sets has at least one edge of *G* between them (all `C(h,2)` required adjacencies are checked). |
| `_validate_iteration_invariants` | Per-iteration (not just final-result) checks: H/S/K pairwise disjointness, branch-set connectivity, branch sets in K are pairwise adjacent in G (K forms an actual clique minor at every iteration, not just at the end), every branch set in K has a neighbor in H, `neighbors_in_H` consistency, and that `\|H\|` strictly decreases every iteration. Also checks, at the point each layer `X` is chosen, that `X` is nonempty and satisfies Equation (3). Gated by a separate flag, `VALIDATE_ITERATIONS` (see below) — this checks internal loop state that is never part of the reported result and cannot be independently re-derived by a reviewer, unlike the two certificate checks above. |

All validators use each vertex's `neighbors_in_G` (the original, unmutated adjacency), never `neighbors_in_H`, since the latter is destructively trimmed as the algorithm progresses and no longer reflects *G* by the time a result is returned.

### What gets checked, and when

`find_largest_minor_h` finds `h*` via exponential + binary search, inspecting only whether each intermediate call returns a minor model or a separator — it does not validate any of those intermediate results. Once `h*` is known, it makes **two extra, dedicated calls** to produce the certificates:

1. One `shallowSeparator` call at `h = h*` (guaranteed by the search to return a minor model) → validated with `validate_minor_model` → logged as a `minor_model_certificate` record.
2. One `shallowSeparator` call at `h = h* + 1` (guaranteed by Lemma 1 to return a separator) → validated with `validate_separator` (construction + size bound) → logged as a `separator_certificate` record.

These are different graph states from different calls, not the same run checked twice — and they are separate JSON records in `validation_evidence.jsonl`, not fields within a single combined record. A `(dataset, const)` row therefore normally has one `separator_certificate` and one `minor_model_certificate` (no `minor_model_certificate` when `h* == 0`, since there's no minor model to check).

The module-level flag

```python
DEBUG_VALIDATE = True   # default
```

at the top of `shallow_separator.py` controls whether:
- both certificates (separator and, when `h* >= 1`, minor model) are validated via the two dedicated extra calls described above, with the reports attached to `find_largest_minor_h`'s return values;
- each input dataset is validated once (self-loops, connectivity, n/m) before the experiment loop begins;
- validation records are appended to `validation_evidence.jsonl` — `input_graph` once per dataset, `separator_certificate` and `minor_model_certificate` once per `(dataset, const)` run. This file is the durable evidence that validation actually ran — it is not gitignored, and is intended to ship as part of the submission package.

A second flag, `VALIDATE_ITERATIONS`, independently controls the per-iteration loop-invariant checks (`_validate_iteration_invariants`) and the layer-`X` validity check. These are development/debugging aids for catching bugs before they manifest as an invalid final certificate — not additional evidence about the reported numbers, since they check internal loop state a reviewer could never independently re-derive. 

**`DEBUG_VALIDATE` defaults to `True`: the code path that produces results is, by default, the same code path that gets validated.** The paper's running-time figures were produced with `DEBUG_VALIDATE = False` — validation adds real overhead (extra calls for the certificate checks, plus per-iteration invariant checks when `VALIDATE_ITERATIONS` is also on) that is not part of the timed algorithm. Set it to `False` only to reproduce those timing numbers; leave it `True` for all other use, including regenerating separator/minor-model results.

**The `find_X_C` "no valid layer found" failure**: at very small `const` values on graphs with large `h*` (mostly social networks), the exponential search can reach an `h` so large that `ℓ` becomes tiny, and — for a small/finite graph — no BFS layer satisfies Equation (3) at all. This was confirmed as a real edge case in the paper's Case-2 layer-existence guarantee (not an implementation bug) via direct inspection of the BFS-layer state at failure: the guarantee holds asymptotically but the specific graph/parameter combination it's checked at can fall outside where it's been shown to hold. These failures are consistently confined to `const <= 3` on social-network datasets in this study's results.

`shallowSeparator` returns a `PRSResult` namedtuple with a `kind` field (`"separator"` or `"minor_model"`) rather than the previous `(SortedSet | "MINOR MODEL", ...)` mixed-type tuple, so results can be checked without relying on positional string comparisons.

## Datasets

The experiments in the paper use:

- **Road networks (Li et al.)** — 5 classical datasets (Oldenburg, San Francisco, San Joaquin, North America, California). Available from the authors of: *F. Li, D. Cheng, M. Hadjieleftheriou, G. Kollios, S. Teng. On Trip Planning Queries in Spatial Databases. SSTD 2005.*
- **Road networks (Network Repository)** — 12 larger OSM-derived datasets. Available at <https://networkrepository.com>.
- **Social networks (SNAP)** — 10 datasets. Available at <https://snap.stanford.edu/data>.

Datasets are not included in this repository.

---

## Reproducing the paper's supporting baselines

The paper compares `ShallowSeparator` against two supporting algorithms (Section 3.1–3.2 of the paper) that are described there in full but are not included in this repository as runnable code:

- **Densest-Subgraph** (Algorithm 2 in the paper) — a linear-time greedy 2-approximation for maximum subgraph density, used to derive a lower-bound estimate on clique minor size via the paper's Equation 4. Repeatedly remove a minimum-degree vertex and track the densest subgraph seen; see the paper for the exact pseudocode and the closed-form solve for `h_δ`.
- **UpperBoundTreewidth** (Algorithm 3 in the paper) — a greedy elimination-ordering heuristic that produces an upper bound on treewidth (and hence on separator size, via `t + 1`), following Maniu, Senellart, and Jog, *An Experimental Study of the Treewidth of Real-World Graph Data*, ICDT 2019 (doi:10.4230/LIPIcs.ICDT.2019.12). Our results used the **DEGREE** variant (always eliminate the minimum-degree vertex in the working graph) with the partial-degree parameter set to `9999`, i.e. large enough that a full elimination ordering is always computed rather than an early-terminated approximation.

Both algorithms are direct implementations of the pseudocode given in the paper (Algorithms 2 and 3); no third-party library was used. The scripts that produced the paper's baseline numbers are not currently in this repository — `Final results - export.csv` records their output per dataset. To reproduce them, implement Algorithms 2 and 3 exactly as specified in the paper, using the DEGREE heuristic and partial-degree cutoff above for the treewidth baseline.

---

## Citation

Citation details are withheld here to preserve double-blind review; see the paper submission for the full reference.
