# Golden-output harness

Behavioural baseline for the bulkiRNA refactor, captured at `752481f` before any
restructuring. **This is Step C's gate** — the refactor is only accepted if
`verify_golden.R` exits 0 or every difference is on the recorded expected-change list.

## Why this exists

64 `source()` call sites across 10 projects depend on current behaviour. The nine legacy
suites in `tests/` assert plot *shape*, not numbers, so they would pass straight through a
silent numerical regression. This harness compares the actual computed values.

## Running

```bash
DK() { docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 "$@"; }

DK Rscript tests/golden/capture_golden.R   # (re)capture   -> data/, manifest.csv
DK Rscript tests/golden/verify_golden.R    # check         -> verify-report.csv, exit 1 on drift
```

`--user` and `HOME` are both mandatory — see `tests/fixtures/README.md`.

## What is stored

`data/<case>.rds`, one per case. Compute results are stored as returned; **plots are stored
as `ggplot_build(p)$data`** — the computed layer data, not pixels. That catches a numerical
change in what is plotted while staying stable across renderer and patch-version churn.
An S4 `gseaResult` is reduced to its `@result` frame; a patchwork is built panel-by-panel.

Comparison is `all.equal(tolerance = 1e-6, check.attributes = FALSE)`.

## Coverage — 20 cases, 0 errors

Captured: `empty_gsea_tibble`, `list_to_term2gene`, `filter_by_size`, `format_pathway_name`,
`ensure_dir`, `build_dge`, `list_reference_dbs`, `load_reference_db` (mitopathways +
transportdb), `run_gsea`, `normalize_gsea_results`, `create_standard_volcano` (×3 — fdr, p,
annotate_counts), `create_MD_plot`, `custom_minimal_theme_with_grid`, `gsea_dotplot`,
`gsea_dotplot_facet`, `gsea_barplot`, `gsea_running_sum_plot`.

Skipped, with reasons in `manifest.csv`: `download_gatom_references` (network),
`run_gsea_analysis` / `plot_all_gsea_results` (pipelines that write files and wrap
already-captured functions), `save_gsea_log` (no return value), `parse_gmx` /
`parse_mitoxplorer` (need raw files excluded from the built package),
`convert_human_to_mouse` (no fixture input).

## The harness self-tests

A gate nobody has seen fail is not a gate. Setting `GOLDEN_SELFTEST_PERTURB=1` wraps
`format_pathway_name` in `toupper()` and re-verifies. Verified behaviour:

```
clean run:      20 pass, 0 fail   -> exit 0
perturbed run:  18 pass, 2 fail   -> exit 1
    format_pathway_name     FAIL  4 string mismatches
    normalize_gsea_results  FAIL  Component "pathway_name": 50 string mismatches
```

It catches the direct change **and** its propagation into a downstream consumer. Never set
that variable in normal use.

## Expected changes during the refactor

Record differences here rather than loosening the comparison:

| Case | Expected change | Why |
|---|---|---|
| ORA (not yet captured) | results will differ | `enrichGO`/`enrichKEGG` → `fgsea::fora`; GO/KEGG now sourced from MSigDB |
| `gsea_running_sum_plot` | **complete redesign** | enrichplot dropped; B4 is a from-scratch three-panel ggplot |
| all cases | function names change | `gs_*` / `de_*` namespace; the deprecation shims keep old names callable |

`run_gsea` must **not** change: `clusterProfiler::GSEA(by = "fgsea")` and
`fgsea::fgseaMultilevel()` are the same computation. A non-empty diff there is a stop
signal, not a tolerance to widen.

## Environment at capture

R 4.5.3 · ggplot2 4.0.3 · fgsea 1.36.2 · `scdock-r-dev:v0.5.11`
