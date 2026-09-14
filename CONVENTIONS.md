# bulkiRNA — internal conventions

Repo-internal. `.Rbuildignore`d, not shipped. Read this before writing a line of
`R/`; it is the contract the parallel module agents (B1–B5) share.

Live state and the plan of record are `bulkiRNA-memory/plans/HANDOFF.md` in the sibling memory
repo ([docs/MEMORY.md](docs/MEMORY.md)), which is the entry point to every other plan document. The original API-design record is outside this repo in
`/data1/users/antonz/pipeline/sciagent-rna/docs/` (`07_api-design.md`); it is history now, not
the contract.

---

## 1. Ownership

| File | Owner |
|---|---|
| `DESCRIPTION` | Step A + integrator |
| `NAMESPACE` | **roxygen only** — never hand-edited, never committed by hand |
| `R/utils.R` | Step A + integrator |
| `R/gs-result.R` | Step A + integrator |
| `R/gs-matrix.R` | Step A + integrator |
| `R/bulkiRNA-package.R` | Step A + integrator |
| `tests/golden/`, `tests/fixtures/` | Step 0 — **do not touch** |
| everything else in `R/` and `tests/testthat/` | the module agent that created it |

Need a new dependency, a new shared helper, or a change to `gs_result`? **Put it
in your handback report.** Do not edit the file. Five agents editing one file is
the failure this table exists to prevent.

## 2. No new shared helpers

There is exactly one place shared helpers live: `R/utils.R`, and it currently
holds three things — `ensure_dir()`, `ensure_parent_dir()`, `%||%`. If you find
yourself writing a fourth, request it in the handback. This is the rule that
stops the `ensure_dir` collision from happening a second time.

Module-private helpers stay in the module file, prefixed `.` and documented
`@keywords internal`.

## 3. Roxygen

- Markdown roxygen (`Roxygen: list(markdown = TRUE)`). Use `[fn()]`, backticks,
  and `#'` throughout — no `\code{}`.
- **Every** function gets `@param` for each formal and an `@return`. Exported
  functions get a runnable `@examples`; internal ones do not need them.
- Non-exported functions: `@keywords internal` and **no `@export`**. They are
  still callable from other files in the package.
- Exported functions: `@export`. The export list is the ~30 names in
  `07_api-design.md` §7. Do not export anything not on it.
- S3 methods: `@export` on the method (roxygen emits `S3method()`), plus
  `@importFrom` for a generic owned by another package (e.g.
  `@importFrom tibble as_tibble`).
- Package imports go in the roxygen block of the function that needs them via
  `@importFrom pkg fn`. **`@import ggplot2` appears in exactly one place —
  `R/bulkiRNA-package.R`** — so no other file needs `ggplot2::` prefixes and no
  two files fight over it. Everything else uses `pkg::fn()` or `@importFrom`.
- `Suggests` packages are reached only behind a guard:

  ```r
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("`method = \"gsva\"` requires the GSVA package. Install it with ",
         "BiocManager::install(\"GSVA\").", call. = FALSE)
  }
  ```

## 4. Naming

- `snake_case` for functions and objects; `.` prefix for module-private helpers.
- Prefixes are the namespace: `gsdb_*` providers, `gs_*` compute/result/persist,
  `gs_plot_*` renderers, `de_*` differential expression, `theme_bulki`.
- File names mirror content: `R/gsdb-load.R`, `R/gs-test.R`, `R/gs-plot-dot.R`.
  Hyphens in file names, underscores in function names.
- Test files mirror the R file: `R/gs-test.R` → `tests/testthat/test-gs-test.R`.

## 5. Errors and messages

- `stop(..., call. = FALSE)` — the internal call stack is noise to a user.
- Say what was wrong, name the offending argument in backticks, and where useful
  say what to do: `` "`species` must be one of \"Mus musculus\", \"Homo sapiens\"; got \"mouse\"." ``
- Never `warning()` for something the user must act on, and never silently
  return an empty result where an error is the truth. An empty `gs_result` with
  zero rows is a valid answer; a `NULL` is not.
- No `print()`/`cat()` for progress in compute functions. Use `message()`, gated
  behind a `verbose` argument that defaults to `FALSE`.

## 6. Compute never plots; viz never computes

Structural, not aspirational:

- **Compute** (`gs_ranks`, `gs_test`, `gs_score`, `gs_leading_edge`) returns
  objects. No file I/O. No plotting. No `ggsave`.
- **Persist** (`gs_write`, `gs_read`, `gs_save`) is the *only* code that touches
  disk.
- **Render** (`gs_plot_*`, `de_*`) takes an object, returns a ggplot. No
  statistics computed inside a renderer beyond what the geom needs — if a plot
  needs a number the object does not have, the object's contract is wrong.
- `gs_save(plot, path)` writes the figure **and** its same-stem source table,
  PDF + PNG from one plot object. That craft rule is now the function's job, not
  an agent's memory.

## 7. The `stat_type` vocabulary

`gs_result$stat` is meaningless without `gs_result$stat_type`. The closed set,
from `gs_stat_types()`:

| `stat_type` | Method | Axis label |
|---|---|---|
| `NES` | fgsea | `NES` |
| `t` | GSVA → limma | `t statistic` |
| `log2_fold_enrichment` | ORA (`fgsea::fora`, `foldEnrichment`) | `log2 fold enrichment` |
| `signed_log10p` | reserved (camera/roast, no effect size) | `signed -log10 p` |
| `pct_var` | GESECA | `% variance explained` |

Rules:

1. Never write `stat` without `stat_type`.
2. Renderers get the axis label from `gs_stat_label(res)` — **never** a literal
   `"NES"`. A GSVA t-statistic mislabelled NES is exactly the class of quiet
   error this column exists to kill.
3. Adding a value means adding it to `gs_stat_types()` in `R/gs-result.R` — a
   Step A file, so it goes through the handback.
4. Signed statistics use `"up"` / `"down"` / `"ns"`, derived from
   `sign(stat)` by `gs_direction()`. The unsigned GESECA `pct_var` statistic
   uses `NA`, so a direction filter cannot silently call every result up. Not
   a factor, not `"UP"`, not `1`/`-1`.

`gs_matrix` scores are *not* a `stat_type`; they carry `score_type` (`"gsva"`,
`"ssgsea"`, `"zscore"`, `"plage"`), and become `stat_type = "t"` only once
`gs_test()` runs limma over them.

## 8. Style

- Two spaces, no tabs. 80-character soft limit. Native pipe `|>`, not `%>%`.
- `TRUE`/`FALSE`, never `T`/`F`. Explicit `return()` only for early returns.
- ggplot2 is **4.0.3**: use `colour = "transparent"`, never `colour = NA` — `NA`
  drops shape-21 points as missing values.
- `seq_len()` / `seq_along()`, never `1:n`.
- Vectorized `vapply(..., FUN.VALUE=)` over `sapply()`.

## 9. Testing

- `tests/testthat/`, edition 3, one file per `R/` file.
- Test the contract, not the implementation: shape, types, error messages,
  boundary cases. Plot tests assert on `ggplot_build(p)$data`, never pixels.
- `tests/testthat/helper-gs.R` holds the small builders (`fake_gs_result()`,
  `fake_gs_matrix()`). Extend it rather than re-rolling a fixture per file.
- **`tests/golden/` and `tests/fixtures/` are off limits.** They are the
  behavioural baseline for the whole refactor; changing them destroys the only
  evidence that the refactor preserved results.
- Before handback: `devtools::load_all()` and your own tests pass.

## 10. Running R

There is no R on the host. Everything runs in a throwaway container, and both
`--user` and `HOME` are mandatory (`saveRDS` permissions; msigdbr's runtime
cache):

```bash
DK() { docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 "$@"; }

DK Rscript -e 'devtools::document()'
DK Rscript -e 'devtools::load_all(); testthat::test_local()'
DK Rscript tests/golden/verify_golden.R      # the Step C gate — must exit 0
```

R 4.5.3 · ggplot2 4.0.3 · fgsea 1.36.2 · msigdbr 26.1.0 (downloads MSigDB at
runtime to `$HOME/.cache/R/msigdbr`; needs network on a cold cache) · GSVA 2.4.9
(param-object API: `gsvaParam()` then `gsva(param)`).

## 11. Reference data

Processed databases ship in `inst/extdata/<dir>/processed/<Species>/*.rds` and
resolve **only** via `system.file("extdata", ..., package = "bulkiRNA")`. Never
walk up from a script's own location; `.resolve_toolkit_dir()` did that with
`sys.frame()$ofile` and worked only during sourcing.

Raw source files stay in `data/references/` in the git checkout and are
`.Rbuildignore`d. Anything requiring them (`rebuild = TRUE`) must fail with an
explicit "source checkout required" message rather than a missing-file error.

### Two tiers, and only two (ADR-004, in `bulkiRNA-memory/adr/`)

**If it can be re-fetched from a versioned source, it is refcache. If it cannot, it
is bundled. Size is the tiebreaker, not the criterion.**

| Tier | Where | Versioned by |
|---|---|---|
| **Bundled** — no API, no fetchable source | `inst/extdata/`, registered in `METADATA.yaml` | the package version |
| **refcache** — has a fetchable versioned source | `$REFCACHE_ROOT/<source>/current` | the snapshot tag |

Large refcache-backed data resolves through **one** internal helper, `.ref_path()`:

```r
.ref_path("coresh")   # explicit arg -> $REFCACHE_ROOT/coresh/current -> error
```

Three rules, all of them things we got wrong at least once:

1. **No new environment variable per source.** `Sys.getenv("CORESH_CHUNKS")` is exactly
   how three mechanisms became three. `REFCACHE_ROOT` is the only one.
2. **No host path in the package** — not in a default, not in an example. `/data2/...`
   is one deployment, bound in at run time.
3. **Resolving `current` must report the tag it landed on**, and callers must put that
   tag in their provenance. Freshness is the right default for a shared cache; a
   *silent* change of bytes under a re-run is not. Do not prevent the drift — make it
   impossible for the drift to be silent.

## 11a. The `gs_db` provider contract

**Frozen. B1 implements it; B2, B3 and B4 consume it.** Every `gsdb_*` provider
returns the same object, whatever its source — MSigDB, a bundled RDS, a GMT, a
user-registered list:

```r
structure(
  sets,                      # named list of character vectors: id -> gene symbols
  pathway_names = <named chr>,   # id -> human-readable label, names() == names(sets)
  database      = <chr(1)>,      # provider label, lands in gs_result$database
  species       = <chr(1)>,      # "Mus musculus" / "Homo sapiens"
  gene_id_type  = <chr(1)>,      # "symbol" (only value today)
  class = "gs_db"
)
```

A **named list of character vectors** is the contract because it is exactly what
`fgsea::fgseaMultilevel()`, `fgsea::fora()` and `GSVA::gsvaParam()` all take.
The old `list(T2G =, T2N =)` pair was a `clusterProfiler::GSEA()` input format
and dies with clusterProfiler. B1 provides an internal `.gsdb_as_t2g(db)`
returning `list(T2G = data.frame(gs_name, gene_symbol), T2N = data.frame(gs_name,
description))` **solely** so Step C's `R/deprecated.R` can keep
`load_reference_db()` returning its old shape. Nothing new calls it.

`database` is the **stable snake_case registry key** — `mitopathways`,
`transportdb`, `mito_unified`, `msigdb_H`, `msigdb_C2_CP_KEGG` — never a
human-readable label. It lands in `gs_result$database`, which is a **join and
filter key**: `gs_filter(res, database == "mitopathways")` and `rbind()` across
databases both need it machine-typeable and stable. Display strings drift; keys
must not. The pretty label rides along as a `database_label` attribute, and
renderers prettify from it — that is the render layer's job, not the data
layer's.

Multiple databases in one call are a **named list of `gs_db`** objects; the list
name becomes `gs_result$database`. Empty sets are dropped at construction, not
at test time. Set names are unique within a `gs_db`; providers prefix
(`MITOPATHWAYS_`, `TRANSPORTDB_`) as they do today.

## 12. Deletion is the job — but the freeze outranks it

The surface goes from 123 definitions to ~30 exports. If a function has no
consumer and no internal caller, **delete it** rather than porting it. Report
what you deleted.

**Consumer evidence** lives at
`/data1/users/antonz/pipeline/sciagent-rna/docs/data/used-functions.tsv` — outside
this repo. There is no `data/used-functions.tsv` here; earlier drafts said there
was.

### The rule when §12 and the freeze collide

`02_api-inventory.md` §5 froze **24 exported names** as of `752481f`. 21 left the public API
in 1.0.0 and survive as non-exported fixtures; 3 remain exported and frozen.
`07_api-design.md` §7 lists the **~30 new exports**. The two lists are not the
same set, and a name can be on the frozen list while having no place in the new
surface. When that happens:

> A name on the frozen 24 may be **renamed** — Step C's `R/deprecated.R` shims
> it. It may **not be deleted without leaving a body the shim can call.**
> §12 governs the ~102 internal functions, not the frozen 24.

The test is **not** "does anything call it" but **"can Step C still make the old
call work"**. Renaming `list_reference_dbs()` → `gsdb_list()` passes: the shim
calls the new name. Keeping `parse_gmx()` as an internal `.gsdb_parse_gmx()`
passes. Deleting a frozen export outright **fails**, because there is then
nothing for the shim to delegate to and 64 call sites in 10 projects break with
"function not found".

Consequence for the deprecation-shim contract: when you rename or internalise a
frozen export, its **old signature must remain reproducible** — same formals, in
order, with the same defaults. You may change the new function freely; the shim
absorbs the difference. Say in your handback which frozen names your module
covers and how the shim reaches them.

If you believe a frozen export genuinely should die, that is a **scope decision
for the user, not a refactor decision.** Report it; do not act on it.

### Documented exceptions

- `msigdbr` prints a once-per-session notice when it maps human sets to mouse
  orthologs. It is correctness-relevant and not suppressible through its API, so
  it stays visible despite §5's "no progress output unless `verbose`". Do not
  write a test asserting that MSigDB providers are silent.

---

## 13. Commit messages and prose

### Commit messages

A commit message says what changed and why, for someone reading this repo later.

Keep deliberation out of it. Reasoning about what we might do next, which layer could
eventually be removed, or how confident we are belongs in the memory repo, not in history.

**Never write a linkable issue reference to another project's repository.** GitHub reads
`owner/repo#123` and full issue URLs in commit messages and cross-posts the whole commit body
onto that issue. A message written for this repo then lands on a stranger's tracker. Write
`issue 123` in plain text instead. This already happened once, on the msigdbr tracker (issue 62).

Our own issues are fine to reference.

### Prose

State what holds, rather than what does not. `install_failures.csv` "keeps its meaning:
requested and failed", not "is deliberately **not** written to". The build "proceeds whatever it
finds", not "never fails".

Short sentences. One idea per line. First person where a person did something.

Things to cut, all of which read as machine-written:

| Cut | Example |
|---|---|
| Aphoristic closers that restate the moral | "A workaround that silently stops working restores exactly the failure it was meant to prevent." |
| Punchy two-word fragments | "Nothing fails." |
| Headings that frame an argument instead of naming a section | "### Why this is easy to miss" |
| `X rather than Y` where the plain statement suffices | "depends on call order rather than call arguments" → "depends on call order" |
| Menus of alternatives nobody asked for | "If that is unattractive: include the queried genes in the cache key, or drop the memoisation." |
| Em-dashes used for drama | use a comma, a full stop, or a colon |

Sparing use of any of these is fine. Three in a paragraph is the tell.
