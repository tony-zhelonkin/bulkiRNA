# `bulkiRNA` tests

Three independent layers. They fail for different reasons, which is the point --
if you change one thing and two layers go red, the third tells you which.

| Directory | What it is | Run it with |
|---|---|---|
| `testthat/` | Unit and contract tests for the package API | `devtools::test(".")` |
| `golden/` | Behavioural baselines captured from the pre-package script library | `Rscript tests/golden/verify_golden.R` |
| `fixtures/` | Shared synthetic and real-symbol inputs used by both | (not run directly) |

`golden/` and `fixtures/` are integrator-owned. Do not modify a stored baseline
to make a test pass -- see "Re-capturing a baseline" below.

## The golden gate

`golden/` holds 20 baselines recorded from the old `scripts/` library before it
was deleted, each keyed to a frozen legacy function name. `verify_golden.R`
calls those names through package internals, so it compares the old recorded
behaviour against the retained implementations in `R/legacy-fixtures-gs.R` and
`R/legacy-fixtures-plot.R`.

That distinction is the whole value of the gate. Until it was migrated, the
harness `source()`d the old scripts, so a green 20/20 proved only that the old
code still worked -- it could not see the new code at all. Pointing it at the
package found five real defects in one sitting, none of which `devtools::test()`
or `R CMD check --as-cran` could detect.

**Known blind spots.** `ggplot_build(p)$data` carries only numeric layer data,
so the baselines cannot see axis, legend, or facet label *text*, and a stored
theme records element names but not their sizes. Two of the five defects lived
exactly there. Assertions about labels and theme sizes belong in `testthat/`,
not here.

### Re-capturing a baseline

Only when a difference is understood, deliberate, and recorded. Re-capture the
specific cases and nothing else:

```bash
Rscript tests/golden/capture_golden.R --cases=gsea_barplot
```

Without `--cases=`, every baseline is overwritten and unrelated drift gets
blessed silently. `--cases=` also leaves `manifest.csv` untouched.

### Checking that the gate is still live

A green suite proves nothing if the gate cannot fail. Perturb one renderer and
confirm `verify_golden.R` reports 19 pass / 1 fail, then revert. The gate is
verified live this way rather than trusted.

## Contracts worth knowing before you touch a plot

Both were bugs once and are now asserted in `testthat/`:

- **ggplot2 4.0+ drops points with `colour = NA`** as missing values. Use
  `colour = "transparent"`. No `colour = NA` remains in `R/`; keep it that way.
- **The volcano's dashed line marks the colour boundary, not the cutoff.** In
  FDR mode it sits at `-log10(max(P.Value[adj.P.Val <= cutoff]))` -- the raw p
  of the last gene that passed -- so the line lands exactly where points stop
  being coloured. In raw-p mode it sits at `-log10(p_cutoff)`. Asserted in
  `testthat/test-de-volcano.R`.

## No R on this host

All R runs in a throwaway container. Both `--user` and `HOME` are required --
`HOME` for `msigdbr`'s runtime cache, `--user` for `saveRDS` permissions:

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 \
  Rscript tests/golden/verify_golden.R
```

The metabolic (`gatom`) tests need the other image, with `HOME=/tmp` and no
cache mount:

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD":/pkg -w /pkg scbio-singleuser:v1.7.1 \
  Rscript -e 'devtools::test(".")'
```

`devtools::test_local()` does not exist in either image. Use `devtools::test(".")`.
Scratch scripts go in `/tmp`, never in the repo.
