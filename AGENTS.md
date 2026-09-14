# Repository Guidelines

> Single source of truth for agent/contributor guidance. `CLAUDE.md` and `GEMINI.md`
> intentionally contain only `@AGENTS.md` so every assistant reads the same instructions.

## What this is

**bulkiRNA is an installed R package** for bulk and pseudo-bulk RNA-seq differential
expression and gene-set analysis. It replaced `RNAseq-toolkit`, a folder of scripts vendored
as a git submodule and pulled in with `source()`.

**There is no `source()` step. Use `library(bulkiRNA)`.** If you find an instruction to source
a script out of `01_modules/RNAseq-toolkit/`, or a claim that there is no `library()` call, it
is pre-package and wrong — including in a consumer project's own documentation.

The package exists because 24 submodule copies drifted, nothing was tested, and no result
could name the code that produced it. Every rule below traces to a defect that actually
happened.

## Before you write code

1. **Read [CONVENTIONS.md](CONVENTIONS.md)** — the code contract: naming, argument validation,
   error-message style, roxygen expectations.
2. **Read the plan of record.** It lives in the memory layer, a sibling repo — see
   [docs/MEMORY.md](docs/MEMORY.md). `bulkiRNA-memory/plans/HANDOFF.md` is the entry point: what
   is in flight, what is parked, and why.
3. **Do not re-derive a convention that is already recorded.**
   `bulkiRNA-memory/plans/2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md`
   names the canonical previous implementation of each analysis, with paths, parameter choices
   and the reasoning behind them. **Read it instead of asking where prior work lives.**

## Running R

**There is no R on the host.** Everything runs in a throwaway container. Both `--user` and
`HOME` are mandatory:

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.13 <command>
```

Add `--network host` only when network access is genuinely needed. Scratch scripts go in
`/tmp`, never in the repo.

## Gates — all of them, before claiming anything works

```r
devtools::document()                    # NAMESPACE and man/ are generated
devtools::test()                        # 0 failures
```
```bash
Rscript tests/golden/verify_golden.R    # must exit 0
```

Plus `rcmdcheck::rcmdcheck()` for anything touching the package surface. It is currently
**0 errors, 0 warnings, 0 notes** and must stay there. It catches what `devtools::test()`
cannot: `tests/fixtures/` is `.Rbuildignore`d, so a fixture-backed test must *skip* in the
built package rather than error. It also catches stale `man/` and a semantic conflict between
two branches that each pass on their own — run it after every merge, not only after an edit.

**The two gates do not cover the same things, in both directions.** The enforcement tests that
parse `R/` — the RNG-ownership test and the name audit — **skip under `R CMD check`**, because an
installed package has no `.R` files. So a green check does not mean the invariants hold; only
`devtools::test()` checks them. Both skips say so out loud. `qs2` is a declared `Suggests` and is
absent from the image, so `R CMD check` needs a scratch library mounted at `R_LIBS`.

## Off-limits

- **`NAMESPACE` and `man/` are roxygen-generated.** Never hand-edit. Write `@export` /
  `@keywords internal` and run `devtools::document()`.
- **`tests/golden/` is read-only** unless deliberately re-capturing. Never run
  `capture_golden.R` without `--cases=<name>`: a bare run rewrites all 20 baselines.
- **`tests/fixtures/` is read-only.** Each fixture has a generating script beside it.
- **The 24 signature-frozen exports may not change formals** — same names, order, defaults.
  `bulkirna_api()` lists them. The freeze covers signatures, **not** behaviour: a corrected
  result is not a breaking change, and `NEWS.md` records it.
- **Consumer projects are read-only.** Mount them `:ro`. Never commit in one — several hold
  live research data with uncommitted work.

## The four architectural rules

1. **Four layers, one job each.** `gsdb_*` providers → `gs_test()`/`gs_score()` compute →
   `gs_result`/`gs_matrix` objects → `gs_plot_*` renderers. **Compute never plots. Renderers
   never compute.**
2. **Results are tibbles**, not S4. No `@result` slot; `dplyr` and `rbind()` work.
3. **The boundary rule.** A table in and a table out is a package function. A result in and a
   *decision* out stays a skill.
4. **Guards assert, they do not merely fix.** A workaround reaching into another package's
   internals must be paired with an independent check on the result. The `msigdbr` ortholog
   cache taught this: the truncation was invisible in the returned object.

## Two invariants, both with enforcement tests

- **Only `R/rng.R` may mutate RNG state.** `set.seed()`, `RNGkind(args)` and writing
  `.Random.seed` belong in `.with_pinned_seed()`. A bare `RNGkind()` is a read and is legal
  anywhere. Reason: `BiocParallel::bplapply()` switches the generator to `L'Ecuyer-CMRG` inside
  the task — even with `SerialParam()` — so seeding without pinning silently changes answers
  under parallelism. `tests/testthat/test-stochastic.R` fails if a fourth mutation site appears.
- **Every stochastic function is declared.** `bulkirna_stochastic()` names them, their seed
  argument and default, and what is random in them. Declaring a sixth fails the suite until it
  is classified as exercised there or covered elsewhere with a reason.

**The seed defaults deliberately disagree** — CoReSh `1L`, GATOM `42`, fgsea `123L`. Each
matches what already produced published numbers; `coresh_*` uses upstream's literal value.
Unifying them would silently move results. They are documented deviations, not untidiness.

## Writing tests for stochastic code

- **Never assert a p-value or a random number.** Assert relations: same seed agrees, different
  seeds differ, parent and `bplapply` worker agree, the caller's stream is untouched.
- **Derive the list under test from the code**, never from a literal, or it goes stale on the
  next addition.
- **A tiny fixture proves nothing.** Saturated p-values compare equal across seeds, so the test
  passes while asserting nothing. Use enough genes and sets to keep the estimator in the range
  where the seed matters.
- **Capture RNG state after building the fixture**, or the fixture's own `set.seed()` is
  mistaken for the function under test disturbing its caller.
- **A hand-built fixture only contains what its author imagined.** `NA` Entrez ids present in
  0.7% of real CoReSh datasets passed the entire suite. Where a real-data fixture exists, use it.
- **Watch for assertions that cannot fail**: a value compared against the function that produced
  it, a bound satisfied by every number (`is.finite(x) || is.infinite(x)`), an `expect_*` on a
  value constructed to satisfy it, or a comparison of two empty vectors.

## Reporting and honesty

- **A green run is not a passing gate.** A cached result, a skipped branch or an empty result
  can all print success. Confirm the code ran: this project has shipped an "exit 0" that never
  executed the migrated path, and a cached `NULL` that made a failure outlive its own fix.
- **An empty result and a broken run must be distinguishable.** Count errors; stop outright when
  everything failed rather than reporting a clean nothing.
- **If a brief's premise is false, stop and report it.** Do not work around it. Agents on this
  package have done that six times and were right every time.

## Commits

- Neutral or first-person, never third person. No "the owner's", no "the integrator's".
- **Never write a linkable issue reference** (`owner/repo#123`, or a full issue URL) in a commit
  message: GitHub cross-posts the whole body onto that project's tracker. Write `issue 123`.
- No internal deliberation in commit bodies. State what changed and the evidence for it.
- Do not rewrite published history without explicit consent.

## Known traps

- **ggplot2 4.0+:** use `color = "transparent"`, not `color = NA`, for shape-21 points; `NA`
  drops them as missing values.
- **`load_or_compute()`-style caches key on filename only.** Any cache whose input is another
  stage's generated artifact can serve results computed against a superseded input. A class
  change forces a cache rename.
- **`msigdbr` keys its ortholog cache on species alone**, so a second collection in one session
  is truncated to the intersection. `gsdb_msigdb()` guards this and asserts on the result.
- **Species strings matter**: `"Mus musculus"` versus `"Homo sapiens"`, and `db_species` is a
  separate argument from the target species.
