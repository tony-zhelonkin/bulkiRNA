# bulkiRNA (development version)

Nothing yet. `1.1.0.9000` names the next release, not this tree.

# bulkiRNA 1.0.0

## API

* The 21 deprecated names have left the public API. Their implementations
  remain as non-exported fixtures for verification against the golden baseline
  captured at `752481f`. Code that calls these names after
  `library(bulkiRNA)` will fail to find them.

* `gs_master_columns()` returns the versioned master-table columns in emission
  order, so consumers stop reading `inst/extdata/master-schema-v1.csv`, whose
  row order is not the column order. It takes `optional`, matching the fact
  that `gs_to_master()` emits `entity_type` only when given one: the default
  describes a table built without it, and `optional = TRUE` puts it first.

* `bulkirna_api()` returns 59 rows, none deprecated. Its `frozen` column lists
  3 names rather than 24, because 21 of the frozen signatures were deprecated
  and left the public API together. `superseded_by` and `removed_in` are
  therefore `NA` throughout; both columns are kept so the registry's shape is
  stable.

* `bulkirna_stochastic()` returns 6 rows. `run_gsea()` leaves it by becoming
  internal. It still consumes randomness, so the registry now describes public
  seed interfaces only.

* Removing `convert_human_to_mouse()` from the registry also removes its
  MSigDB-only `superseded_by` claim, which did not cover custom gene sets.

## Documentation

* `README.md`, `MIGRATION.md` and `CONVENTIONS.md` no longer say the old names
  work. The 21 help pages are marked internal, so the public index stops
  documenting functions that cannot be called.

## Packaging

* `withr` is declared in `Suggests`. Two test files have used it through `::`
  since `ba910b2` without declaring it, which `R CMD check` reported as an
  unstated dependency. It is classified dev-only, so `bulkirna_check_deps()`
  does not report it as a requirement for using the package.

# bulkiRNA 0.6.0

## API

* `gs_write()` now returns its directory invisibly, as documented. Its value
  and `files` manifest attribute are unchanged.

* Validation errors consistently omit the internal call, name caller
  arguments in backticks, and render observed values with stable ASCII quotes.
  Errors caused by corrupt bundled data or malformed third-party results name
  that object instead of inventing a caller argument.

* `gsdb_coresh()` composes CoReSh search and set construction into an
  experimental provider returning `database = "coresh"`. It records the
  resolved reference snapshot, compact query names and sizes, and the
  set-building controls in database provenance. Exact Entrez IDs remain in
  the caller's query definitions instead of being serialized into the object.

* Named calls using the result-limit arguments `top` or `n_top` must now use
  `top_n`. This breaking rename affects `gs_leading_edge()`, `gs_plot_bar()`,
  `gs_plot_dot()`, `gs_plot_heatmap()`, `gs_plot_running()`,
  `coresh_loadings()`, and `coresh_sets()`. The database-provenance field
  written by `coresh_sets()` is likewise renamed from `n_top` to `top_n`.

* Named calls using the colour-palette arguments `color_palette` or `colours`
  must now use `palette`. This breaking rename affects `de_bfc_plot()`,
  `de_md_plot()`, `de_pca()`, `de_volcano()`, `gs_plot_bar()`,
  `gs_plot_dot()`, and `gs_plot_heatmap()`. No argument-level compatibility
  aliases are provided for either rename.

* `gs_coregulation()` runs public `fgsea::geseca()` on a general genes x
  samples expression matrix and returns a `gs_result` with
  `stat_type = "pct_var"`. It is experimental and is separate from
  `gs_test()` because it takes a matrix and has no contrast.

* `gs_stat_types()` gains `pct_var = "% variance explained"`, the plot-axis
  label for GESECA's unsigned statistic.

* `gs_coregulation()` records `direction = NA` on every row because percentage
  of variance explained has no up/down sign. Consumers must rank or filter it
  by `stat`, `p_value`, or `padj`; filtering its `direction` for `"up"`,
  `"down"`, or `"ns"` returns no rows. `gs_top(by_direction = TRUE)` now keeps
  an allowed missing direction as one group instead of silently dropping it.

* `gs_to_master()` deliberately rejects a `pct_var` result by default because
  its `nes` column is reserved for NES. `stat_as_nes = TRUE` remains the
  explicit override and still returns the fixed 14-column ADR-002 schema.

* `gatom_download_refs(dir = )` is the layer-prefixed GATOM reference
  downloader. The frozen `download_gatom_references(dest_dir = )` name remains
  as a deprecated compatibility shim until 1.0.0.

* Human and mouse species aliases now resolve through one internal contract.
  All species-taking functions accept `human`, `Homo sapiens`, `hsa`, `mouse`,
  `Mus musculus`, and `mmu`, including separator and case variants.
  `annotate_genes()` retains its historical partial scientific-name matching,
  while user-supplied custom `gs_db` species labels remain supported.

* A recognised alias is now normalised in the `gs_db` `species` attribute:
  `gsdb_register(species = "mouse")` records `"Mus musculus"` where it
  previously recorded `"mouse"`. An unrecognised label is still recorded as
  given, with underscores replaced by spaces.

* `annotate_genes(species = NULL)` now errors instead of silently selecting
  `"Mus musculus"`, which was the first choice returned by its historical
  `match.arg()` call. This catches unset configuration values before they can
  annotate human genes through the mouse annotation package.

## Reproducibility

* A release version now identifies exactly one commit, enforced by a test
  rather than a release checklist. Between releases the version is
  `A.B.C.9000`, where `A.B.C` names the *next* release, so a development tree
  cannot claim a released number. `0.5.0` previously described both the
  `v0.5.0` tag and 50 later commits carrying 15 additional exports.

* `gsdb_coresh()` pins the compendium search with CoReSh's historical
  `seed = 1L`. Its returned database records the snapshot tag even though the
  shared refcache remains free to advance `current`.

* `gs_coregulation()` pins GESECA's stochastic work through the package RNG
  owner and is declared by `bulkirna_stochastic()`. Its default seed is `123L`,
  matching the general fgsea adapter; CoReSh retains its historical `1L`.

## CoReSh

* CoReSh dataset identity is now the `(gse, gpl)` pair. `coresh_loadings()`
  gains optional `gpl =`; when it is omitted for a multi-platform accession,
  the function warns with all available platforms and uses the radix-first GPL
  deterministically. Its result now carries `gse` and the GPL actually used.

* `coresh_sets()` threads platform identity from search hits into loading
  extraction and set provenance. Historical hit tables without a `gpl` column
  remain usable: ambiguous accessions warn and select the radix-first platform,
  independently of chunk-file or object ordering. If two platforms for one
  query/accession are both requested, their set names receive GPL suffixes so
  neither dataset is discarded as a duplicate.

* `coresh_convergence()` now requires and groups by both `gse` and `gpl`, and
  returns both columns. `coresh_search()` already scored and retained every
  platform as a separate row; `coresh_validate()` performs no accession lookup
  or aggregation and needed no behavioural change.

* `coresh_loadings()` projects a CoReSh expression object onto a normalized
  query profile and returns the strongest absolute gene loadings.

* `coresh_sets()` turns ranked CoReSh hits into a `gs_db`, filters and
  deterministically de-duplicates the derived sets, distinguishes a valid
  empty result from a completely failed extraction, and records both database-
  and set-level provenance.

* Jaccard de-duplication in `coresh_sets()` competes across queries in global
  CoReSh rank order; the reference implementation competed in caller order,
  which accounts for the one-missing/one-extra difference from the consumer's
  stored GMT.

## Gene-set databases

* `gs_db()` accepts generic database provenance and set-keyed provenance.
  Both survive subsetting, while set-level rows are restricted and reordered
  to the retained sets. `gsdb_info()` surfaces both records for a `gs_db`.

# bulkiRNA 0.5.0

The package version is the unit of reproducibility (ADR-001), so the changes
below get a version of their own: they moved numbers in a consumer project.

## API

* `bulkirna_api()` adds a machine-readable registry of all public functions,
  keeping lifecycle status and the historical signature freeze as independent
  axes.

* `bulkirna_stochastic()` declares the five public functions that consume
  randomness, how each accepts its seed, its versioned default, and the
  operation that is stochastic. `bulkirna_api()` marks the same functions in a
  new `stochastic` column derived from that registry.

## Reproducibility

* `write_session_provenance()` now records all three `RNGkind()` components and
  this package version's stochastic defaults. The existing defaults remain
  unchanged: CoReSh retains upstream's literal `1L`, GATOM retains `42`, the
  `gs_test()` fgsea adapter retains `123L`, and signature-frozen `run_gsea()`
  retains `123`.

## Gene-set testing

* `gs_test(method = "fgsea")` now retains fgsea's `log2err` p-value uncertainty
  estimate in `gs_result`. Infinite values are preserved because they identify
  estimates past reliable resolution. Methods without a Monte-Carlo component
  omit the column, and `gs_to_master()` continues to return the fixed ADR-002
  schema without it.

## CoReSh

* `pvalues = TRUE` now works for CoReSh matching and search. Every CoReSh result
  gains a `log2err` column; this shifts positional column indices for consumers
  that do not select columns by name. Duplicate query Entrez IDs are now
  de-duplicated once with a message so `size`, `pct_var`, and `p_value` describe
  the same gene set.

## GATOM

* `gatom_module()` gains a final argument, `gene2reaction_extra`. It reaches
  `gatom::makeMetabolicGraph(gene2reaction.extra = )` and is what the combined
  network needs to see all of its edges. On the 14839 combined graph it takes
  the edge count from 7,608 to 8,179 against an unchanged 6,495 vertices. A
  value that is neither `NULL` nor a data frame with `gene` and `reaction`
  columns is rejected with a message naming the fix.

* `gatom_module()` records `attr(m, "solution_weight")`, the weight
  `mwcsr::solve_mwcsp()` reports for the module it returned. The three-call
  sequence discarded it, so consumers reaching for the solver's own objective
  had nothing to read.

* `gatom_refs()` accepts `network = "combined"` alongside `"kegg"`, resolves
  `network.combined.rds` and `met.combined.db.rds`, passes the flavour through
  to `download_gatom_references()`, and names the network in both the argument
  check and the not-found message.

## Formatting

* `format_pathway_name()` capitalizes a function word when it opens a label:
  `VIA_NEGATIVE_REGULATION` reads `Via Negative Regulation` rather than
  `via Negative Regulation`. Mid-label function words stay lowercase, so
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` still reads
  `TNF-alpha Signaling via NF-kappaB`. Only the nine function words are
  affected; Greek letters, chemical prefixes, roman numerals and lowercase gene
  symbols such as `mTOR` are untouched.

  This is a behaviour change inside a frozen signature. The export freeze
  covers formals, not output.

## Reproducibility note for existing GATOM results

`gatom::scoreGraph()` is stochastic. `gatom_module()` has seeded it since it was
written, but hand-rolled three-call sequences generally do not, and results
produced that way are one unrepeatable draw — the same input scored 165.85 and
then 169.32 in a single session. Modules recomputed through `gatom_module()`
will differ from such caches. The new values are the reproducible ones.

# bulkiRNA 0.4.0 and earlier

These releases predate this file. See the git tags and `MIGRATION.md` for the
legacy-name shims and the 24 frozen exports.
