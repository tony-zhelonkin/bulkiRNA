# Golden-harness fixture

Synthetic bulk RNA-seq data with **planted gene-set signal**, used to capture golden
outputs before the bulkiRNA refactor and to verify them after. Fully synthetic — no real
biology, no identities.

Regenerate deterministically (`set.seed(20260807)`):

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg \
  scdock-r-dev:v0.5.11 Rscript tests/fixtures/make_fixture.R
```

Two container flags are mandatory, and both cost a cycle to discover:

- **`--user "$(id -u):$(id -g)"`** — without it the container writes as root and `saveRDS`
  fails with a permission error.
- **`-e HOME=/cache` + a mounted cache** — **msigdbr 26.1.0 no longer bundles MSigDB**; it
  downloads the release archive to `$HOME/.cache/R/msigdbr` on first use. With the default
  container `HOME=/` the path resolves to `//.cache/...` and the download fails. Network is
  required on a cold cache.

## Files (272 KB total)

| File | Contents |
|---|---|
| `counts.rds` | integer matrix, 2000 genes × 8 samples, Ensembl-style rownames |
| `metadata.rds` | `sample`, `group` (WT/KO ×4), `batch` (A/B) |
| `toptable.rds` | limma-voom `topTable`, `~ 0 + group + batch`, contrast `KO − WT`; carries `ensembl_id` + `gene_name` |
| `gene_annotation.rds` | Ensembl ID → symbol map |
| `genesets.rds` | `list(alpha = <15 sets>, beta = <8 sets>)`, keyed on symbol |
| `de_intermediates.rds` | filtered counts, voom `logcpm`, `design`, `contrast` |
| `de_real_symbols.rds` | DE frame over **real mouse symbols** (`t`, `logFC`, `P.Value`, `adj.P.Val`) |
| `ranks_real.rds` | named ranked `t` vector over the same real symbols |

### Why a real-symbol companion

`run_gsea()` and `load_reference_db()` resolve sets by **real gene symbol** (msigdbr,
MitoCarta), so synthetic `Gsym####` names return empty results and a worthless golden.
`de_real_symbols.rds` / `ranks_real.rds` cover a 4393-symbol MSigDB-mouse universe with
signal planted in two real Hallmark sets:

| Planted | Recovered by fgsea |
|---|---|
| `HALLMARK_INTERFERON_ALPHA_RESPONSE` (100) up | NES **+4.06**, padj 2e-62 |
| `HALLMARK_MYC_TARGETS_V1` (201) down | NES **−4.17**, padj 9e-128 |

`HALLMARK_INTERFERON_GAMMA_RESPONSE` also comes out strongly — expected, since it shares
genes with the alpha set. That overlap is realistic and left in.

**Reproducibility caveat:** this half depends on the MSigDB release msigdbr downloads
(2026.1 at capture) and on the human→mouse ortholog mapping msigdbr applies by default.
The synthetic half has no such dependency.

## What is planted, and why

| Set | Design intent |
|---|---|
| `ALPHA_SIGNAL_UP` (50) | Strong positive. Must come out significant and **up** |
| `ALPHA_SIGNAL_DOWN` (50) | Strong negative. Must come out significant and **down** |
| `ALPHA_NULL_SET` (50) | **Negative control** — a real set with no planted effect. Must stay non-significant |
| `ALPHA_BACKGROUND_01..12` | Random sets, so FDR has something to correct against |
| `BETA_PARTIAL_UP` | 25 genes shared with `SIGNAL_UP` + 15 random — partial overlap across databases |
| `BETA_MIXED` | 12 up + 12 down — **should cancel** and stay non-significant |
| `BETA_BACKGROUND_01..06` | Random |

`beta` deliberately shares genes with `alpha`. Independent databases are the unrealistic
case, and the overlap is what exercises cross-database plotting.

Also planted: **20 duplicated gene symbols** (unique Ensembl IDs) to exercise duplicate-ID
collapsing, and 80 scattered non-set DE genes so the volcano is not two clean blobs.

## Verified behaviour at `752481f` / `scdock-r-dev:v0.5.11`

1872 of 2000 genes survive `filterByExpr`; 27 DE at `adj.P.Val < 0.05` (16 up / 11 down).

| Engine | `SIGNAL_UP` | `SIGNAL_DOWN` | `NULL_SET` | `BETA_MIXED` |
|---|---|---|---|---|
| fgsea | NES +3.32, padj 9e-30 | NES −2.97, padj 1e-15 | padj 0.66 | padj 0.21 |
| fora | FE 7.56, padj 7e-31 | FE 5.88, padj 2e-15 | ns | — |
| GSVA → limma | t +10.9, padj 2e-10 | t −10.5, padj 3e-10 | ns | — |

All three engines agree on direction and significance while reporting different
statistics — which is what the `stat` / `stat_type` contract exists to carry.

Per-gene DE is deliberately modest (27 genes). GSEA is more sensitive than per-gene
testing, and a fixture where everything is significant would not exercise the
non-significant rendering paths.

## `coresh-chunk-micro.rds`

Four real CoReSh dataset objects, subset to 300 genes each, from mouse snapshot
`syn66227307_20260721`. 40 KB. Public GEO-derived data.

**Columns are deliberately not subset.** They are at most 20 wide, and cutting them would
destroy the `ncol(E1024) < nsamples` property that marks a PCA-reduced matrix — which is one of
the four structures this fixture exists to carry. `samples` and `nsamples` are as the source
recorded them.

| element | accession | what it is here for |
|---|---|---|
| `plain` | GSE100219 | the ordinary case |
| `duplicate_ids` | GSE10000 | repeated Entrez ids, which `fgsea::geseca()` rejects unless uniquified |
| `na_ids` | GSE101177 | 60 `NA` Entrez ids, and duplicates as well |
| `pca_reduced` | GSE100012 | a matrix whose columns are principal components, not samples |

**Why this exists.** Every CoReSh test before it built its object by hand, and a hand-built
object contains only what its author thought to put there. The `NA` Entrez ids in 0.7% of real
datasets passed the entire suite and were found by a compendium sweep instead. A fixture cut
from the real thing cannot have that blind spot.

`totalVar` is recomputed for the subset, since the stored value describes the full matrix.
`E1024` keeps its integer quantization, so the `/1024` path is exercised as it is in production.
Read it with `readRDS()` — no `qs2` needed, so these tests never skip.

Regenerate with `make_coresh_micro.R`, which needs the refcache mounted read-only and `qs2` in a
scratch library. It reproduces this file byte-for-byte.
