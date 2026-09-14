# Deterministic synthetic fixture for the bulkiRNA golden-output harness.
#
# Generates a small bulk RNA-seq dataset with planted gene-set signal, so golden
# outputs exercise significance, direction, and highlighting rather than a
# degenerate all-null case. No real biology or identities.
#
# Run:  docker run --rm -v <repo>:/pkg -w /pkg scdock-r-dev:v0.5.11 \
#         Rscript tests/fixtures/make_fixture.R

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

set.seed(20260807)
out_dir <- "tests/fixtures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_genes <- 2000L
n_per   <- 4L

# ---- gene identity ----------------------------------------------------------
# Symbols are synthetic. 20 symbols are duplicated on purpose so the fixture
# exercises duplicate-ID collapsing; Ensembl-style IDs stay unique.
symbols <- sprintf("Gsym%04d", seq_len(n_genes))
dup_idx <- sample(seq_len(n_genes), 20L)
symbols[dup_idx] <- symbols[dup_idx - 1L]          # each duplicates its neighbour
ensembl <- sprintf("ENSFAKE%08d", seq_len(n_genes))

# ---- design -----------------------------------------------------------------
meta <- data.frame(
  sample = sprintf("S%02d", seq_len(2L * n_per)),
  group  = factor(rep(c("WT", "KO"), each = n_per), levels = c("WT", "KO")),
  batch  = factor(rep(c("A", "B"), times = n_per)),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample

# ---- planted gene sets ------------------------------------------------------
set_up   <- sample(seq_len(n_genes), 50L)
set_down <- sample(setdiff(seq_len(n_genes), set_up), 50L)
set_null <- sample(setdiff(seq_len(n_genes), c(set_up, set_down)), 50L)

lfc <- numeric(n_genes)
lfc[set_up]   <-  rnorm(50L, mean =  1.6, sd = 0.4)
lfc[set_down] <- -rnorm(50L, mean =  1.6, sd = 0.4)
# a little unstructured signal so the volcano is not two clean blobs
scattered <- sample(setdiff(seq_len(n_genes), c(set_up, set_down, set_null)), 80L)
lfc[scattered] <- rnorm(80L, 0, 1.1)

base_mu   <- 2^runif(n_genes, 3, 11)
dispersion <- 0.15 + 1 / sqrt(base_mu)

counts <- matrix(0L, nrow = n_genes, ncol = nrow(meta),
                 dimnames = list(ensembl, meta$sample))
for (j in seq_len(nrow(meta))) {
  eff <- if (meta$group[j] == "KO") 2^lfc else rep(1, n_genes)
  bat <- if (meta$batch[j] == "B") 2^rnorm(n_genes, 0, 0.15) else rep(1, n_genes)
  counts[, j] <- rnbinom(n_genes, mu = base_mu * eff * bat, size = 1 / dispersion)
}
storage.mode(counts) <- "integer"

# ---- DE: limma-voom ---------------------------------------------------------
dge  <- edgeR::DGEList(counts = counts, samples = meta, group = meta$group)
keep <- edgeR::filterByExpr(dge, group = meta$group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- edgeR::calcNormFactors(dge)

design <- model.matrix(~ 0 + group + batch, data = meta)
colnames(design) <- make.names(colnames(design))
v    <- limma::voom(dge, design)
fit  <- limma::lmFit(v, design)
cm   <- limma::makeContrasts(KO_vs_WT = groupKO - groupWT, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))

tt <- limma::topTable(fit2, coef = "KO_vs_WT", number = Inf, sort.by = "none")
tt$ensembl_id <- rownames(tt)
tt$gene_name  <- symbols[match(rownames(tt), ensembl)]
tt <- tt[order(tt$P.Value), ]

# ---- gene-set databases -----------------------------------------------------
# Two "databases" of synthetic sets, both keyed on gene SYMBOL (the ranking key).
sym_of <- function(i) unique(symbols[i])

mk_bg <- function(prefix, n_sets, sizes) {
  stats::setNames(
    lapply(seq_len(n_sets), function(k) sym_of(sample(seq_len(n_genes), sizes[k]))),
    sprintf("%s_BACKGROUND_%02d", prefix, seq_len(n_sets))
  )
}

db_alpha <- c(
  list(
    ALPHA_SIGNAL_UP   = sym_of(set_up),
    ALPHA_SIGNAL_DOWN = sym_of(set_down),
    ALPHA_NULL_SET    = sym_of(set_null)
  ),
  mk_bg("ALPHA", 12L, sample(20:80, 12L, replace = TRUE))
)

# db_beta shares some genes with the planted sets so the two databases are not
# independent — that is the realistic case and it exercises cross-db plotting.
db_beta <- c(
  list(
    BETA_PARTIAL_UP = sym_of(c(sample(set_up, 25L), sample(seq_len(n_genes), 15L))),
    BETA_MIXED      = sym_of(c(sample(set_up, 12L), sample(set_down, 12L)))
  ),
  mk_bg("BETA", 6L, sample(15:60, 6L, replace = TRUE))
)

gene_annotation <- data.frame(
  ensembl_id = ensembl,
  gene_name  = symbols,
  stringsAsFactors = FALSE
)

# ---- write ------------------------------------------------------------------
saveRDS(counts,          file.path(out_dir, "counts.rds"))
saveRDS(meta,            file.path(out_dir, "metadata.rds"))
saveRDS(tt,              file.path(out_dir, "toptable.rds"))
saveRDS(gene_annotation, file.path(out_dir, "gene_annotation.rds"))
saveRDS(list(alpha = db_alpha, beta = db_beta),
        file.path(out_dir, "genesets.rds"))
saveRDS(list(counts_filtered = dge$counts, logcpm = v$E,
             design = design, contrast = cm),
        file.path(out_dir, "de_intermediates.rds"))

# ---- report -----------------------------------------------------------------
cat("genes (raw / kept):    ", n_genes, "/", nrow(dge), "\n")
cat("samples:               ", nrow(meta), "(", paste(levels(meta$group), collapse = "/"), ")\n")
cat("duplicated symbols:    ", sum(duplicated(symbols)), "\n")
cat("DE genes (adj.P < .05):", sum(tt$adj.P.Val < 0.05), "\n")
cat("  up / down:           ", sum(tt$adj.P.Val < .05 & tt$logFC > 0), "/",
                               sum(tt$adj.P.Val < .05 & tt$logFC < 0), "\n")
cat("db_alpha sets:         ", length(db_alpha), "\n")
cat("db_beta sets:          ", length(db_beta), "\n")

# ---- real-symbol companion --------------------------------------------------
# `run_gsea()` and `load_reference_db()` resolve gene sets by real gene symbol
# (msigdbr / MitoCarta), so synthetic symbols would return empty results. This
# builds a second ranked list over REAL mouse symbols with planted signal in two
# real Hallmark sets, so those functions have a non-degenerate golden.
#
# Reproducibility note: depends on msigdbr's bundled data (26.1.0 at capture).

set.seed(20260807)
hs <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
sets <- lapply(split(hs$gene_symbol, hs$gs_name)[c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_MYC_TARGETS_V1"
)], unique)
universe <- sort(unique(hs$gene_symbol))

t_stat <- stats::setNames(rnorm(length(universe)), universe)
t_stat[intersect(names(t_stat), sets[[1]])] <-
  rnorm(length(intersect(names(t_stat), sets[[1]])), mean =  2.6, sd = 0.8)
t_stat[intersect(names(t_stat), sets[[2]])] <-
  rnorm(length(intersect(names(t_stat), sets[[2]])), mean = -2.6, sd = 0.8)

de_real <- data.frame(
  t         = t_stat,
  logFC     = t_stat / 2.5,
  P.Value   = 2 * stats::pnorm(-abs(t_stat)),
  row.names = names(t_stat)
)
de_real$adj.P.Val <- stats::p.adjust(de_real$P.Value, "BH")
de_real <- de_real[order(de_real$P.Value), ]

saveRDS(de_real, file.path(out_dir, "de_real_symbols.rds"))
saveRDS(sort(t_stat, decreasing = TRUE), file.path(out_dir, "ranks_real.rds"))

cat("real-symbol universe:  ", length(universe), "\n")
cat("  planted UP set:      ", names(sets)[1], "(", length(sets[[1]]), ")\n")
cat("  planted DOWN set:    ", names(sets)[2], "(", length(sets[[2]]), ")\n")
