library(cicero)
library(monocle3)
library(Matrix)
library(irlba)
library(SingleCellExperiment)

## ── Parse arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript step1_preprocess.R <atac_csv> <out_prefix> <chrom_sizes>")
}

input_file   <- args[1]
output_file  <- args[2]
chrom_sizes  <- args[3]

cat("Input file    :", input_file,  "\n")
cat("Output prefix :", output_file, "\n")
cat("Chrom sizes   :", chrom_sizes, "\n")

## ── Load ATAC data ────────────────────────────────────────────
indata <- read.csv(input_file, header = TRUE, row.names = 1, sep = ",")
cat("Dimensions:", dim(indata), "\n")

## ── Build sparse peak x cell matrix ───────────────────────────
sparse_matrix <- as(as.matrix(indata), "dgCMatrix")

## ── Cell metadata ─────────────────────────────────────────────
cellinfo <- data.frame(
  cells = colnames(sparse_matrix),
  row.names = colnames(sparse_matrix)
)

## ── Peak metadata ─────────────────────────────────────────────
peak_names <- rownames(sparse_matrix)

peakinfo <- data.frame(site_name = peak_names)

peak_coords <- do.call(rbind, strsplit(peak_names, "[:-]"))
peakinfo <- cbind(peak_coords, peakinfo)

names(peakinfo)[1:3] <- c("chr", "bp1", "bp2")

peakinfo$bp1 <- as.numeric(peakinfo$bp1)
peakinfo$bp2 <- as.numeric(peakinfo$bp2)

# monocle3 requires gene_short_name in gene_metadata
peakinfo$gene_short_name <- peakinfo$site_name

rownames(peakinfo) <- peakinfo$site_name

## Sanity checks ────────────────────────────────────────────────
stopifnot(identical(rownames(cellinfo), colnames(sparse_matrix)))
stopifnot(identical(rownames(peakinfo), rownames(sparse_matrix)))

## ── Build monocle3 CDS ────────────────────────────────────────
input_cds <- monocle3::new_cell_data_set(
  expression_data = sparse_matrix,
  cell_metadata   = cellinfo,
  gene_metadata   = peakinfo
)

## ── QC filtering ──────────────────────────────────────────────
counts_mat <- SingleCellExperiment::counts(input_cds)

# Remove peaks with zero counts
keep_peaks <- Matrix::rowSums(counts_mat) != 0
input_cds <- input_cds[keep_peaks, ]

# Recalculate after peak filtering
counts_mat <- SingleCellExperiment::counts(input_cds)

max_count <- 50000
min_count <- 2000

cell_counts <- Matrix::colSums(counts_mat)
keep_cells <- cell_counts >= min_count & cell_counts <= max_count

input_cds <- input_cds[, keep_cells]

cat("After QC dimensions:", dim(input_cds), "\n")

## ── Monocle3 preprocessing ────────────────────────────────────
set.seed(2017)

# monocle3 equivalent preprocessing
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- monocle3::estimate_size_factors(input_cds)

## Cicero requires a reduced coordinate matrix whose rownames are cell IDs.
counts_mat <- SingleCellExperiment::counts(input_cds)
cell_by_peak <- Matrix::t(counts_mat)

svd_out <- irlba::irlba(cell_by_peak, nv = 2)

reduced_coords <- svd_out$u %*% diag(svd_out$d)
rownames(reduced_coords) <- rownames(cell_by_peak)
colnames(reduced_coords) <- c("Dim1", "Dim2")

## ── Build Cicero CDS ──────────────────────────────────────────
cicero_cds <- cicero::make_cicero_cds(
  input_cds,
  reduced_coordinates = reduced_coords
)

## Save Cicero CDS object
cicero_cds_file <- paste0(output_file, "_cicero_cds.Rds")
saveRDS(cicero_cds, cicero_cds_file)
cat("Saved cicero CDS to:", cicero_cds_file, "\n")

## ── Load chromosome lengths ───────────────────────────────────
chromosome_length <- read.table(chrom_sizes)
cat("Chromosome length file loaded:", chrom_sizes, "\n")

## ── Run Cicero ────────────────────────────────────────────────
conns <- cicero::run_cicero(cicero_cds, chromosome_length)

## Save connections
connfile <- paste0(output_file, "_cicero_connections.Rds")
saveRDS(conns, connfile)
cat("Saved connections to:", connfile, "\n")

## ── Save peak list and connections as CSV ─────────────────────
peakfile <- paste0(output_file, "_all_peaks.csv")

all_peaks <- rownames(input_cds)
all_peaks_cleaned <- gsub("[:\\-]", "_", all_peaks)

write.csv(
  x = data.frame(peak = all_peaks_cleaned),
  file = peakfile,
  row.names = FALSE
)

cat("Saved all peaks to:", peakfile, "\n")

cicero_conn <- paste0(output_file, "_cicero_connections.csv")
write.csv(x = conns, file = cicero_conn, row.names = FALSE)
cat("Saved cicero connections CSV to:", cicero_conn, "\n")

cat("=== step1_preprocess.R complete ===\n")