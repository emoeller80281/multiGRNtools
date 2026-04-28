## ============================================================
##  run_directnet.R — unified Human / Mouse DIRECTNET script
##  Usage:
##    Rscript run_directnet.R <rna_file> <atac_file> <out_dir> \
##                            <sample_name> <genome> <script_dir> \
##                            <tss_file> <gtf_file>
##
##  Arguments:
##    rna_file    : full path to RNA counts CSV
##    atac_file   : full path to ATAC counts CSV
##    out_dir     : full path to output directory
##    sample_name : label used as cell identity (e.g. K562, mESC_E7.5_rep1)
##    genome      : "hg38" for Human, "mm10" for Mouse
##    script_dir  : directory containing reference files
##                  (Mouse annotation RDS)
##    tss_file    : full path to gene TSS file
##    gtf_file    : full path to gene annotation GTF/GTF.GZ file
## ============================================================

library(DIRECTNET)
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

## ── Parse arguments ───────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop(paste(
    "Usage: Rscript run_directnet.R",
    "<rna_file> <atac_file> <out_dir>",
    "<sample_name> <genome> <script_dir>",
    "<tss_file> <gtf_file> <num_cpus>"
  ))
}

rna_file    <- args[1]   # RNA counts CSV
atac_file   <- args[2]   # ATAC counts CSV
out_dir     <- args[3]   # Output directory
sample_name <- args[4]   # Cell identity label (e.g. K562)
genome      <- args[5]   # "hg38" or "mm10"
script_dir  <- args[6]   # Directory with reference files
tss_file    <- args[7]   # Gene TSS path
gtf_file    <- args[8]   # Gene annotation path
num_cpus    <- ifelse(length(args) >= 9, args[9], 8)  # Optional: number of CPUs for parallel processing

cat("======================================\n")
cat("  DIRECTNET run\n")
cat("  RNA file    :", rna_file,    "\n")
cat("  ATAC file   :", atac_file,   "\n")
cat("  Out dir     :", out_dir,     "\n")
cat("  Sample name :", sample_name, "\n")
cat("  Genome      :", genome,      "\n")
cat("  Script dir  :", script_dir,  "\n")
cat("  TSS file    :", tss_file,    "\n")
cat("  GTF file    :", gtf_file,    "\n")
cat("  Num CPUs    :", num_cpus,    "\n")
cat("======================================\n")

## ── Validate inputs ───────────────────────────────────────────
if (!file.exists(rna_file))  stop("RNA file not found: ",  rna_file)
if (!file.exists(atac_file)) stop("ATAC file not found: ", atac_file)
if (!genome %in% c("hg38", "mm10")) stop("genome must be 'hg38' or 'mm10'. Got: ", genome)
if (!file.exists(tss_file)) stop("TSS file not found: ", tss_file)
if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── Genome-specific settings ──────────────────────────────────
if (genome == "hg38") {

  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  ensdb      <- EnsDb.Hsapiens.v86
  bsgenome   <- BSgenome.Hsapiens.UCSC.hg38
  species    <- "Homo sapiens"

  annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

} else {  # mm10

  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)

  ensdb      <- EnsDb.Mmusculus.v79
  bsgenome   <- BSgenome.Mmusculus.UCSC.mm10
  species    <- "Mus musculus"

  ## Mouse annotations loaded from saved RDS (avoids EnsDb query overhead)
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
}

## ── Validate reference files ──────────────────────────────────
## Supports either DIRECTNET-style tabular files (with a genes column)
## or BED-like gene_tss files from shared reference genome directories.
read_genome_info <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  is_bed_like <- grepl("\\t", first_line) && !grepl("genes", first_line, ignore.case = TRUE)

  if (is_bed_like) {
    genome_info <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(genome_info) < 4) {
      stop("BED-like TSS file must have at least 4 columns: ", path)
    }
    colnames(genome_info)[1:4] <- c("Chrom", "Starts", "Ends", "genes")
  } else {
    genome_info <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
    if (!"genes" %in% colnames(genome_info)) {
      stop("TSS file must contain a 'genes' column: ", path)
    }

    colnames(genome_info) <- sub("^chr$", "Chrom", colnames(genome_info), ignore.case = TRUE)
    colnames(genome_info) <- sub("^chrom$", "Chrom", colnames(genome_info), ignore.case = TRUE)
    colnames(genome_info) <- sub("^starts?$", "Starts", colnames(genome_info), ignore.case = TRUE)
    colnames(genome_info) <- sub("^ends?$", "Ends", colnames(genome_info), ignore.case = TRUE)
  }

  genome_info
}

## ── Load genome info ──────────────────────────────────────────
genome.info <- read_genome_info(tss_file)
required_columns <- c("Chrom", "Starts", "Ends", "genes")
missing_columns <- setdiff(required_columns, colnames(genome.info))
if (length(missing_columns) > 0) {
  stop("TSS file is missing required columns after normalization: ", paste(missing_columns, collapse = ", "))
}

cat("Checking genome.info...\n")

# Ensure numeric
genome.info$Starts <- as.numeric(genome.info$Starts)
genome.info$Ends   <- as.numeric(genome.info$Ends)

cat("NA starts:", sum(is.na(genome.info$Starts)), "\n")
cat("NA ends  :", sum(is.na(genome.info$Ends)), "\n")

# Which genes have NA coords?
na_genes <- genome.info$genes[is.na(genome.info$Starts) | is.na(genome.info$Ends)]
head(na_genes)

genome.info <- genome.info[!is.na(genome.info$Starts) & !is.na(genome.info$Ends), ]

unik        <- !duplicated(genome.info$genes)
genome.info <- genome.info[unik, ]
cat("Genome info rows (deduplicated):", nrow(genome.info), "\n")

## ── Load counts ───────────────────────────────────────────────
cat("Loading RNA counts...\n")
rna_counts  <- read.table(rna_file,  header = TRUE, row.names = 1, sep = ",", comment.char = "")

cat("Loading ATAC counts...\n")
atac_counts <- read.table(atac_file, header = TRUE, row.names = 1, sep = ",", comment.char = "")

## ── Build Seurat object ───────────────────────────────────────
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## ── Add ATAC assay ────────────────────────────────────────────
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))

# Remove malformed peaks
valid_peaks <- 
  !is.na(start(grange.counts)) &
  !is.na(end(grange.counts)) &
  start(grange.counts) < end(grange.counts)

cat("Removing", sum(!valid_peaks), "malformed peaks\n")

grange.counts <- grange.counts[valid_peaks]
atac_counts   <- atac_counts[valid_peaks, ]

# Keep only standard chromosomes
grange.use <- as.vector(
  seqnames(grange.counts) %in% standardChromosomes(grange.counts)
)

grange.counts <- grange.counts[grange.use]
atac_counts   <- atac_counts[grange.use, ]

chrom_assay <- CreateChromatinAssay(
  counts     = atac_counts,
  sep        = c(":", "-"),
  genome     = genome,
  min.cells  = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

## ── RNA processing ────────────────────────────────────────────
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## ── ATAC processing + WNN integration ────────────────────────
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 1)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc <- FindMultiModalNeighbors(pbmc,
                                reduction.list = list("pca", "lsi"),
                                dims.list      = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

## ── Assign cell identity from sample_name arg ────────────────
## Replaces the hardcoded "Buffer1" in the original scripts
Idents(pbmc)    <- sample_name
pbmc$celltype   <- Idents(pbmc)

DefaultAssay(pbmc) <- "RNA"

##########################
## ── Seurat v5 compatibility fix ───────────────────────────────
## DIRECTNET was written for Seurat v4 and calls @counts directly on
## the RNA assay. In Seurat v5 the assay class changed to "Assay5"
## which uses a different internal structure and breaks that access.
## Downcasting back to the v4 "Assay" class fixes the error:
##   "no slot of name 'counts' for this object of class 'Assay5'"
if (inherits(pbmc[["RNA"]], "Assay5")) {
  cat("Seurat v5 detected — downcasting RNA assay to v4 Assay for DIRECTNET compatibility...\n")
  pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")
}


########################
## ── Run DIRECTNET ─────────────────────────────────────────────
markers <- row.names(pbmc)
cat("Running DIRECTNET on", length(markers), "markers...\n")

# Original markers
all_markers <- row.names(pbmc)

# Keep only genes with valid TSS
markers <- intersect(all_markers, genome.info$genes)

cat("Total genes:", length(all_markers), "\n")
cat("Markers with valid TSS:", length(markers), "\n")
cat("Genes removed:", length(setdiff(all_markers, markers)), "\n")

source(file.path(script_dir, "DIRECTNET_parallel.R"))

pbmc <- Run_DIRECT_NET_parallel(
  pbmc,
  peakcalling = FALSE,
  k_neigh = 50,
  atacbinary = TRUE,
  max_overlap = 0.5,
  size_factor_normalize = FALSE,
  genome.info = genome.info,
  focus_markers = markers,
  progress_every = 100,
  n_workers = as.integer(num_cpus),
  nthread = 1,
  coordinate_format = "directnet_original",
  parallel_backend = "mclapply"
)

direct.net_result <- Misc(pbmc, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind, direct.net_result))

direct.net_result$function_type <- gsub("HF",   "HC", direct.net_result$function_type)
direct.net_result$function_type <- gsub("Rest",  "MC", direct.net_result$function_type)
direct.net_result$function_type <- gsub("LF",   "LC", direct.net_result$function_type)

## ── Load GTF annotation ───────────────────────────────────────
cat("Loading GTF annotation from:", gtf_file, "\n")
gene_anno             <- rtracklayer::readGFF(gtf_file)
gene_anno$chromosome  <- paste0("chr", gene_anno$seqid)
gene_anno$gene        <- gene_anno$gene_id
gene_anno$transcript  <- gene_anno$transcript_id
gene_anno$symbol      <- gene_anno$gene_name

## ── Build focused markers dataframe ───────────────────────────
focused_markers <- data.frame(
  gene  = markers,
  group = rep(sample_name, length(markers))   # uses sample_name, not hardcoded
)

## ── CRE-Gene links ────────────────────────────────────────────
CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = focused_markers)

cat("CREs_Gene names:\n")
print(names(CREs_Gene))

cat("CREs_Gene distal:\n")
print(class(CREs_Gene$distal))
print(dim(as.data.frame(CREs_Gene$distal)))
print(head(as.data.frame(CREs_Gene$distal)))

cat("CREs_Gene promoter:\n")
print(class(CREs_Gene$promoter))
print(dim(as.data.frame(CREs_Gene$promoter)))
print(head(as.data.frame(CREs_Gene$promoter)))

cat("direct.net_result columns:\n")
print(colnames(direct.net_result))
print(head(direct.net_result))

## ── Variable peaks ────────────────────────────────────────────
DefaultAssay(pbmc) <- "ATAC"
variable_peaks <- VariableFeatures(pbmc)

cat("Variable peaks length:\n")
print(length(variable_peaks))
print(head(variable_peaks))

L_G_record_list <- list()
P_L_G_record_list <- list()

L_G_record_list[[sample_name]] <- as.data.frame(CREs_Gene$distal)
P_L_G_record_list[[sample_name]] <- as.data.frame(CREs_Gene$promoter)

# DIRECTNET generate_CRE expects underscore peak format
L_G_record_list[[sample_name]]$loci <- gsub("-", "_", L_G_record_list[[sample_name]]$loci)
P_L_G_record_list[[sample_name]]$loci <- gsub("-", "_", P_L_G_record_list[[sample_name]]$loci)

variable_peaks1 <- list()
variable_peaks1[[sample_name]] <- gsub("-", "_", variable_peaks)

cat("Distal overlap after format normalization:\n")
print(sum(L_G_record_list[[sample_name]]$loci %in% variable_peaks1[[sample_name]]))

cat("Promoter overlap after format normalization:\n")
print(sum(P_L_G_record_list[[sample_name]]$loci %in% variable_peaks1[[sample_name]]))

if (sum(L_G_record_list[[sample_name]]$loci %in% variable_peaks1[[sample_name]]) == 0 &&
    sum(P_L_G_record_list[[sample_name]]$loci %in% variable_peaks1[[sample_name]]) == 0) {
  stop("No CRE loci overlap variable_peaks after format normalization.")
}

## ── Focused CREs ──────────────────────────────────────────────
Focused_CREs <- generate_CRE(
  L_G_record   = L_G_record_list,
  P_L_G_record = P_L_G_record_list,
  da_peaks_list = variable_peaks1
)
cat("Focused_CREs names:\n")
print(names(Focused_CREs))

names(Focused_CREs$distal) <- sample_name
names(Focused_CREs$promoter) <- sample_name
names(Focused_CREs$L_G_record) <- sample_name
names(Focused_CREs$P_L_G_record) <- sample_name

cat("Distal CRE object:\n")
print(class(Focused_CREs$distal))
print(length(Focused_CREs$distal))
print(names(Focused_CREs$distal))
print(str(Focused_CREs$distal, max.level = 2))

cat("Promoter CRE object:\n")
print(class(Focused_CREs$promoter))
print(length(Focused_CREs$promoter))
print(names(Focused_CREs$promoter))
print(str(Focused_CREs$promoter, max.level = 2))

safe_generate_peak_TF_links <- function(peaks_bed_list, species, genome, markers, label) {
  if (is.null(peaks_bed_list) || length(peaks_bed_list) == 0) {
    warning(label, " peaks_bed_list is empty. Returning empty TF record.")
    return(data.frame())
  }

  # Drop empty list elements
  nonempty <- vapply(peaks_bed_list, function(x) {
    !is.null(x) && NROW(x) > 0
  }, logical(1))

  peaks_bed_list <- peaks_bed_list[nonempty]

  if (length(peaks_bed_list) == 0) {
    warning(label, " peaks_bed_list has no non-empty elements. Returning empty TF record.")
    return(data.frame())
  }

  # Keep markers only for groups present in CRE list, if names exist
  if (!is.null(names(peaks_bed_list)) && all(names(peaks_bed_list) != "")) {
    markers <- markers[markers$group %in% names(peaks_bed_list), , drop = FALSE]
  }

  if (nrow(markers) == 0) {
    warning(label, " markers has no groups matching peaks_bed_list. Returning empty TF record.")
    return(data.frame())
  }

  generate_peak_TF_links(
    peaks_bed_list = peaks_bed_list,
    species        = species,
    genome         = genome,
    markers        = markers
  )
}

cat("Detecting TFs for distal CREs...\n")
L_TF_record <- safe_generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$distal,
  species        = species,
  genome         = bsgenome,
  markers        = focused_markers,
  label          = "distal"
)
cat("  - Done\n")

cat("Detecting TFs for promoter CREs...\n")
P_L_TF_record <- safe_generate_peak_TF_links(
  peaks_bed_list = Focused_CREs$promoter,
  species        = species,
  genome         = bsgenome,
  markers        = focused_markers,
  label          = "promoter"
)
cat("  - Done\n")

## ── Generate network links ────────────────────────────────────
cat("Generating network links...\n")
groups <- pbmc$celltype
network_links <- generate_links_for_Cytoscape(
  L_G_record   = Focused_CREs$L_G_record,
  L_TF_record,
  P_L_G_record = Focused_CREs$P_L_G_record,
  P_L_TF_record,
  groups
)
cat("  - Done\n")

## ── Save outputs ──────────────────────────────────────────────
cat("Saving outputs...\n")
outfile <- file.path(out_dir, paste0(sample_name, "_Network_links.csv"))
write.csv(network_links, outfile, row.names = FALSE, quote = FALSE)
cat("Saved network links to:", outfile, "\n")

Node_attribute <- generate_node_for_Cytoscape(network_links, markers = focused_markers)
node_outfile   <- file.path(out_dir, paste0(sample_name, "_Node_attributes.csv"))
write.csv(Node_attribute, node_outfile, row.names = FALSE, quote = FALSE)
cat("Saved node attributes to:", node_outfile, "\n")

cat("=== DIRECTNET.R complete ===\n")
