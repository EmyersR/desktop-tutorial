# 01_compare_bed.R
# First analysis: compare MCF7 (cancer) vs MCF10A (normal) BED files

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("rtracklayer", "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
library(rtracklayer)
library(GenomicRanges)

# File paths
bed_mcf7   <- "data/MCF7.final.transcripts.bed"
bed_mcf10a <- "data/MCF10A.final.transcripts.bed"

# Import as GRanges
mcf7   <- import(bed_mcf7, format = "BED")
mcf10a <- import(bed_mcf10a, format = "BED")

# Compare overlaps
hits <- findOverlaps(mcf7, mcf10a, ignore.strand = TRUE)
mcf7_unique   <- mcf7[-unique(queryHits(hits))]
mcf10a_unique <- mcf10a[-unique(subjectHits(hits))]

# Output summary
cat("MCF7 total:", length(mcf7), "\n")
cat("MCF10A total:", length(mcf10a), "\n")
cat("MCF7 unique:", length(mcf7_unique), "\n")
cat("MCF10A unique:", length(mcf10a_unique), "\n")

# Save unique transcripts
dir.create("results", showWarnings = FALSE)
export(mcf7_unique, "results/MCF7_unique.bed", format = "BED")
export(mcf10a_unique, "results/MCF10A_unique.bed", format = "BED")