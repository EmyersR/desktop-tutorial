# 02_chr_summary.R
# Summarize transcripts per chromosome

# Load required library
library(GenomicRanges)

# File paths
bed_mcf7   <- "data/MCF7.final.transcripts.bed"
bed_mcf10a <- "data/MCF10A.final.transcripts.bed"

# Import as GRanges
mcf7   <- import(bed_mcf7, format = "BED")
mcf10a <- import(bed_mcf10a, format = "BED")

# Count transcripts per chromosome
mcf7_counts   <- table(seqnames(mcf7))
mcf10a_counts <- table(seqnames(mcf10a))

# Combine into summary table
chr_summary <- data.frame(
  Chromosome = union(names(mcf7_counts), names(mcf10a_counts)),
  MCF7_Counts = as.numeric(mcf7_counts[union(names(mcf7_counts), names(mcf10a_counts))]),
  MCF10A_Counts = as.numeric(mcf10a_counts[union(names(mcf7_counts), names(mcf10a_counts))])
)

# Replace NAs with 0
chr_summary[is.na(chr_summary)] <- 0

# Save output
dir.create("results", showWarnings = FALSE)
write.csv(chr_summary, "results/chr_summary.csv", row.names = FALSE)

# Print preview
print(head(chr_summary))