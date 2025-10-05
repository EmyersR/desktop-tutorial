# --- packages ---
# If needed once:
# install.packages("BiocManager")
# BiocManager::install(c("rtracklayer","GenomicRanges","TxDb.Hsapiens.UCSC.hg38.knownGene"))

library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# --- paths to your four bigWig files ---
mcf7_plus  <- "data/GROseq/MCF7/GSE96859_MCF7_GROseq_Untreated_PlusStrand.bw"
mcf7_minus <- "data/GROseq/MCF7/GSE96859_MCF7_GROseq_Untreated_MinusStrand.bw"
m10_plus   <- "data/GROseq/MCF10A/GSE96859_MCF10A_GROseq_Untreated_PlusStrand.bw"
m10_minus  <- "data/GROseq/MCF10A/GSE96859_MCF10A_GROseq_Untreated_MinusStrand.bw"

# --- import bigWigs ---
bw_import <- function(f) import(f, format = "BigWig")
mcf7_bw_plus  <- bw_import(mcf7_plus)
mcf7_bw_minus <- bw_import(mcf7_minus)
m10_bw_plus   <- bw_import(m10_plus)
m10_bw_minus  <- bw_import(m10_minus)

# --- quick QC: total signal by chromosome (writes a small CSV) ---
sum_by_chr <- function(gr) tapply(gr$score * width(gr), seqnames(gr), sum)
qc_mcf7 <- data.frame(
  chr        = names(sum_by_chr(mcf7_bw_plus)),
  mcf7_plus  = unlist(sum_by_chr(mcf7_bw_plus)),
  mcf7_minus = unlist(sum_by_chr(mcf7_bw_minus))
  
)
dir.create("results", showWarnings = FALSE)
write.csv(qc_mcf7, "results/qc_groseq_signal_by_chr_MCF7.csv", row.names = FALSE)

# --- save compact R objects for faster reload later (kept locally, not pushed) ---
saveRDS(list(plus=mcf7_bw_plus, minus=mcf7_bw_minus), "results/mcf7_bw.rds")
saveRDS(list(plus=m10_bw_plus,  minus=m10_bw_minus),  "results/mcf10a_bw.rds")

# --- build promoter mask (±2 kb from TSS) for later enhancer filtering ---
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
prom <- promoters(txdb, upstream = 2000, downstream = 500)
saveRDS(prom, "results/hg38_promoters_2kb_mask.rds")