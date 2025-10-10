# 03_groseq_coverage.R — GRO-seq Coverage Quantification

# Load libraries
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# Load processed GRanges objects
mcf7 <- readRDS("data/processed/mcf7_transcripts.rds")
mcf10a <- readRDS("data/processed/mcf10a_transcripts.rds")

# Compute per-chromosome coverage (how many reads per region)
cov_mcf7 <- coverage(mcf7, weight = "score")
cov_mcf10a <- coverage(mcf10a, weight = "score")

# Summarize total coverage by chromosome
summary_df <- data.frame(
  Chromosome = names(cov_mcf7),
  MCF7_Coverage = sapply(cov_mcf7, sum),
  MCF10A_Coverage = sapply(cov_mcf10a, sum)
)

# Normalize (for comparison)
summary_df <- summary_df %>%
  mutate(Ratio_MCF7_to_MCF10A = MCF7_Coverage / MCF10A_Coverage)

# Save coverage summary
write.csv(summary_df, "data/processed/groseq_coverage_summary.csv", row.names = FALSE)

# Plot coverage ratio
ggplot(summary_df, aes(x = Chromosome, y = Ratio_MCF7_to_MCF10A)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "MCF7 vs MCF10A GRO-seq Coverage Ratio",
       x = "Chromosome", y = "Coverage Ratio (MCF7 / MCF10A)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/groseq_coverage_ratio.png", width = 7, height = 5)