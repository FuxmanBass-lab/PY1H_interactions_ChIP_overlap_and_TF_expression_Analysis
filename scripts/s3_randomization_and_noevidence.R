# Runs a degree-preserving randomization test on the pY1H vs ChIP-seq
# network, and builds the list of PDIs with no ChIP-seq evidence.
#
# Inputs (from s2)
#   outputs/py1h_chipseq_evidence.xlsx
#   outputs/py1h_macs2_peaks_full.xlsx
#
# Outputs
#   outputs/py1h_randomization_results.csv
#   outputs/py1h_randomization_histogram.png
#   outputs/py1h_noevidence_TF_available.xlsx

library(dplyr)
library(readxl)
library(writexl)
library(igraph)
library(reshape2)
library(ggplot2)
library(readr)

# Finds outputs/ relative to this script's own location (scripts/../outputs),
# so it works no matter what directory you run Rscript from.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(".")
}
outputs_dir <- file.path(dirname(get_script_dir()), "outputs")

evidence_path <- file.path(outputs_dir, "py1h_chipseq_evidence.xlsx")
peaks_path <- file.path(outputs_dir, "py1h_macs2_peaks_full.xlsx")

if (!file.exists(evidence_path) || !file.exists(peaks_path)) {
  stop("Run s1_generate_all_tfs_chipseq.py and s2_extract_chipseq_evidence.py first - missing:\n  ",
       evidence_path, "\n  ", peaks_path)
}

## Real (reported) py1h PDIs, restricted to TFs with GTRD data - this is the network we test
pdis <- read_excel(evidence_path, sheet = "ChIP-seq results")
pdis$name <- paste0(pdis$tf, "-", pdis$cytokine)

## Full TF x cytokine ChIP-seq evidence network (fixed reference / "true" edges)
peaks <- read_excel(peaks_path, sheet = "peaks")
reference_edges <- peaks %>% dplyr::select(tf, cytokine) %>% unique()
reference_edges$name <- paste0(reference_edges$tf, "-", reference_edges$cytokine)

cat("py1h PDIs (TF has GTRD data):", nrow(pdis), "\n")
cat("Reference ChIP-seq network edges (full TF x cytokine cross-product):",
    nrow(reference_edges), "\n")

## ---- Build networks -------------------------------------------------

test_network <- graph_from_data_frame(pdis[, c("tf", "cytokine")])
reference_network <- graph_from_data_frame(reference_edges[, c("tf", "cytokine")])

true <- as_adjacency_matrix(reference_network, sparse = FALSE)
true <- melt(true)
true <- true[true$value == 1, ]
true <- paste0(true$Var1, true$Var2)

real_info <- as_adjacency_matrix(test_network, sparse = FALSE)
real_info <- melt(real_info)
real_info <- real_info[real_info$value == 1, ]
real_info <- paste0(real_info$Var1, real_info$Var2)

real_number <- length(intersect(true, real_info))
cat("True overlap (real py1h PDIs with ChIP-seq evidence):", real_number, "\n")

## ---- Degree-preserving randomization test ----------------------------

numberofedgeswitching <- 20000
iteration <- 100000

results <- c()
for (i in 1:iteration) {
  randomized_network <- rewire(test_network, keeping_degseq(niter = numberofedgeswitching))
  randomized <- as_adjacency_matrix(randomized_network, sparse = FALSE)
  randomized <- melt(randomized)
  randomized <- randomized[randomized$value == 1, ]
  randomized <- paste0(randomized$Var1, randomized$Var2)
  results <- c(results, length(intersect(true, randomized)))
}

results_df <- data.frame(
  overlap = c(results, real_number),
  identity = c(rep("Randomized", iteration), "True")
)

mean_results <- results_df %>% filter(identity == "Randomized") %>%
  dplyr::select(overlap) %>% unlist() %>% mean()
sd_results <- results_df %>% filter(identity == "Randomized") %>%
  dplyr::select(overlap) %>% unlist() %>% sd()
zscore <- (real_number - mean_results) / sd_results
empirical_p <- mean(results >= real_number)

cat(sprintf("Randomized overlap: mean=%.2f sd=%.2f\n", mean_results, sd_results))
cat(sprintf("z-score: %.3f\n", zscore))
cat(sprintf("Empirical p (randomized >= true): %.5f (%d / %d)\n",
            empirical_p, sum(results >= real_number), iteration))

write_csv(results_df, file.path(outputs_dir, "py1h_randomization_results.csv"))

ggplot(results_df, aes(x = overlap)) +
  geom_histogram(binwidth = 1, color = "black") +
  geom_vline(xintercept = real_number, colour = "red") +
  theme_bw() +
  ggtitle("py1h PDIs vs ChIP-seq evidence - degree-preserving randomization")
ggsave(file.path(outputs_dir, "py1h_randomization_histogram.png"), width = 7, height = 5)

## ---- No-evidence set: TF has GTRD data, but this specific PDI has no peak ----

noevidence <- pdis[!(pdis$name %in% reference_edges$name), c("tf", "cytokine", "name")]
cat("No-evidence PDIs (TF has GTRD data, 0 peaks for this pair):", nrow(noevidence), "\n")

write_xlsx(noevidence, file.path(outputs_dir, "py1h_noevidence_TF_available.xlsx"))
