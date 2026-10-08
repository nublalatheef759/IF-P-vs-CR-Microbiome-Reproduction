# Create phyloseq object from QIIME2 exports
# Run this FIRST before other R analysis scripts

library(phyloseq)
library(tidyverse)

cat("=== Creating phyloseq object from QIIME2 exports ===\n\n")

# Set paths to QIIME2 exported files
base_path <- "qiime2_results/exported-for-R"

# Read ASV table
asv_table <- read.delim(file.path(base_path, "asv-table.tsv"),
                        skip = 1, row.names = 1, check.names = FALSE)
asv_mat <- as.matrix(asv_table)
cat("ASV table:", nrow(asv_mat), "ASVs x", ncol(asv_mat), "samples\n")

# Read taxonomy
taxonomy <- read.delim(file.path(base_path, "taxonomy/taxonomy.tsv"),
                       row.names = 1)
# Parse taxonomy string into columns
tax_split <- strsplit(as.character(taxonomy$Taxon), "; ")
tax_mat <- do.call(rbind, lapply(tax_split, function(x) {
  x <- gsub("^[dpcofgs]__", "", x)  # Remove prefixes
  length(x) <- 7  # Pad to 7 levels
  x
}))
colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_mat) <- rownames(taxonomy)
cat("Taxonomy:", nrow(tax_mat), "taxa\n")

# Read metadata
metadata <- read.delim(file.path(base_path, "metadata.tsv"),
                       row.names = 1, check.names = FALSE)
cat("Metadata:", nrow(metadata), "samples x", ncol(metadata), "variables\n")

# Create phyloseq object
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(metadata)
)

cat("\nPhyloseq object created:\n")
print(ps)

# Save phyloseq object
saveRDS(ps, "phyloseq_object.rds")
cat("\n✅ Saved: phyloseq_object.rds\n")

# Also export alpha diversity metrics for figure_1e_alpha.R
alpha_obs <- read.delim(file.path(base_path, "observed/alpha-diversity.tsv"),
                        check.names = FALSE)
write_tsv(alpha_obs, "alpha-observed.tsv")
cat("✅ Saved: alpha-observed.tsv\n")

# Export distance matrix for figure_1g_beta.R
dm <- read.delim(file.path(base_path, "bray_curtis/distance-matrix.tsv"),
                 check.names = FALSE)
write_tsv(dm, "distance-matrix.tsv")
cat("✅ Saved: distance-matrix.tsv\n")

# Export metadata in expected format
meta_out <- metadata %>%
  rownames_to_column("sample-id")
write_tsv(meta_out, "metadata_merged.tsv")
write_tsv(meta_out, "metadata_merged_qiime_FINAL.tsv")
cat("✅ Saved: metadata_merged.tsv\n")
cat("✅ Saved: metadata_merged_qiime_FINAL.tsv\n")

cat("\n=== All input files created. Ready to run analysis scripts. ===\n")
