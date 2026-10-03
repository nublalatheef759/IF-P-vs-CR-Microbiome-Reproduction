library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)

# ============================================
# FIGURE 1f: Alpha Diversity (Faith's PD)
# ============================================

# Load Faith's PD from QIIME2
faith_pd <- read.delim("exported_for_R/faith_pd.tsv", row.names = 1)
metadata <- read.delim("metadata.tsv", row.names = 1)

# Merge
alpha_data <- merge(faith_pd, metadata, by = "row.names")

# Box plot
p_alpha <- ggplot(alpha_data, aes(x = Group, y = faith_pd, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    y = "Faith's Phylogenetic Diversity",
    x = "",
    title = "Figure 1f: Alpha Diversity"
  ) +
  scale_fill_manual(values = c("#5DADE2", "#AF7AC5"))

ggsave("Figure_1f_FaithPD.pdf", plot = p_alpha, width = 6, height = 6, dpi = 300)
ggsave("Figure_1f_FaithPD.png", plot = p_alpha, width = 6, height = 6, dpi = 300)
print("✅ Figure 1f created!")

# ============================================
# FIGURE 1g: Beta Diversity (PCoA)
# ============================================

# Load distance matrix
bray_dist <- read.delim("exported_for_R/bray_curtis_distance.tsv",
                        row.names = 1, check.names = FALSE)

# Convert to distance object
bray_dist_matrix <- as.dist(bray_dist)

# Run PCoA
pcoa <- cmdscale(bray_dist_matrix, k = 2, eig = TRUE)

# Calculate variance explained
var_exp <- pcoa$eig / sum(pcoa$eig) * 100

# Create data frame
pcoa_data <- data.frame(
  Sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

# Merge with metadata
pcoa_data <- merge(pcoa_data, metadata, by.x = "Sample", by.y = "row.names")

# Plot
p_beta <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    title = "Figure 1g: Beta Diversity (PCoA)"
  ) +
  scale_color_manual(values = c("#5DADE2", "#AF7AC5"))

ggsave("Figure_1g_beta_diversity.pdf", plot = p_beta, width = 8, height = 6, dpi = 300)
ggsave("Figure_1g_beta_diversity.png", plot = p_beta, width = 8, height = 6, dpi = 300)

# PERMANOVA test
adonis_result <- adonis2(bray_dist_matrix ~ Group, data = metadata)
print(adonis_result)
print("✅ Figure 1g created!")
