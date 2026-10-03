# IF-P vs CR Microbiome Reproduction Study

[![QIIME2](https://img.shields.io/badge/QIIME2-2025.4-blue)](https://qiime2.org/)
[![R](https://img.shields.io/badge/R-4.0+-green)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Educational-yellow)](LICENSE)

## 📋 Overview

This repository contains the complete bioinformatics pipeline and analysis code for reproducing the gut microbiome analysis from:

> **Mohr AE, et al. (2024).** "Gut microbiome remodeling and metabolomic profile improves in response to protein pacing with intermittent fasting versus continuous caloric restriction." *Nature Communications*, 15, 4155.

### Study Summary

The original study compared two calorie-matched dietary interventions over 8 weeks in 41 participants with overweight/obesity:
- **IF-P (Intermittent Fasting with Protein Pacing):** 5:2 fasting pattern with 4-5 protein-rich meals/day
- **CR (Continuous Caloric Restriction):** Heart-healthy diet with 25% caloric deficit

**Key Finding:** IF-P produced distinct gut microbiome signatures independent of caloric restriction, with enrichment of beneficial bacteria like *Christensenellaceae*.

---

## 🔗 Links

| Resource | Link |
|----------|------|
| **Original Paper** | [Nature Communications](https://doi.org/10.1038/s41467-024-48355-5) |
| **Raw Data (SRA)** | [PRJNA847971](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA847971) |
| **Original Code** | [GitHub](https://github.com/Alex-E-Mohr/GM-Remodeling-IF-ProteinPacing-vs-CaloricRestriction) |

---

## 📁 Repository Structure

```
IF-P-vs-CR-Microbiome-Reproduction/
│
├── scripts/
│   ├── bash/                           # QIIME2 pipeline scripts
│   │   ├── download.sh                 # Download data from ENA
│   │   ├── fastqc.sh                   # Run FastQC quality control
│   │   ├── multiqc.sh                  # Aggregate QC reports
│   │   ├── qiime2_dada2.sh             # DADA2 denoising
│   │   ├── qiime2_taxo.sh              # Taxonomy classification (SILVA 138)
│   │   ├── qiime2_train_silva_138_v4.sh # Train V4-specific classifier
│   │   ├── qiime2_sepp.sh              # Phylogenetic tree (SEPP)
│   │   ├── qiime2_tree.sh              # Alternative tree building
│   │   ├── qiime2_alpha.sh             # Alpha diversity analysis
│   │   └── taxa_collapse.sh            # Collapse taxa at different levels
│   │
│   └── R/                              # R analysis scripts
│       ├── figure_generation.R         # Figures 1f (Faith's PD), 1g (PCoA), 1h (Taxa barplot)
│       ├── beneficial_bacteria.R       # Christensenellaceae/Rikenellaceae/Ruminococcaceae boxplots
│       ├── beta_trajectory.R           # PCoA with individual trajectories over time
│       └── machine_learning.R          # Random Forest classification (Top 10 features)
│
├── metadata/
│   ├── metadata.tsv                    # Original sample metadata
│   ├── metadata_merged.tsv             # Merged with ENA data
│   ├── metadata_merged_qiime_FINAL.tsv # QIIME2-formatted final metadata
│   ├── ena_file_info.tsv               # ENA file information
│   ├── checksums.md5                   # File integrity checksums
│   ├── SRR_Acc_List.txt                # SRA accession numbers
│   ├── common_samples_with_header.txt  # Sample filtering list
│   ├── common_samples.txt              # Common sample IDs
│   ├── download_urls.txt               # ENA download URLs
│   └── srr_sorted.txt                  # Sorted SRR accessions
│
├── results/
│   └── figures/                        # Generated figures
│
├── docs/
│   └── pipeline_workflow.md            # Detailed workflow documentation
│
├── .gitignore                          # Git ignore rules
└── README.md                           # Project documentation
```

---

## 🔬 Pipeline Overview

### Stage 1: Data Preparation

```mermaid
graph LR
    A[ENA Database] --> B[Download FASTQ]
    B --> C[FastQC/MultiQC]
    C --> D[Import to QIIME2]
    D --> E[demux.qza]
```

**Steps:**
1. Download 120 paired-end samples (240 FASTQ files) from ENA
2. Run quality control with FastQC
3. Aggregate reports with MultiQC
4. Import into QIIME2 using paired-end manifest

### Stage 2: QIIME2 Processing

```mermaid
graph LR
    A[demux.qza] --> B[DADA2]
    B --> C[Feature Table]
    B --> D[Rep Seqs]
    D --> E[Taxonomy<br>SILVA 138]
    D --> F[Phylogeny<br>SEPP]
    C --> G[Rarefaction<br>6,500 reads]
```

**Steps:**
1. **DADA2 Denoising:** Quality filter, error correct, merge pairs, remove chimeras
   - Trim: 20bp (primer removal)
   - Truncate: 240bp (quality threshold)
2. **Taxonomy Classification:** SILVA 138 database with V4-trained classifier
3. **Phylogenetic Tree:** SEPP fragment insertion
4. **Rarefaction:** Normalise to 6,500 reads/sample

### Stage 3: Diversity & Export

```mermaid
graph LR
    A[Rarefied Table] --> B[Alpha Diversity]
    A --> C[Beta Diversity]
    B --> D[Export to R]
    C --> D
    D --> E[Statistical Tests]
    E --> F[Figures]
```

**Steps:**
1. Calculate alpha diversity (Observed ASVs, Shannon, Faith's PD)
2. Calculate beta diversity (Bray-Curtis, UniFrac)
3. Export to R: ASV table, taxonomy, metadata
4. Statistical analysis and figure generation

---

## 💻 Requirements

### HPC Environment (BlueBear)

| Software | Version | Module |
|----------|---------|--------|
| QIIME2 | 2025.4 | `module load QIIME2/2025.4` |
| FastQC | 0.11.9+ | `module load FastQC` |
| MultiQC | 1.12+ | `module load MultiQC` |

### R Environment

```r
# Required packages
install.packages(c("tidyverse", "vegan", "ggplot2", "RColorBrewer"))

# Bioconductor packages
BiocManager::install(c("phyloseq", "DESeq2"))

# Optional for machine learning
install.packages(c("randomForest", "caret"))
```

---

## 🚀 Usage

### Running the QIIME2 Pipeline

```bash
# Navigate to project directory
cd /rds/projects/e/elhamsak-group6

# Load QIIME2 module
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

# Run pipeline scripts in order
sbatch scripts/bash/download.sh          # Day 1: Download data
sbatch scripts/bash/fastqc.sh            # Day 1: Quality control
sbatch scripts/bash/qiime2_dada2.sh      # Days 2-3: DADA2 (48 hours)
sbatch scripts/bash/qiime2_taxo.sh       # Day 4: Taxonomy
sbatch scripts/bash/qiime2_sepp.sh       # Day 4: Phylogenetic tree
sbatch scripts/bash/qiime2_alpha.sh      # Day 5: Diversity metrics
sbatch scripts/bash/export_for_R.sh      # Day 5: Export for R
```

### R Analysis

```r
# Figure Generation (Figure 1f, 1g, 1h)
# Requires: exported_for_R/ directory with faith_pd.tsv, bray_curtis_distance.tsv, 
#           feature_table.tsv, taxonomy.tsv, and metadata.tsv
source("scripts/R/figure_generation.R")

# Machine Learning (Random Forest)
# Requires: phyloseq_object.rds
source("scripts/R/machine_learning.R")
```

**Output Figures:**
- `Figure_1f_FaithPD.pdf/png` - Alpha diversity boxplot
- `Figure_1g_beta_diversity.pdf/png` - PCoA with 95% confidence ellipses
- `Figure_1h_taxa_barplot.pdf/png` - Top 10 genera relative abundance
- `Beta_Trajectory_PCoA.pdf/png` - PCoA with individual trajectories (WK0→WK4→WK8)
- `RF_Top10Features.pdf/png` - Random Forest feature importance
- `BeneficialBacteria_Boxplot.pdf/png` - Christensenellaceae/Rikenellaceae/Ruminococcaceae with p-values

---

## 📊 Key Results

| Metric | Value |
|--------|-------|
| **Samples processed** | 117 (after quality filtering) |
| **Read retention** | ~60-70% after DADA2 |
| **Rarefaction depth** | 6,500 reads/sample |
| **ASVs identified** | ~X,XXX unique sequences |

### Key Findings Reproduced:
- ✅ Beta diversity separation between IF-P and CR groups
- ✅ *Christensenellaceae* enrichment in IF-P
- ✅ Machine learning classification accuracy: ~79.5%

---

## 👥 Team Contributions

**Collaborative Project by:** Nubla Latheef, Rayane Baroudi, Farah Malwan

All team members contributed equally to:
- QIIME2 pipeline development and execution
- R statistical analysis and visualisation
- Machine learning validation
- Documentation and presentation

---

## 📚 References

1. Mohr AE, et al. (2024). Gut microbiome remodeling and metabolomic profile improves in response to protein pacing with intermittent fasting versus continuous caloric restriction. *Nature Communications*, 15, 4155.

2. Bolyen E, et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37, 852-857.

3. Callahan BJ, et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581-583.

4. Quast C, et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41, D590-D596.

---

## 🎓 Acknowledgments

- **Original Study Authors:** Mohr AE and colleagues
- **Institution:** University of Birmingham, Dubai
- **Course:** MSc Bioinformatics
- **Computing Resources:** BlueBear HPC, University of Birmingham

---

## 📄 License

This project is for educational purposes as part of coursework at the University of Birmingham. The original data and methodology belong to Mohr et al. (2024).

---

## 📧 Contact

For questions about this reproduction study:
- fnl759@student.bham.ac.uk

For questions about the original study:
- See original publication contact information
