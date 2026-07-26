#!/bin/bash
#SBATCH --job-name=export_data
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/export_%j.out
#SBATCH --error=logs/export_%j.err

set -euo pipefail

echo "=== Data Export for R Started: $(date) ==="
echo "Job ID: $SLURM_JOB_ID"

cd /rds/projects/e/elhamsak-group6/qiime2_results || exit 1

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

qiime --version

mkdir -p exported-for-R

############################################
# 1. Export ASV table (rarefied + matched to metadata)
############################################

echo ""
echo "=== Exporting ASV table ==="

qiime tools export \
  --input-path table-rarefied-6500-matched.qza \
  --output-path exported-for-R

echo "Converting BIOM to TSV..."

biom convert \
  -i exported-for-R/feature-table.biom \
  -o exported-for-R/asv-table.tsv \
  --to-tsv

echo "✅ ASV table exported"

############################################
# 2. Export taxonomy
############################################

echo ""
echo "=== Exporting taxonomy ==="

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-for-R/taxonomy

cp exported-for-R/taxonomy/taxonomy.tsv exported-for-R/taxonomy.tsv

echo "✅ Taxonomy exported"

############################################
# 3. Export rooted phylogenetic tree
############################################

echo ""
echo "=== Exporting phylogenetic tree ==="

qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-for-R/tree

cp exported-for-R/tree/tree.nwk exported-for-R/tree.nwk

echo "✅ Tree exported"

############################################
# 4. Export representative sequences
############################################

echo ""
echo "=== Exporting representative sequences ==="

qiime tools export \
  --input-path rep-seqs.qza \
  --output-path exported-for-R/rep-seqs

cp exported-for-R/rep-seqs/dna-sequences.fasta exported-for-R/rep-seqs.fasta

echo "✅ Representative sequences exported"

############################################
# 5. Export alpha diversity
############################################

echo ""
echo "=== Exporting alpha diversity ==="

# Faith's PD
qiime tools export \
  --input-path diversity-results/faith_pd_vector.qza \
  --output-path exported-for-R/alpha

cp exported-for-R/alpha/alpha-diversity.tsv exported-for-R/faith_pd.tsv

# Shannon
qiime tools export \
  --input-path diversity-results/shannon_vector.qza \
  --output-path exported-for-R/shannon

cp exported-for-R/shannon/alpha-diversity.tsv exported-for-R/shannon.tsv

# Observed features
qiime tools export \
  --input-path diversity-results/observed_features_vector.qza \
  --output-path exported-for-R/observed

cp exported-for-R/observed/alpha-diversity.tsv exported-for-R/observed_features.tsv

# Evenness
qiime tools export \
  --input-path diversity-results/evenness_vector.qza \
  --output-path exported-for-R/evenness

cp exported-for-R/evenness/alpha-diversity.tsv exported-for-R/evenness.tsv

echo "✅ Alpha diversity exported"

############################################
# 6. Export beta diversity matrices
############################################

echo ""
echo "=== Exporting beta diversity matrices ==="

# Unweighted UniFrac
qiime tools export \
  --input-path diversity-results/unweighted_unifrac_distance_matrix.qza \
  --output-path exported-for-R/unweighted_unifrac

cp exported-for-R/unweighted_unifrac/distance-matrix.tsv \
   exported-for-R/unweighted_unifrac_distance.tsv

# Weighted UniFrac
qiime tools export \
  --input-path diversity-results/weighted_unifrac_distance_matrix.qza \
  --output-path exported-for-R/weighted_unifrac

cp exported-for-R/weighted_unifrac/distance-matrix.tsv \
   exported-for-R/weighted_unifrac_distance.tsv

# Bray-Curtis
qiime tools export \
  --input-path diversity-results/bray_curtis_distance_matrix.qza \
  --output-path exported-for-R/bray_curtis

cp exported-for-R/bray_curtis/distance-matrix.tsv \
   exported-for-R/bray_curtis_distance.tsv

# Jaccard
qiime tools export \
  --input-path diversity-results/jaccard_distance_matrix.qza \
  --output-path exported-for-R/jaccard

cp exported-for-R/jaccard/distance-matrix.tsv \
   exported-for-R/jaccard_distance.tsv

echo "✅ Beta diversity matrices exported"

############################################
# 7. Export PCoA results
############################################

echo ""
echo "=== Exporting PCoA results ==="

qiime tools export \
  --input-path diversity-results/unweighted_unifrac_pcoa_results.qza \
  --output-path exported-for-R/pcoa

cp exported-for-R/pcoa/ordination.txt exported-for-R/pcoa_results.txt

echo "✅ PCoA results exported"

############################################
# 8. Copy metadata
############################################

echo ""
echo "=== Copying metadata ==="

cp ../metadata.tsv exported-for-R/metadata.tsv

echo "✅ Metadata copied"

############################################
# Summary
############################################

echo ""
echo "=== Export Complete: $(date) ==="
echo ""
echo "Files created in exported-for-R/:"
echo ""

ls -lh exported-for-R/*.tsv exported-for-R/*.nwk exported-for-R/*.fasta exported-for-R/*.txt 2>/dev/null | head -20

echo ""
echo "=== File count ==="
ls exported-for-R/*.tsv exported-for-R/*.nwk exported-for-R/*.fasta exported-for-R/*.txt 2>/dev/null | wc -l

echo ""
echo "✅ All data exported successfully!"
echo ""
echo "⚠️  NEXT STEPS:"
echo "1. Download to your laptop:"
echo "   scp -r fn1759@bluebear.bham.ac.uk:/rds/projects/e/elhamsak-group6/qiime2_results/exported-for-R ."
echo ""
echo "2. Optional - compress for faster download:"
echo "   cd /rds/projects/e/elhamsak-group6/qiime2_results"
echo "   tar -czf exported-for-R.tar.gz exported-for-R/"
echo "   # Then download:"
echo "   scp fn1759@bluebear.bham.ac.uk:/rds/projects/e/elhamsak-group6/qiime2_results/exported-for-R.tar.gz ."
echo "   tar -xzf exported-for-R.tar.gz"
