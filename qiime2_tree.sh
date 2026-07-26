#!/bin/bash
#SBATCH --job-name=qiime_make_tree
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/make_tree_%j.out
#SBATCH --error=logs/make_tree_%j.err

set -euo pipefail

echo "=============================================="
echo " QIIME2 PHYLOGENETIC TREE BUILD STARTED"
echo " Date: $(date)"
echo " Job ID: ${SLURM_JOB_ID:-NA}"
echo " Host: $(hostname)"
echo "=============================================="

# Load QIIME2
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

# Temp dir (good practice on BEAR)
export TMPDIR=/rds/projects/e/elhamsak-group6/tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR

# Go to your QIIME2 results folder
cd /rds/projects/e/elhamsak-group6/qiime2_results || exit 1
mkdir -p logs

# Sanity check
if [ ! -f rep-seqs.qza ]; then
  echo "❌ ERROR: rep-seqs.qza not found in $(pwd)"
  echo "Tree building needs rep-seqs.qza (representative sequences)."
  exit 1
fi

# Skip if already exists
if [ -f rooted-tree.qza ]; then
  echo "✅ rooted-tree.qza already exists. Skipping tree building."
  ls -lh rooted-tree.qza
  exit 0
fi

echo "🌳 Building phylogenetic tree from rep-seqs.qza ..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

echo "✅ Tree building complete. Outputs:"
ls -lh aligned-rep-seqs.qza masked-aligned-rep-seqs.qza unrooted-tree.qza rooted-tree.qza

echo "=============================================="
echo " QIIME2 PHYLOGENETIC TREE BUILD FINISHED"
echo " Date: $(date)"
echo "=============================================="


