#!/bin/bash
# =============================================================
# preprocessing.sh
# ChIP-seq Signal Extraction Pipeline — Watanabe Lab Case Study
# Author: Addison Yam | March 2026
#
# Prepares genome-wide signal tables for LOESS normalization:
#   1. Filter genome reference to canonical chromosomes
#   2. Tile genome into 100bp windows
#   3. Extract ChIP and input signal per window
#   4. Annotate windows overlapping peaks or blacklist regions
#
# REQUIRES: bedtools/2.31.0, ucscutils/2018-11-27
# =============================================================

module load bedtools/2.31.0
module load ucscutils/2018-11-27

echo "Starting preprocessing pipeline..."


# --- STEP 1: Filter to canonical chromosomes ---
# The original file has 93 entries including alt haplotypes and unplaced
# scaffolds. We keep only chr1-22, chrX, chrY, chrM (25 total).

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y|M)$/ && $1 !~ /[_.]/' \
    hg19.chrom.sizes > hg19.canonical.chrom.sizes

echo "  Canonical chromosomes kept: $(wc -l < hg19.canonical.chrom.sizes)"


# --- STEP 2: Tile genome into 100bp windows ---

bedtools makewindows \
    -g hg19.canonical.chrom.sizes \
    -w 100 \
    > 100bp.windows.bed

echo "  Windows generated: $(wc -l < 100bp.windows.bed)"


# --- STEP 3: Add coordinate-based names to each window ---
# bigWigAverageOverBed requires a name column. We use coordinates
# (e.g. chr1_0_100) so names are unique and parseable downstream.

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1"_"$2"_"$3}' \
    100bp.windows.bed > 100bp.windows.named.bed


# --- STEP 4: Extract signal from bigWig files ---
# mean0 = average over all 100 bases with zeros included.
# Preferred over mean because uncovered bases are true zeros, not missing data.

echo "Extracting ChIP signal..."
time bigWigAverageOverBed \
    SampleA_ChIP.bw \
    100bp.windows.named.bed \
    SampleA_ChIP.tab

echo "Extracting input signal..."
time bigWigAverageOverBed \
    SampleA_input.bw \
    100bp.windows.named.bed \
    SampleA_input.tab

# Verify row counts match before passing to R
chip_rows=$(wc -l < SampleA_ChIP.tab)
input_rows=$(wc -l < SampleA_input.tab)

if [ "$chip_rows" -eq "$input_rows" ]; then
    echo "  Row count check: PASSED ($chip_rows rows)"
else
    echo "  ERROR: Row counts do not match. Exiting."
    exit 1
fi


# --- STEP 5: Annotate peak and blacklist windows ---
# The LOESS model trains on clean background only — peaks are excluded
# so real enrichment doesn't get absorbed into the background model.
# -u reports each window once if it overlaps anything in -b.

bedtools intersect \
    -a 100bp.windows.named.bed \
    -b SampleA.narrowPeak.bed \
    -u > windows_in_peaks.bed

bedtools intersect \
    -a 100bp.windows.named.bed \
    -b blacklist_hg19.bed \
    -u > windows_in_blacklist.bed

echo "  Peak windows:      $(wc -l < windows_in_peaks.bed)"
echo "  Blacklist windows: $(wc -l < windows_in_blacklist.bed)"
echo "Preprocessing complete."
