# ChIP-seq LOESS Normalization with Copy Number Correction

## Overview
This pipeline normalizes ChIP-seq signal using LOESS regression trained on input control, then applies **input-based copy number normalization** to correct for oncogene amplification artifacts in cancer cell lines.

## Why copy number matters
In SCLC cell line H69, MYCN is amplified ~70x. Without correction, ChIP signal at this locus is inflated by copy number, not biology. This pipeline estimates copy number from input signal and scales ChIP signal accordingly.

## Key results
- MYCN residuals reduced by **87%** after CN normalization
- MYC residuals reduced by **90%**
- Diploid genome changed by only **4%** (correction is targeted)

## Files
| File | Purpose |
|------|---------|
| `Preprocessing.sh` | 100bp windows, signal extraction from bigWigs |
| `LOESS_Normalization_Analysis.R` | LOESS fit, residual calculation |
| `Anomaly2_CN_Normalization.R` | Copy number estimation and correction |

## Dependencies
- bedtools v2.31.0
- UCSC bigWigAverageOverBed
- R 4.4.0 (data.table, ggplot2)

## Citation note
This method is an approximation. For precise copy number correction, use WGS data with tools like CNVkit or sequenza.
