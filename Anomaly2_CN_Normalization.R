# =============================================================
# Anomaly2_CN_Normalization.R
# Copy Number Normalization & Amplicon Flagging — Watanabe Lab Case Study
# Author: Addison Yam | March 2026
#
# Addresses extreme residuals at MYCN and MYC loci caused by oncogene
# amplification in the H69 SCLC cell line. Amplified regions have more
# DNA copies, which inflates ChIP signal independently of histone biology.
#
# This script:
#   1. Flags windows within known amplicon boundaries
#   2. Estimates copy number from input signal relative to diploid baseline
#   3. Divides ChIP signal by CN factor and recalculates residuals
#
# Note: CN estimation from input is an approximation. Input is also affected
# by chromatin accessibility and GC bias, not just copy number. WGS-based
# copy number calls would give more precise correction.
#
# INPUT:  SampleA_residuals.tsv (output of loess_analysis.R)
# OUTPUT: SampleA_residuals_cn_normalized.tsv, three diagnostic plots
# =============================================================

library(data.table)
library(ggplot2)


# --- LOAD DATA ------------------------------------------------

message("Loading residuals table...")
dt <- fread("SampleA_residuals.tsv")
message(sprintf("Loaded %d windows.", nrow(dt)))


# --- FLAG KNOWN AMPLICON WINDOWS ------------------------------
# MYCN (chr2:15.5-18Mb) and MYC region (chr2:134-136Mb) were identified
# as extreme outliers in the LOESS output — top residuals in the genome.
# H69 is a documented MYCN-amplified SCLC cell line.

amplicons <- data.table(
  chr   = c("chr2",          "chr2"),
  start = c(15500000,        134000000),
  end   = c(18000000,        136000000),
  label = c("MYCN_amplicon", "MYC_amplicon")
)

dt[, amplicon := "none"]

for (i in 1:nrow(amplicons)) {
  dt[chr == amplicons$chr[i] &
     start >= amplicons$start[i] &
     end   <= amplicons$end[i],
     amplicon := amplicons$label[i]]
}

message(sprintf("MYCN windows flagged: %d", sum(dt$amplicon == "MYCN_amplicon")))
message(sprintf("MYC windows flagged:  %d", sum(dt$amplicon == "MYC_amplicon")))


# --- ESTIMATE COPY NUMBER FROM INPUT --------------------------
# Diploid baseline = median non-zero input signal across the genome.
# CN factor = input / baseline, floored at 1.0 (no downward correction).
# Amplified regions will have input signal well above the baseline.

diploid_baseline <- median(dt$input[dt$input > 0])
message(sprintf("Diploid baseline: %.4f", diploid_baseline))

dt[, cn_factor := pmax(1.0, input / diploid_baseline)]

message("\nCN factor summary — MYCN amplicon:")
print(summary(dt[amplicon == "MYCN_amplicon"]$cn_factor))

message("\nCN factor summary — MYC amplicon:")
print(summary(dt[amplicon == "MYC_amplicon"]$cn_factor))


# --- APPLY CN NORMALIZATION ----------------------------------
# Divide ChIP by CN factor to scale signal back to the diploid equivalent.
# Recalculate residual using CN-normalized ChIP against the original
# LOESS-predicted background (which was trained on diploid regions).

dt[, chip_cn_normalized := chip / cn_factor]
dt[, residual_cn := pmax(0, chip_cn_normalized - predicted)]

# Diploid genome should change very little (<5%) — larger change would
# indicate the correction is distorting normal signal.
message("\nBefore vs after — MYCN:")
message(sprintf("  Mean residual: %.2f -> %.2f",
                mean(dt[amplicon == "MYCN_amplicon"]$residual),
                mean(dt[amplicon == "MYCN_amplicon"]$residual_cn)))

message("Before vs after — MYC:")
message(sprintf("  Mean residual: %.2f -> %.2f",
                mean(dt[amplicon == "MYC_amplicon"]$residual),
                mean(dt[amplicon == "MYC_amplicon"]$residual_cn)))

message("Before vs after — diploid genome:")
message(sprintf("  Mean residual: %.4f -> %.4f",
                mean(dt[amplicon == "none"]$residual),
                mean(dt[amplicon == "none"]$residual_cn)))


# --- SAVE OUTPUT ---------------------------------------------

fwrite(dt, "SampleA_residuals_cn_normalized.tsv", sep = "\t")
message("Saved: SampleA_residuals_cn_normalized.tsv")


# --- DIAGNOSTIC PLOTS ----------------------------------------

# Plot 1 — CN factor distribution: most genome diploid (tall bar at 1),
# long right tail = amplified regions. Red line at 2x for reference.
p1 <- ggplot(dt[cn_factor > 1], aes(x = cn_factor)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white", linewidth = 0.2) +
  scale_x_log10() +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of CN Factors (windows above diploid baseline)",
       subtitle = "Red dashed line = 2x (single-copy amplification)",
       x = "CN Factor (log10 scale)", y = "Count") +
  theme_bw()

ggsave("plot_cn_factor_distribution.png", p1, width = 8, height = 6, dpi = 150)

# Plot 2 — Before vs after at each amplicon: red right tail should shrink after correction
amplicon_long <- melt(
  dt[amplicon != "none", .(residual_before = residual, residual_after = residual_cn, amplicon)],
  id.vars = "amplicon", measure.vars = c("residual_before", "residual_after"),
  variable.name = "normalization", value.name = "residual"
)
amplicon_long[, normalization := ifelse(normalization == "residual_before",
                                        "Before", "After")]

p2 <- ggplot(amplicon_long[residual > 0], aes(x = residual, fill = normalization)) +
  geom_density(alpha = 0.6) +
  scale_x_log10() +
  scale_fill_manual(values = c("Before" = "firebrick", "After" = "steelblue")) +
  facet_wrap(~amplicon) +
  labs(title = "Residuals at Amplicon Loci Before vs After CN Normalization",
       x = "Residual (log10 scale)", y = "Density", fill = "") +
  theme_bw()

ggsave("plot_amplicon_before_after.png", p2, width = 10, height = 6, dpi = 150)

# Plot 3 — MYCN vs diploid: after correction, MYCN (salmon) should shift toward
# the diploid range (blue). Diploid distributions should barely move.
set.seed(42)
sample_idx <- sample(which(dt$amplicon == "none" & dt$residual > 0), 50000)

comparison_dt <- rbind(
  data.table(residual = dt$residual[sample_idx],    group = "Diploid (before)"),
  data.table(residual = dt$residual_cn[sample_idx], group = "Diploid (after)"),
  data.table(residual = dt[amplicon == "MYCN_amplicon"]$residual,    group = "MYCN (before)"),
  data.table(residual = dt[amplicon == "MYCN_amplicon"]$residual_cn, group = "MYCN (after)")
)

p3 <- ggplot(comparison_dt[residual > 0], aes(x = residual, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = c(
    "Diploid (before)" = "steelblue", "Diploid (after)"  = "cornflowerblue",
    "MYCN (before)"    = "firebrick", "MYCN (after)"     = "salmon"
  )) +
  labs(title = "Effect of CN Normalization: MYCN vs Diploid Regions",
       x = "Residual (log10 scale)", y = "Density", fill = "") +
  theme_bw()

ggsave("plot_mycn_vs_diploid.png", p3, width = 9, height = 6, dpi = 150)

message("All done.")
