# =============================================================
# LOESS_Normalization_Analysis.R
# ChIP-seq LOESS Normalization — Watanabe Lab Case Study
# Author: Addison Yam | March 2026
#
# Fits a LOESS model to the input-ChIP relationship in clean background
# windows, then subtracts predicted background from observed ChIP signal.
# The residual represents genuine histone enrichment above background.
#
# INPUT:  SampleA_ChIP.tab, SampleA_input.tab,
#         windows_in_peaks.bed, windows_in_blacklist.bed
# OUTPUT: SampleA_residuals.tsv, four diagnostic plots
# =============================================================

library(data.table)
library(ggplot2)


# --- LOAD DATA ------------------------------------------------

message("Loading ChIP and input tab files...")

chip <- fread("SampleA_ChIP.tab",
              header = FALSE,
              col.names = c("name", "size", "covered", "sum", "mean0", "mean"))

input <- fread("SampleA_input.tab",
               header = FALSE,
               col.names = c("name", "size", "covered", "sum", "mean0", "mean"))

message(sprintf("Windows loaded: %d", nrow(chip)))

# Confirm both files are aligned row-for-row before any calculations
stopifnot(all(chip$name == input$name))
message("Row alignment check passed.")


# --- ANNOTATE PEAKS AND BLACKLIST -----------------------------

peaks_windows     <- fread("windows_in_peaks.bed",     header = FALSE,
                           col.names = c("chr", "start", "end", "name"))
blacklist_windows <- fread("windows_in_blacklist.bed",  header = FALSE,
                           col.names = c("chr", "start", "end", "name"))

is_peak      <- chip$name %in% peaks_windows$name
is_blacklist <- chip$name %in% blacklist_windows$name

message(sprintf("Peak windows:      %d (%.2f%%)", sum(is_peak),
                100 * sum(is_peak) / nrow(chip)))
message(sprintf("Blacklist windows: %d (%.2f%%)", sum(is_blacklist),
                100 * sum(is_blacklist) / nrow(chip)))
message(sprintf("Clean background:  %d (%.2f%%)", sum(!is_peak & !is_blacklist),
                100 * sum(!is_peak & !is_blacklist) / nrow(chip)))


# --- BUILD LOESS TRAINING SET --------------------------------
# Train only on clean background (no peaks, no blacklist).
# log1p applied to input to compress right-skewed signal range.
# Residuals stay in original ChIP space so log1p is not applied to ChIP.

bg_indices   <- which(!is_peak & !is_blacklist)
chip_signal  <- chip$mean0[bg_indices]
input_signal <- input$mean0[bg_indices]
log_input    <- log1p(input_signal)

# Collapse to unique (chip, input) pairs weighted by frequency.
# Reduces 30M rows to ~22K unique pairs — equivalent fit, much faster.
bg_dt       <- data.table(ch = chip_signal, inpt = log_input)
weighted_bg <- bg_dt[, .(w = .N), by = .(ch, inpt)]

message(sprintf("Unique training pairs after collapsing: %d", nrow(weighted_bg)))


# --- FIT LOESS MODEL -----------------------------------------
# span=1.0 : smooth global fit using all data for each local estimate
# degree=1 : linear local fits (faster, appropriate here)
# surface="interpolate" : speeds up prediction at new points

message("Fitting LOESS model...")

model <- loess(ch ~ inpt,
               data    = weighted_bg,
               weights = w,
               span    = 1.0,
               degree  = 1,
               control = loess.control(surface = "interpolate"))

message(sprintf("Residual standard error: %.6f", sqrt(model$s)))


# --- PREDICT AND CALCULATE RESIDUALS -------------------------

log_input_all <- log1p(input$mean0)
preds         <- predict(model, newdata = data.frame(inpt = log_input_all))
residuals     <- pmax(0, chip$mean0 - preds)   # floor at zero; negatives are noise

message(sprintf("Windows with residual > 0: %d (%.2f%%)",
                sum(residuals > 0), 100 * sum(residuals > 0) / length(residuals)))
message(sprintf("Max residual: %.4f | Mean non-zero: %.4f",
                max(residuals, na.rm = TRUE),
                mean(residuals[residuals > 0], na.rm = TRUE)))


# --- SAVE OUTPUT ---------------------------------------------

output_dt <- data.table(
  name      = chip$name,
  chr       = sub("_[0-9]+_[0-9]+$", "", chip$name),
  start     = as.integer(sub(".*_([0-9]+)_[0-9]+$", "\\1", chip$name)),
  end       = as.integer(sub(".*_[0-9]+_([0-9]+)$", "\\1", chip$name)),
  chip      = chip$mean0,
  input     = input$mean0,
  log_input = log_input_all,
  predicted = preds,
  residual  = residuals
)

fwrite(output_dt, "SampleA_residuals.tsv", sep = "\t")
message("Saved: SampleA_residuals.tsv")


# --- DIAGNOSTIC PLOTS ----------------------------------------
# Sample 50k background windows for visualization; 30M points is impractical.

set.seed(42)
plot_idx <- sample(which(!is_peak & !is_blacklist), 50000)

plot_dt <- data.table(
  log_input = log_input_all[plot_idx],
  chip      = chip$mean0[plot_idx],
  predicted = preds[plot_idx],
  residual  = residuals[plot_idx]
)

# Plot 1 — LOESS fit: curve should sit near zero and rise slightly at high input
curve_x  <- seq(min(plot_dt$log_input), max(plot_dt$log_input), length.out = 500)
curve_y  <- predict(model, newdata = data.frame(inpt = curve_x))
curve_dt <- data.table(log_input = curve_x, predicted = curve_y)

p1 <- ggplot(plot_dt, aes(x = log_input, y = chip)) +
  geom_point(alpha = 0.05, size = 0.3, color = "steelblue") +
  geom_line(data = curve_dt, aes(x = log_input, y = predicted),
            color = "red", linewidth = 1.2) +
  labs(title = "LOESS Fit: ChIP Signal vs log1p(Input)",
       subtitle = "Red = LOESS curve | Blue = background windows (n=50,000 sampled)",
       x = "log1p(Input Signal)", y = "ChIP Signal (mean0)") +
  theme_bw()

ggsave("plot1_loess_fit.png", p1, width = 8, height = 6, dpi = 150)

# Plot 2 — Residual distribution: should be smooth; spikes indicate discretization
p2 <- ggplot(plot_dt[residual > 0], aes(x = residual)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white", linewidth = 0.2) +
  scale_x_log10() +
  labs(title = "Distribution of Positive Residuals (Background Windows)",
       x = "Residual (log10 scale)", y = "Count") +
  theme_bw()

ggsave("plot2_residual_distribution.png", p2, width = 8, height = 6, dpi = 150)

# Plot 3 — Predicted vs actual: most points near zero; vertical streaks = real peaks
p3 <- ggplot(plot_dt, aes(x = predicted, y = chip)) +
  geom_point(alpha = 0.05, size = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1.0, linetype = "dashed") +
  labs(title = "Predicted vs Actual ChIP Signal (Background Windows)",
       x = "LOESS Predicted Signal", y = "Observed ChIP Signal") +
  theme_bw()

ggsave("plot3_predicted_vs_actual.png", p3, width = 8, height = 6, dpi = 150)

# Plot 4 — Peak vs background: peak distribution must sit clearly right of background
comparison_dt <- data.table(
  residual = c(
    residuals[sample(which(!is_peak & !is_blacklist), 50000)],
    residuals[sample(which(is_peak), min(50000, sum(is_peak)))]
  ),
  group = c(rep("Background", 50000), rep("Peak", min(50000, sum(is_peak))))
)

p4 <- ggplot(comparison_dt, aes(x = residual, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = c("Background" = "steelblue", "Peak" = "firebrick")) +
  labs(title = "Residual Distribution: Peak vs Background Windows",
       x = "Residual (log10 scale)", y = "Density", fill = "Region") +
  theme_bw()

ggsave("plot4_peak_vs_background.png", p4, width = 8, height = 6, dpi = 150)

message("All plots saved. Pipeline complete.")
