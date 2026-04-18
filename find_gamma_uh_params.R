#!/usr/bin/env Rscript
# find_gamma_uh_params.R
# Fit gamma distribution parameters to unit hydrograph ordinates from a FEWS XML file.
# Usage: Rscript find_gamma_uh_params.R <input_xml>

box::use(
  fitdistrplus[fitdist],
  nwsrfsr[uh2p_cfs_in]
)
source("fews_uh_pars.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript find_gamma_uh_params.R <input_xml>")
}

input_file <- args[1]
uhg <- get_uhg_params(input_file)

# Extract ordinates, interval, and zone area from the returned list
y_values     <- uhg[[5]]$ordinates
interval_val <- uhg$uhg_interval
zone_area    <- uhg$drainage_area

# Compute interval series
interval_series <- seq_len(length(y_values)) * interval_val

cat("Original ordinates:\n")
print(data.frame(interval = interval_series, ordinate = y_values))
cat("\n")

# Scale ordinates down by the cfs/in conversion factor
y_values <- y_values / (5280^2 * 24 / 12 / 86400 / interval_val * zone_area)

# Inverse transform sampling: convert PDF ordinates to a gamma-distributed sample dataset
# Step 1: normalize ordinates to a proper PDF (area = 1)
pdf_vals <- y_values / sum(y_values)

# Step 2: compute CDF and ensure last value is exactly 1
cdf_vals <- cumsum(pdf_vals)
cdf_vals <- cdf_vals / max(cdf_vals)

# Step 3: inverse transform sampling — map uniform random numbers through the CDF
set.seed(123)
n_samples <- 10000
u       <- runif(n_samples)
samples <- approx(cdf_vals, interval_series, xout = u)$y
samples <- samples[!is.na(samples)]

# Step 4: fit gamma distribution to the generated sample dataset using MLE
fit   <- fitdist(samples, "gamma", method = "mle")
shape <- unname(fit$estimate["shape"])   # n (shape parameter)
rate  <- unname(fit$estimate["rate"])
scale <- 1 / rate                        # K (scale parameter)


# Compute UH ordinates from fitted parameters and derive time of concentration
# scale is in hours from the fit; uh2p_cfs_in expects scale in days
uh_ords <- uh2p_cfs_in(shape, scale / 24, interval_val, zone_area)
cat("Fitted UH ordinates:\n")
print(data.frame(interval = seq_len(length(uh_ords)) * interval_val, ordinate = uh_ords))
cat("\n")
n_ords  <- length(uh_ords)
toc     <- (n_ords - 1) * interval_val

cat("Fitted gamma unit hydrograph parameters:\n")
cat(sprintf("  Shape (n):                %.6f\n", shape))
cat(sprintf("  Scale (K) (days):         %.6f\n", scale / 24))
cat(sprintf("  Number of UH ordinates:   %d\n",   n_ords))
cat(sprintf("  Time of concentration:    %.2f hours\n", toc))

# Write output CSV next to the input file
output_file <- file.path(
  dirname(input_file),
  sub("\\.xml$", "_gamma_params.csv", basename(input_file))
)
output_df <- data.frame(
  parameter = c("shape_n", "scale_K", "toc"),
  value     = c(shape, scale / 24, toc)
)
write.csv(output_df, output_file, row.names = FALSE)
cat(sprintf("\nParameters saved to: %s\n", output_file))
