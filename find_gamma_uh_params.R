#!/usr/bin/env Rscript
# find_gamma_uh_params.R
# Fit gamma distribution parameters to unit hydrograph ordinates from a pars CSV file.
# Usage: Rscript find_gamma_uh_params.R <input_csv>

box::use(
  fitdistrplus[fitdist],
  nwsrfsr[uh2p_cfs_in]
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript find_gamma_uh_params.R <input_csv>")
}

input_file <- args[1]
pars <- read.csv(input_file, stringsAsFactors = FALSE)

# Extract UH ordinates, ordered by their sequence number (last digits of name)
uh_rows <- pars[grepl("^unit_ord_", pars$name) & pars$type == "uh", ]
seq_nums <- as.integer(sub("^unit_ord_", "", uh_rows$name))
uh_rows  <- uh_rows[order(seq_nums), ]
y_values <- uh_rows$value


# Extract interval and compute interval series
interval_val <- pars$value[pars$name == "interval" & pars$type == "uh"]
interval_series <- seq_len(length(y_values)) * interval_val

cat("interval values:", interval_series, "\n")
cat("y values:", y_values, "\n\n")

# Read zone area
zone_area <- pars$value[pars$name == "zone_area" & pars$type == "uh"]

# Scale ordinates down by the cfs/in conversion factor
y_values <- y_values / (5280^2 * 24 / 12 / 86400 / interval_val * zone_area)

# Fit gamma distribution using MLE on the UH ordinates (density values) directly
# replace any remaining zeros with 0.1 (fitdist requires strictly positive data)
s   <- ifelse(y_values == 0, 0.1, y_values)
fit <- fitdist(s, "gamma", method="mle")
shape <- unname(fit$estimate["shape"])   # n (shape parameter)
rate  <- unname(fit$estimate["rate"])
scale <- 1 / rate                        # K (scale parameter)
plot(fit)

# Compute UH ordinates from fitted parameters and derive time of concentration
uh_ords <- uh2p_cfs_in(shape, scale, interval_val, zone_area)
cat("UH ordinates:", uh_ords, "\n\n")
n_ords  <- length(uh_ords)
toc     <- (n_ords - 1) * interval_val

cat("Fitted gamma unit hydrograph parameters:\n")
cat(sprintf("  Shape (n):                %.6f\n", shape))
cat(sprintf("  Scale (K):                %.6f\n", scale))
cat(sprintf("  Number of UH ordinates:   %d\n",   n_ords))
cat(sprintf("  Time of concentration:    %.2f hours\n", toc))

# Write output CSV next to the input file
output_file <- file.path(
  dirname(input_file),
  sub("\\.csv$", "_gamma_params.csv", basename(input_file))
)
output_df <- data.frame(
  parameter = c("shape_n", "scale_K", "toc"),
  value     = c(shape, scale, toc)
)
write.csv(output_df, output_file, row.names = FALSE)
cat(sprintf("\nParameters saved to: %s\n", output_file))
