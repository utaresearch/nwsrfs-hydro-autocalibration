#!/usr/bin/env Rscript

box::use(
  ./fews_lagk_pars[get_lagk_params]
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript find_lagk_params.R <pixml_file>")
}

parpixmlfile <- args[1]
if (!file.exists(parpixmlfile)) {
  stop(paste("File not found:", parpixmlfile))
}

# Read the PIXML file
lagk      <- get_lagk_params(parpixmlfile)
lagq_pairs <- lagk$lagq_pairs
kq_pairs   <- lagk$kq_pairs

# ---------------------------------------------------------------------------
# Lag-Q fitting
# ---------------------------------------------------------------------------
Lag_max <- max(lagq_pairs$lagq_pairs_lags)
Lag_min <- min(lagq_pairs$lagq_pairs_lags)

Q_max <- max(lagq_pairs$lagq_pairs_qs)
Q_min <- min(lagq_pairs$lagq_pairs_qs)

# Normalize Lag and Q values
normalized_lag <- data.frame(
  x = (lagq_pairs$lagq_pairs_qs  - Q_min) / (Q_max - Q_min),
  y = (lagq_pairs$lagq_pairs_lags - Lag_min) / (Lag_max - Lag_min)
)

lag_fit   <- lm(y ~ x + I(x^2), data = normalized_lag)
lag_coefs <- coef(lag_fit)
lag_c     <- lag_coefs["(Intercept)"]
lag_b     <- lag_coefs["x"]
lag_a     <- lag_coefs["I(x^2)"]

# ---------------------------------------------------------------------------
# K-Q fitting
# ---------------------------------------------------------------------------
K_max <- max(kq_pairs$kq_pairs_ks)
K_min <- min(kq_pairs$kq_pairs_ks)

kq_Q_max <- max(kq_pairs$kq_pairs_qs)
kq_Q_min <- min(kq_pairs$kq_pairs_qs)

# Normalize K and Q values
normalized_k <- data.frame(
  x = (kq_pairs$kq_pairs_qs - kq_Q_min) / (kq_Q_max - kq_Q_min),
  y = (kq_pairs$kq_pairs_ks  - K_min)   / (K_max   - K_min)
)

k_fit   <- lm(y ~ x + I(x^2), data = normalized_k)
k_coefs <- coef(k_fit)
k_c     <- k_coefs["(Intercept)"]
k_b     <- k_coefs["x"]
k_a     <- k_coefs["I(x^2)"]

# ---------------------------------------------------------------------------
# Print results
# ---------------------------------------------------------------------------
cat("=== lagq_pairs ===\n")
print(lagq_pairs)
cat("\n=== kq_pairs ===\n")
print(kq_pairs)

cat("\n=== Lag quadratic fit (normalized Lag ~ normalized Q) ===\n")
cat(sprintf("  lag_a   = %g\n", lag_a))
cat(sprintf("  lag_b   = %g\n", lag_b))
cat(sprintf("  lag_c   = %g\n", lag_c))
cat(sprintf("  lag_d   = 0\n"))
cat(sprintf("  Lag_max = %g\n", Lag_max))
cat(sprintf("  Lag_min = %g\n", Lag_min))
cat(sprintf("  Q_max   = %g\n", Q_max))
cat(sprintf("  Q_min   = %g\n", Q_min))
cat(sprintf("  R-squared = %g\n", summary(lag_fit)$r.squared))

cat("\n=== K quadratic fit (normalized K ~ normalized Q) ===\n")
cat(sprintf("  k_a      = %g\n", k_a))
cat(sprintf("  k_b      = %g\n", k_b))
cat(sprintf("  k_c      = %g\n", k_c))
cat(sprintf("  k_d      = 0\n"))
cat(sprintf("  K_max    = %g\n", K_max))
cat(sprintf("  K_min    = %g\n", K_min))
cat(sprintf("  kq_Q_max = %g\n", kq_Q_max))
cat(sprintf("  kq_Q_min = %g\n", kq_Q_min))
cat(sprintf("  R-squared = %g\n", summary(k_fit)$r.squared))

# ---------------------------------------------------------------------------
# Save to CSV
# ---------------------------------------------------------------------------
results <- data.frame(
  lag_a    = lag_a,
  lag_b    = lag_b,
  lag_c    = lag_c,
  lag_d    = 0,
  Lag_max  = Lag_max,
  Lag_min  = Lag_min,
  Q_max    = Q_max,
  Q_min    = Q_min,
  k_a      = k_a,
  k_b      = k_b,
  k_c      = k_c,
  k_d      = 0,
  K_max    = K_max,
  K_min    = K_min,
  kq_Q_max = kq_Q_max,
  kq_Q_min = kq_Q_min,
  row.names = NULL
)

# Add lagq_pairs and kq_pairs as numbered columns
n_lagq <- nrow(lagq_pairs)
n_kq   <- nrow(kq_pairs)
lagq_lags_cols <- setNames(as.list(lagq_pairs$lagq_pairs_lags), paste0("lagq_pairs_lags_", seq_len(n_lagq)))
lagq_qs_cols   <- setNames(as.list(lagq_pairs$lagq_pairs_qs),   paste0("lagq_pairs_qs_",   seq_len(n_lagq)))
kq_ks_cols     <- setNames(as.list(kq_pairs$kq_pairs_ks),       paste0("kq_pairs_ks_",     seq_len(n_kq)))
kq_qs_cols     <- setNames(as.list(kq_pairs$kq_pairs_qs),       paste0("kq_pairs_qs_",     seq_len(n_kq)))
results <- cbind(results, as.data.frame(lagq_lags_cols),
                           as.data.frame(lagq_qs_cols),
                           as.data.frame(kq_ks_cols),
                           as.data.frame(kq_qs_cols))

outfile <- sub("\\.xml$", "_lagk_params.csv", parpixmlfile, ignore.case = TRUE)
if (outfile == parpixmlfile) {
  outfile <- paste0(parpixmlfile, "_lagk_params.csv")
}

write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nParameters saved to: %s\n", outfile))
