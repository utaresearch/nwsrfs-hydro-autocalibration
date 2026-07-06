#!/usr/bin/env Rscript

# Creates cross-validation fold forcing CSV files from a period-of-record forcing XML.
#
# Usage (4 folds auto-split):
#   Rscript create_cvfold_forcing_csv.R <por_forcing.xml>
#
# Usage (single fold from explicit date range):
#   Rscript create_cvfold_forcing_csv.R <por_forcing.xml> <start_date> <end_date> <fold_num>
#   dates in YYYY-MM-DD format
#
# Output: forcing_validation_cv_<fold>_<LID>-<zone>.csv  (same directory as input)

box::use(
  data.table[fwrite, as.data.table],
  stringr[str_match],
  tools[file_path_sans_ext]
)

box::use(
  ./get_fews_forcing[get_fews_forcing]
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
  por_file  <- args[1]
  start_dt  <- NULL
  end_dt    <- NULL
  fold_num  <- NULL
} else if (length(args) == 4) {
  por_file  <- args[1]
  start_dt  <- as.Date(args[2])
  end_dt    <- as.Date(args[3])
  fold_num  <- as.integer(args[4])
} else {
  stop(
    "Expected 1 or 4 arguments:\n",
    "  1 arg : <por_forcing.xml>\n",
    "  4 args: <por_forcing.xml> <start_date YYYY-MM-DD> <end_date YYYY-MM-DD> <fold_num>"
  )
}

# Extract LID and zone number from filename (pattern: forcing_por_<LID>-<zone>.xml)
bn    <- file_path_sans_ext(basename(por_file))
m     <- str_match(bn, "^forcing_por_(.+)-([0-9]+)$")
if (is.na(m[1, 1])) {
  stop("Cannot parse LID and zone from filename: ", por_file,
       "\nExpected pattern: forcing_por_<LID>-<zone>.xml")
}
lid      <- m[1, 2]
zone_num <- m[1, 3]
out_dir  <- dirname(normalizePath(por_file))

cat("Reading forcing:", por_file, "\n")
forcing <- get_fews_forcing(por_file)
forcing <- as.data.table(forcing)[, .(year, month, day, hour, map_mm, mpe_mm)]

write_fold <- function(fold_data, fold_n) {
  out_file <- file.path(out_dir,
    sprintf("forcing_validation_cv_%d_%s-%s.csv", fold_n, lid, zone_num))
  fwrite(fold_data, out_file)
  cat("Written:", out_file, sprintf("(%d rows)\n", nrow(fold_data)))
}

if (is.null(start_dt)) {
  # Auto-split: divide rows evenly into 4 folds
  n         <- nrow(forcing)
  fold_size <- floor(n / 4)
  for (k in 1:4) {
    idx_start <- (k - 1) * fold_size + 1
    idx_end   <- if (k == 4) n else k * fold_size
    write_fold(forcing[idx_start:idx_end], k)
  }
} else {
  # Single fold from explicit date range
  forcing[, date_ := as.Date(ISOdate(year, month, day))]
  fold_data <- forcing[date_ >= start_dt & date_ <= end_dt,
                       .(year, month, day, hour, map_mm, mpe_mm)]
  if (nrow(fold_data) == 0) {
    stop("No data found for date range ", start_dt, " to ", end_dt)
  }
  write_fold(fold_data, fold_num)
}
