#!/usr/bin/env Rscript

# Creates a CV validation forcing file for event-based calibration.
# The output contains forcing timesteps that are OUTSIDE the identified
# streamflow events, so the optimizer only calibrates against event periods.
#
# Usage:
#   Rscript create_cvfold_for_calibration.R <flow_instantaneous.csv> <forcing_por_LID-zone.xml> <fold_num>
#
# Arguments:
#   flow_instantaneous.csv  : observed instantaneous flow CSV
#   forcing_por_LID-zone.xml: period-of-record forcing PIXML file
#   fold_num                : integer fold number for the output filename

box::use(
  data.table[fread, fwrite, as.data.table],
  stringr[str_match],
  tools[file_path_sans_ext],
  hydroEvents[baseflowB, eventPOT]
)

box::use(
  ./get_fews_forcing[get_fews_forcing]
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Expected 3 arguments:\n",
    "  <flow_instantaneous.csv> <forcing_por_LID-zone.xml> <fold_num>"
  )
}

flow_file  <- args[1]
por_file   <- args[2]
fold_num   <- as.integer(args[3])

# Parse LID and zone from forcing filename
bn <- file_path_sans_ext(basename(por_file))
m  <- str_match(bn, "^forcing_por_(.+)-([0-9]+)$")
if (is.na(m[1, 1])) {
  stop("Cannot parse LID and zone from filename: ", por_file,
       "\nExpected pattern: forcing_por_<LID>-<zone>.xml")
}
lid      <- m[1, 2]
zone_num <- m[1, 3]
out_dir  <- dirname(normalizePath(por_file))

# Read observed instantaneous flow
cat("Reading flow:", flow_file, "\n")
obs <- fread(flow_file)
obs[, date_ := as.Date(ISOdate(year, month, day))]

# Read forcing
cat("Reading forcing:", por_file, "\n")
forcing <- get_fews_forcing(por_file)
forcing <- as.data.table(forcing)[, .(year, month, day, hour, map_mm, mpe_mm)]
forcing[, date_ := as.Date(ISOdate(year, month, day))]

# Extract flow vector; replace NA with 0 so baseflow separation doesn't fail
flow_vec <- obs$flow_cfs
flow_vec[is.na(flow_vec)] <- 0

# Identify events using baseflow separation + peak-over-threshold
cat("Identifying events...\n")
bf     <- baseflowB(flow_vec, alpha = 0.925)
qf     <- flow_vec - bf$bf
events <- eventPOT(data=qf, threshold=80, min.diff=8 )

cat("Found", nrow(events), "events\n")

if (nrow(events) > 0) {
  event_datetimes <- data.frame(
    start = ISOdatetime(obs$year[events$srt], obs$month[events$srt],
                        obs$day[events$srt],  obs$hour[events$srt], 0, 0, tz = "UTC"),
    end   = ISOdatetime(obs$year[events$end], obs$month[events$end],
                        obs$day[events$end],  obs$hour[events$end], 0, 0, tz = "UTC")
  )
  print(event_datetimes)
}

if (nrow(events) == 0) {
  warning("No events identified; output will contain all forcing timesteps.")
  in_event <- rep(FALSE, nrow(forcing))
} else {
  # Map event row indices back to dates in the obs data.table
  event_ranges <- data.frame(
    start_date = obs$date_[events$srt],
    end_date   = obs$date_[events$end]
  )

  # Mark forcing timesteps that fall inside any event period
  in_event <- rep(FALSE, nrow(forcing))
  for (i in seq_len(nrow(event_ranges))) {
    in_event <- in_event |
      (forcing$date_ >= event_ranges$start_date[i] &
       forcing$date_ <= event_ranges$end_date[i])
  }
}

# Keep only forcing rows OUTSIDE events (these will be masked in calibration)
out_forcing <- forcing[!in_event, .(year, month, day, hour, map_mm, mpe_mm)]

out_file <- file.path(out_dir,
  sprintf("forcing_validation_cv_%d_%s-%s.csv", fold_num, lid, zone_num))

fwrite(out_forcing, out_file)
cat(sprintf("Written: %s (%d of %d forcing rows outside events)\n",
            out_file, nrow(out_forcing), nrow(forcing)))
