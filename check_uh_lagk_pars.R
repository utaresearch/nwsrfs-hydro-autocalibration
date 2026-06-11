
box::use(
  data.table[is.data.table]
)

  #' Match rows in a data.table using one of two pattern lists.
  #'
  #' @param dt        A data.table
  #' @param pattern1  A named list of column = value pairs for the first pattern
  #' @param pattern2  A named list of column = value pairs for the second pattern
  #' @return          The matched rows as a data.table
  match_rows_exclusive <- function(dt, pattern1, pattern2) {
    stopifnot(is.data.table(dt))

     apply_pattern <- function(dt, pattern) {
        idx <- rep(TRUE, nrow(dt))
        for (col in names(pattern)) {
          val <- pattern[[col]]
          if (is.function(val)) {
            idx <- idx & val(dt[[col]])
          } else {
            idx <- idx & (dt[[col]] %in% val)
          }
        }
        dt[idx]
      }

    rows1 <- apply_pattern(dt, pattern1)
    rows2 <- apply_pattern(dt, pattern2)

    has1 <- nrow(rows1) > 0
    has2 <- nrow(rows2) > 0

    if (has1 && has2) {
      stop("Ambiguous match: rows found for both pattern1 and pattern2.")
    }
    if (!has1 && !has2) {
      stop("No match: no rows found for either pattern1 or pattern2.")
    }

    if (has1) "has1" else "has2"
  }


calibrate_uh <- function( default_pars ){
 
 uh_type <-  match_rows_exclusive(
    default_pars,
    pattern1 = list(type = "uh", name = function(x) 
		    grepl("^(duration|interval|baseflow|unit_ord)", x)),
    pattern2 = list(type = "uh", name = function(x)
		    grepl("^(unit_shape|unit_toc|unit_toc_adj|unit_scale)", x)
		    )
  )

  if ( uh_type == "has1" ) FALSE else TRUE
}

calibrate_lagk <- function( default_pars ){
 
 lagk_type <-  match_rows_exclusive(
    default_pars,
    pattern1 = list(type = "lagk", name = function(x) 
		    grepl("^(init_co_pars|init_co_q_|init_co_lag_)", x)),
    pattern2 = list(type = "lagk", name = function(x)
		    grepl(paste0("^(lagtbl_a|lagtbl_b|lagtbl_c|lagtbl_d|",
				 "ktbl_a|ktbl_b|ktbl_c|ktbl_d|",
				 "lagk_lagmax|lagk_kmax|lagk_qmax|",
				 "lagk_lagmin|lagk_kmin|lagk_qmin)"), x)
		    )
  )

  if ( lagk_type == "has1" ) FALSE else TRUE
}


check_uh_pars <- function(default_pars) {
  do_calibrate <- tryCatch(
    calibrate_uh(default_pars),
    error = function(e) {
      message("Error: Could not determine UH calibration mode: ", conditionMessage(e))
      quit(status = 1)
    }
  )

  uh_names <- default_pars[default_pars[["type"]] == "uh", ][["name"]]

  if (do_calibrate) {
    has <- function(r) any(grepl(paste0("^", r), uh_names))
    has_scale   <- has("unit_scale")
    has_toc     <- has("unit_toc")
    has_toc_adj <- has("unit_toc_adj")

    if (!has("unit_shape")) {
      message("Error: calibrate_uh is TRUE but required parameter 'unit_shape' is missing.")
      quit(status = 1)
    }

    if (has_scale && (has_toc || has_toc_adj)) {
      message("Error: calibrate_uh is TRUE but 'unit_scale' and 'unit_toc'/'unit_toc_adj' are mutually exclusive.")
      quit(status = 1)
    }

    if (!has_scale && !has_toc) {
      message("Error: calibrate_uh is TRUE but neither 'unit_scale' nor 'unit_toc' is present; one is required.")
      quit(status = 1)
    }
  } else {
    required <- c("duration", "interval", "baseflow", "unit_ord")
    missing_pars <- required[!sapply(required, function(r) any(grepl(paste0("^", r), uh_names)))]
    if (length(missing_pars) > 0) {
      message("Error: The following required UH parameters are missing (calibrate_uh = FALSE):")
      message(paste(" -", missing_pars, collapse = "\n"))
      quit(status = 1)
    }
  }

  invisible(NULL)
}


check_lagk_pars <- function(default_pars) {
  do_calibrate <- tryCatch(
    calibrate_lagk(default_pars),
    error = function(e) {
      message("Error: Could not determine lag-k calibration mode: ", conditionMessage(e))
      quit(status = 1)
    }
  )

  lagk_names <- default_pars[default_pars[["type"]] == "lagk", ][["name"]]

  if (do_calibrate) {
    required <- c("lagtbl_a", "lagtbl_b", "lagtbl_c", "lagtbl_d",
                  "ktbl_a",   "ktbl_b",   "ktbl_c",   "ktbl_d",
                  "lagk_lagmax", "lagk_kmax", "lagk_qmax",
                  "lagk_lagmin", "lagk_kmin", "lagk_qmin")
  } else {
    required <- c("init_co_pairs", "init_co_q_", "init_co_lag_")
  }

  missing_pars <- required[!sapply(required, function(r) any(grepl(paste0("^", r), lagk_names)))]

  if (length(missing_pars) > 0) {
    message("Error: The following required lag-k parameters are missing (calibrate_lagk = ", do_calibrate, "):")
    message(paste(" -", missing_pars, collapse = "\n"))
    quit(status = 1)
  }

  invisible(NULL)
}
