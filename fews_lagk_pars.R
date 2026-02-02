#!/usr/bin/env Rscript

#.libPaths(c('/opt/RLibs/nwrfc', .libPaths()))
.libPaths(c('/data/zcui/bin/R/library', .libPaths()))

library(xml2)
library(XML)
library(data.table)


get_lagk_params <- function( parpixmlfile ){

  pixml <- read_xml(parpixmlfile)

  #xml_structure(pixml)

  ns <- xml_ns(pixml)  # get namespaces

  params <- xml_find_all( pixml, "/d1:parameters/d1:group/d1:parameter", ns )

  lagq_pairs_lags <- c()
  lagq_pairs_qs <- c()
  kq_pairs_ks <- c()
  kq_pairs_qs <- c()
  for ( p in params){
    id <- xml_attr(p, "id")
    if ( id == "LAGQ_PAIRS" ){
      rows <- xml_find_all( p, "./d1:table/d1:row", ns )
      for ( r in rows ){
        lag = as.numeric( xml_attr(r, "A") )
        lagq_pairs_lags <- c( lagq_pairs_lags, lag )
        q = as.numeric( xml_attr(r, "B") )
        lagq_pairs_qs <- c( lagq_pairs_qs, q )
      }
      lagq_pairs_df <- data.frame( lagq_pairs_lags, lagq_pairs_qs)
    }
    else if ( id == "KQ_PAIRS"  ){
      rows <- xml_find_all( p, "./d1:table/d1:row", ns )
      for ( r in rows ){
        k = as.numeric( xml_attr(r, "A") )
        kq_pairs_ks <- c( kq_pairs_ks, k )
        q = as.numeric( xml_attr(r, "B") )
        kq_pairs_qs <- c( kq_pairs_qs, q )
      }
      kq_pairs_df <- data.frame( kq_pairs_ks, kq_pairs_qs)
    }
    else if ( id == "NUMBER_OF_KQ_PAIRS" ){
      number_of_kq_pairs <- 
          as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "INFLOW_TS_DATA_TYPE" ){
       inflow_ts_data_type <-
         xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "OUTFLOW_TS_DATA_TYPE" ){
       outflow_ts_data_type <-
         xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "TRANSMISSION_LOSS_THRESHOLD_FLOW" ){
       transmission_loss_threshold_flow <-
          as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
    else if ( id == "METR_OR_ENGL_UNITS" ){
       metr_or_engl_units <-
         xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "INFLOW_TS_ID" ){
       inflow_ts_id <-
         xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "OUTFLOW_TS_INTERVAL" ){
       outflow_ts_interval <-
          as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "INFLOW_TS_INTERVAL" ){
       inflow_ts_interval <-
          as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "TRANSMISSION_LOSS_COEFFICIENT" ){
       transmission_loss_coefficient <-
          as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
    else if ( id == "OUTFLOW_TS_ID" ){
       outflow_ts_id <-
         xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "NUMBER_OF_LAGQ_PAIRS" ){
      number_of_lagq_pairs <- 
          as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "CONSTANT_K_VALUE" ){
      constant_k_value <-
          as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
    else if ( id == "CONSTANT_LAG_VALUE" ){
      constant_lag_value <-
          as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
  }
  return( list( number_of_kq_pairs=number_of_kq_pairs,
#	        inflow_ts_data_type=inflow_ts_data_type,
#	        outflow_ts_data_type=outflow_ts_data_type,
	    transmission_loss_threshold_flow=transmission_loss_threshold_flow,
#	    metr_or_engl_units=metr_or_engl_units,
#	    inlfow_ts_id=inflow_ts_id,
	    outflow_ts_interval=outflow_ts_interval,
	    inflow_ts_interval=inflow_ts_interval,
	    transmission_loss_coefficient=transmission_loss_coefficient,
	    #outflow_ts_id=outflow_ts_id,
            number_of_lagq_pairs=number_of_lagq_pairs,
	    constant_k_value=constant_k_value,
	    constant_lag_value=constant_lag_value,
	       lagq_pairs=lagq_pairs_df,
	       kq_pairs=kq_pairs_df )
        )
}

get_all_lagk_pars <- function( dir, basin) {
   basin_dir <- file.path(dir, basin)
   pifiles <- list.files( path=basin_dir, 
			 pattern=paste0("LAGK_", basin, "_.+\\.xml$"),
			 full.names=TRUE )

   return( lapply(pifiles, get_lagk_params ) )
}

add_lagk_pars_to_default_pars <- function( lagk_pars, default_pars) {

   ids <- names( lagk_pars )
   for ( id in ids ){
      lp <- lagk_pars[[id]] 
      selected_pars <- lp[!names(lp) %in% c("lagq_pairs", "kq_pairs")]

      dt <- data.table( name = names( selected_pars ),
		        type = rep( "lagk", times = length(selected_pars)),
			zone = rep(id, times = length(selected_pars)),
			value = unlist( selected_pars, use.names=FALSE) )

      default_pars <- rbindlist( list( default_pars, dt ), use.names=TRUE )

      lagq_pairs_qs <- paste0("lagq_pairs_qs_", 1:nrow(lp[["lagq_pairs"]]))
      lagq_pairs_lags <- paste0("lagq_pairs_lags_", 1:nrow(lp[["lagq_pairs"]]))
      kq_pairs_qs <- paste0("kq_pairs_qs_", 1:nrow(lp[["kq_pairs"]]))
      kq_pairs_ks <- paste0("kq_pairs_ks_", 1:nrow(lp[["kq_pairs"]]))
      
      dt <- data.table( name = lagq_pairs_qs,
		        type = rep( "lagk", times = length(lagq_pairs_qs)),
			zone = rep(id, times = length(lagq_pairs_qs)),
			value = lp[["lagq_pairs"]]$lagq_pairs_qs )

      default_pars <- rbindlist( list( default_pars, dt ), use.names=TRUE )

      dt <- data.table( name = lagq_pairs_lags,
		        type = rep( "lagk", times = length(lagq_pairs_lags)),
			zone = rep(id, times = length(lagq_pairs_lags)),
			value = lp[["lagq_pairs"]]$lagq_pairs_lags )

      default_pars <- rbindlist( list( default_pars, dt ), use.names=TRUE )

      dt <- data.table( name = kq_pairs_qs,
		        type = rep( "lagk", times = length(kq_pairs_qs)),
			zone = rep(id, times = length(kq_pairs_qs)),
			value = lp[["kq_pairs"]]$kq_pairs_qs )

      default_pars <- rbindlist( list( default_pars, dt ), use.names=TRUE )

      dt <- data.table( name = kq_pairs_ks,
		        type = rep( "lagk", times = length(kq_pairs_ks)),
			zone = rep(id, times = length(kq_pairs_ks)),
			value = lp[["kq_pairs"]]$kq_pairs_ks )

      default_pars <- rbindlist( list( default_pars, dt ), use.names=TRUE )
   } 
   return( default_pars )
}


#parfile <- "runs/1zone/GNGT2/LAGK_GNGT2_GEOT2_UpdateStates.xml"
#
#lagk_pars <- get_lagk_params(parfile)
#print(lagk_pars )
#
#all_lagk_pars <- get_all_lagk_pars( "runs/1zone", "GNGT2" )
#print( all_lagk_pars )
