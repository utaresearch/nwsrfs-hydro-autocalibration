
get_fews_forcing <- function( pixmlfile ){

  pixml <- read_xml(pixmlfile)

  #xml_structure(pixml)

  ns <- xml_ns(pixml)  # get namespaces

  timezonenode <- xml_find_all( pixml, "/d1:TimeSeries/d1:timeZone", ns )

  timezone = 0
  if ( length(timezonenode) > 0  ){
     
     timezone = as.numeric( xml_text( timezonenode ) )
  }
  print( timezone )
  series <- xml_find_all( pixml, "/d1:TimeSeries/d1:series", ns )

  all_series_dfs <- list()
  ids <- list()
  for ( s in series){
    #print( s )
    id <- xml_find_first( s, "./d1:header/d1:locationId", ns )
    id_text <- xml_text(id)          # Get text content
    par_id <- xml_find_first( s, "./d1:header/d1:parameterId", ns )
    par_id_text <- xml_text(par_id)          # Get text content
    #  cat(sprintf("Item ID: %s, Value: %s\n", id, value))
    events <- xml_find_all( s, "./d1:event", ns )

    timestamp <- c()
    values <- c()
    locationids <- c()
    parids <- c()
    ids <- c( ids, paste(id_text, par_id_text, sep="_" ) )
    if ( length (events) > 0 ){
      for ( e in events ){
        d <- xml_attr(e, "date")
        t <- xml_attr(e, "time")
        v <- xml_attr(e, "value")
	timestr = format(as.POSIXct( paste(d,t), format = "%Y-%m-%d %H:%M:%S") 
			 + abs(timezone)*60*60, "%Y%m%d%H%M%S"  ) 
        timestamp <- c(timestamp, timestr  ) 
        values <- append( values, as.numeric( v ) )
	locationids <- append( locationids, id_text )
	parids <- append( parids, par_id_text )
      }

      #lseries <- list("date_time" = timestamp )
      #lseries[[ paste(id_text, par_id_text, sep='_')]] = values
      #lseries <- data.frame(date_time=timestamp, valuename = values ) 
      lseries <- as.data.frame( cbind(timestamp, values, parids, locationids ) ) 
      #names(lseries)[names(lseries) == "values"] <- paste(id_text, par_id_text, sep='_')
    }
    all_series_dfs <- append( all_series_dfs, list(lseries) )
  }

  dt_map <- setDT( all_series_dfs[[1]] )
  dt_mpe <- setDT( all_series_dfs[[2]] )

  mapids = dt_map[, unique(locationids)]
  mappids = dt_map[, unique(parids)]

  if ( length( mapids )!= 1  || length( mappids ) != 1 ){
    stop("There is more then one location or parameter, exiting")
  } 

  mpeids = dt_mpe[, unique(locationids)]
  mpepids = dt_mpe[, unique(parids)]

  if ( length( mpeids )!= 1  || length( mpepids ) != 1 ){
    stop("There is more then one location or parameter, exiting")
  } 

  if ( ! identical( mapids, mpeids) ){
    stop("Locations are not the same, exiting")
  } 

  setnames(dt_map, old = c("values"), new = c("map_mm"))
  setnames(dt_mpe, old = c("values"), new = c("mpe_mm"))

  merged_dt <- merge(dt_map, dt_mpe, by = "timestamp", all = TRUE)

  merged_dt[is.na(mpe_mm), mpe_mm := "-999" ]
  merged_dt[, map_mm := as.numeric(map_mm)]
  merged_dt[, mpe_mm := as.numeric(mpe_mm)]
  str=( merged_dt )
  merged_dt = merged_dt[, ':=' 
   	(year=as.integer(substring(timestamp, first=1, last = 4 )), 
	     month=as.integer(substring(timestamp,first=5, last=6 )), 
	     day=as.integer(substring(timestamp, first=7, last=8 )), 
	     hour=as.integer(substring(timestamp, first=9, last=10)))]
  return( merged_dt )
}

