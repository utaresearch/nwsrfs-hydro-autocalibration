#!/usr/bin/env Rscript

#.libPaths(c('/opt/RLibs/nwrfc', .libPaths()))
.libPaths(c('/data/zcui/bin/R/library', .libPaths()))

library(xml2)
library(XML)
library(data.table)


get_uhg_params <- function( parpixmlfile ){

  pixml <- read_xml(parpixmlfile)

  #xml_structure(pixml)

  ns <- xml_ns(pixml)  # get namespaces

  params <- xml_find_all( pixml, "/d1:parameters/d1:group/d1:parameter", ns )

  params_names <- c()
  params_values <- c()
  uhg_ordinates <- c()
  for ( p in params){
    id <- xml_attr(p, "id")
    #value <- as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    print( id )
    if ( id == "UHG_ORDINATES" ){
      ords <- xml_find_all( p, "./d1:table/d1:row", ns )
      for ( o in ords ){
        ord = as.numeric( xml_attr(o, "A") )
        print( ord )
	uhg_ordinates <- c( uhg_ordinates, ord )
      }
    }
    else if ( id == "DRAINAGE_AREA"  ){
	    drainage_area <- 
		    as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
    else if ( id == "UHG_DURATION" ) {
	    uhg_duration <- 
		    as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "UNIT" ) {
        unit <- xml_text( xml_find_first( p, "./d1:stringValue", ns ) )
    }
    else if ( id == "UHG_INTERVAL" ){
	 uhg_interval <- 
		    as.integer( xml_text( xml_find_first( p, "./d1:intValue", ns ) ) )
    }
    else if ( id == "CONSTANT_BASE_FLOW" ){
	    baseflow <- 
		    as.numeric( xml_text( xml_find_first( p, "./d1:dblValue", ns ) ) )
    }
  }
  if ( unit == "ENGLISH" ){
     # square mile to square km
     drainage_area <- drainage_area * 2.58998811 # 1 square mile = 2.59 square kilometers
     #convert ordinates from CFS per INCH to CMS per MM
     uhg_ordinates <- uhg_ordinates / 35.3147  # 35.3147 FT3 = 1 M3
     uhg_ordinates <- uhg_ordinates / 25.4     # 1 IN = 25.4 MM
  }


  return( list( constant_base_flow = baseflow, 
                uhg_interval = uhg_interval,
                uhg_duration = uhg_duration,
                drainage_area = drainage_area,
                data.table( ordinates = uhg_ordinates ) ) )
}
