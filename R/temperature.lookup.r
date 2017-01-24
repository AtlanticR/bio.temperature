
temperature.lookup = function( p, locs, timestamp ) {

  domain = bathymetry.db(p=p, DS="baseline" )
  domain_map = lbm::array_map( "xy->1", domain, gridparams=p$gridparams )

  locs_map = lbm::array_map( "xy->1", locs, gridparams=p$gridparams )
  locs_index = match( locs_map, domain_map )
  yrs = lubridate::year(timestamp)
  yrs_index = match( yrs, p$yrs )
  dyear = lubridate::decimal_date( timestamp ) - yrs
  dyear_index = as.numeric( cut( dyear, breaks=c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] ) , include.lowest=T, ordered_result=TRUE ) )
  dindex = cbind(locs_index, yrs_index, dyear_index ) # check this
  temp = temperature.db( p=p, DS="spatial.annual.seasonal")
  out = temp[dindex]
  return(out)
}


