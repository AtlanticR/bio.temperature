
  temperature.timeseries.interpolate = function( ip=NULL, p ) {

    if (exists( "init.files", p)) LoadFiles( p$init.files ) 
    if (exists( "libs", p)) RLibrary( p$libs ) 
    if (is.null(ip)) ip = 1:p$nruns

    # default output grid
    z0 = expand.grid( dyear=1:p$nw, yr=p$tyears )
    attr( z0, "out.attrs" ) = NULL
    z0$fit = NA  # these will be filled in with predicted fits and se's
    z0$se  = NA
    z0$tiyr = z0$yr + (z0$dyear-0.5) / p$nw # mid-points
    z0 = z0[ order(z0$tiyr), ]

    if ( p$tsmethod %in% c("annual", "seasonal.basic", "seasonal.smoothed", "harmonics.1", "harmonics.2", "harmonics.3" ) ) {
        interpolate.ts = temperature.timeseries.interpolate.gam
    }
    if (p$tsmethod %in% c("inla.ts.simple" ) ) {
        interpolate.ts = temperature.timeseries.interpolate.inla
    }
    
    if ( p$tsmethod %in% c("ar" ) ) {
        interpolate.ts = temperature.timeseries.interpolate.spectral  ## to do..
    }
    

    B = hydro.db( p=p, DS="bottom.gridded.all"  )
    B$tiyr = lubridate::decimal_date ( B$date )
 
    # globally remove all unrealistic data  
    keep = which( B$t >= -3 & B$t <= 25 ) # hard limits
    if (length(keep) > 0 ) B = B[ keep, ]
    TR = quantile(B$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    keep = which( B$t >=  TR[1] & B$t <=  TR[2] )
    if (length(keep) > 0 ) B = B[ keep, ]

    P = bathymetry.db( p=p, DS="baseline" )
    tbot <- bigmemory::attach.big.matrix( p$descriptorfile.tbot  )
    tbot.se <- bigmemory::attach.big.matrix( p$descriptorfile.tbotse  )

    for ( iip in ip ) {
      mm = p$runs[iip,"loc"]
      if ( is.nan( tbot[mm,1] )) next() #  this location is problematic .. skip
      if ( !is.na( tbot[mm,1] )) next() # has a solution from previous run .. skip
      tbot[mm,1] = NaN # flag as being operated upon .. in case a restart is needed
      res = NULL
      res = try( interpolate.ts ( p=p, bb=B, pp=P[mm,], zz=z0 ), silent=TRUE )
      if ( class(res) %in% "try-error" ) next()
      if ( any(is.finite(res$fit)) ) {
        print (mm)			
        tbot[ mm,] = res$fit
        tbot.se[mm,] = res$se
      }
    } # end each point
    
    return( "completed")
  }



