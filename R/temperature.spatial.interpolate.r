
  temperature.spatial.interpolate = function( method="kernel.density", p, z=NULL, fill.missing=TRUE ) {
    #\\ interpolate spatial data and optionally fill in missing data
    #\\ it is "simple" in that there are no considerations for covariates

    out = NULL

    tofill = which( ! is.finite( z) )
    if (length( tofill) == 0 ) return(NULL)
    
    if (method=="kernel.density" ) {
      #\\ method="kernel.density" : gaussian kernel density method using fields, input has same dimensions as output
      require(fields)
      Z = matrix( NA, nrow=p$nplons, ncol=p$nplats)
      Z[p$O2M] = z
      kd = try( fields::image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z )
      z[tofill] = kd[p$O2M][ tofill]
 
      # update and repeat
      tofill = which( ! is.finite( z) )
      if (length( tofill) > 0 ) {
        Z[p$O2M] = z # update matrix
        kd = try( fields::image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z )
        z[tofill] = kd[p$O2M][ tofill]
      }
      return(z)
    }

    if (method=="gam" ) {
      #\\ method="kernel.density" : gaussian kernel density method using fields, input has same dimensions as output
      require(mgcv)
      mod = try( gam( z ~ s(plon, plat) + s(plon) +s(plat), data=p$O) )
      if ( ! (class(mod) %in% "try-error") ) {
        preds = try( predict( mod, newdata=p$O[tofill,]  ))
        if ( ! (class(preds) %in% "try-error") ) z[tofill] = preds
      }
      return (z)
    }


    if (method=="inverse.distance" ) {
      dta = cbind( p$O, z)
      names(dta) = c("plon", "plat", "z")
      gs = try( gstat( id="z", formula=z~1, locations=~plon+plat, data=dta[-tofill,], maxdist=p$maxdist, set=list(idp=.5)), silent=TRUE ) 
      if ( ! "try-error" %in% class(gs) ) {
        preds = try( predict( object=gs, newdata=dta[tofill,]  ) )
        if ( ! (class(preds) %in% "try-error") ) z[tofill] = preds[,3]
      }
      return(z)
    }
  
  }

