 
temperature.timeseries.interpolate.inla = function(p, bb, pp, zz ) {

  # clean data first
  require(INLA)

  dm = max(p$dist.km)  # km ... ie 40 km X  40 km  = 160 km^2 total ! 

  drange = c(-1,1) * dm
  plon0 = pp$plon + drange
  plat0 = pp$plat + drange
  i = which( bb$plon > plon0[1] & bb$plon < plon0[2] & bb$plat > plat0[1] & bb$plat < plat0[2] )
  
  if (length(i) < p$nMin.tbot ) return (zz) # zz is already filled with missing data  

  # only attempt interpolation if we have enough data (nMin.tbot)
  x = bb[i,] # faster to reduce the size of bb here
  # remove potentially noisy/erroneous data --- they are highly influential when there is little data 
  x$w = 1 / (( pp$plon - x$plon)**2 + (pp$plat - x$plat)**2 )# weight data in space: inverse distance squared
  x$w[ which( x$w < 1e-3 ) ] = 1e-3
  x$w[ which( x$w > 1 ) ] = 1

  xnames = c( "yr", "t", "dyear", "w") 
  x = x[ , xnames ]
  x$dataid = 1:length(i)
  x$predid = NA

  zz$w = 1
  zz$dataid = NA
  zz$predid = 1:nrow(zz)
  zz$t = NA

  x = rbind( x, zz[, colnames(x)] )  # combine predictor and estimator dframes .. this is how inla/Bayesian methods operate
  
  preds = which( is.finite( x$predid ))

  # data transformations and creation of new variables where required for raw data 
  nn = abs( diff( x$tiyr ) )
  dd = median( nn[nn>0], na.rm=TRUE )
  x$tiyr = jitter( x$tiyr, amount=dd / 20 ) # add noise as inla seems unhappy with duplicates in x?

  x$pryr = x$dyear
  x$b0 = 1

  x$cos.w  = cos( x$tiyr )
  x$sin.w  = sin( x$tiyr )

  mf = formula( t ~ 0 + b0 
               + f( tiyr, model='ar1',  diagonal=0.02 )
               + f( pryr, model='ar1', cyclic=TRUE, diagonal=0.02 ) 
  )

  tsmodel = try( inla( mf, family='gaussian', data=x, 
#                     control.inla=list(h=0.01) , 
#                     control.compute=list(dic=FALSE), 
                      control.predictor=list( compute=TRUE)), silent=TRUE )

  zz$fit = NA
  zz$fit.se = NA
  if ( ! "try-error" %in% class(tsmodel) ) {
   zz$fit = tsmodel$summary.fitted.values$mean [preds] 
   zz$fit.sd = tsmodel$summary.fitted.values$sd [preds] 
  }

  if (FALSE) {
    summary(tsmodel)
    plot( t~ tiyr, x, pch= 20)
    lines( fit~ tiyr, zz )
  }

  return( zz )

}


