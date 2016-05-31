 
temperature.timeseries.interpolate.gam = function(p, bb, pp, zz ) {

  # harmonic method explanation:
  # old method use 1 harmonic ... forcing sinusoid as a seasonal component
  # to add an offset to a trig function (b) must add cos to a sin function
  # y ~ a + c*sin(x+b)
  # y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)  
  #   .. as C*sin(x+b) = C*( cos(b) * sin(x) + sin(b) * cos(x) )
  # y ~ b0 + b1*x1 + b2*x2
  # where: 
  #   a = b0
  #   c^2 = b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b))
  #   c = sqrt(b1^2 + b2^2)
  #   b1/b2 = tan(b)  
  #   b = arctan(b1/b2)

# choose model formula for GAM-based models
  mf = switch( p$tsmethod ,
    annual = ' t ~ s(yr) ',
    seasonal.basic = ' t ~ s(yr) + s(dyear, bs="cc") ', 
    seasonal.smoothed = ' t ~ s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ', 
    harmonics.1 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ', 
    harmonics.2 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' , 
    harmonics.3 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) '
  )

  mf = formula(mf)
  
  for ( dm in p$dist.km ) { 
      drange = c(-1,1) * dm
      plon0 = pp$plon + drange
      plat0 = pp$plat + drange
      i = which( bb$plon > plon0[1] & bb$plon < plon0[2] & bb$plat > plat0[1] & bb$plat < plat0[2] )
      if (length(i) > p$nMin.tbot ) {  
        
      #  browser()

        # only attempt interpolation if we have enough data (nMin.tbot)
        x = bb[i,] # faster to reduce the size of bb here
        # remove potentially noisy/erroneous data --- they are highly influential when there is little data 
        # xt = quantile( x$t, probs=c(0.005, 0.995) )
        # xi = which( x$t >= xt[1] & x$t <= xt[2] ) 
        # if (length(xi) < p$nMin.tbot ) next()
        # x = x[xi, ] 
        
        x$w = 1 / (( pp$plon - x$plon)**2 + (pp$plat - x$plat)**2 )# weight data in space: inverse distance squared
        x$w[ which( x$w < 1e-3 ) ] = 1e-3
        x$w[ which( x$w > 1 ) ] = 1

        # data transformations and creation of new variables where required for raw data 
        if ( p$tsmethod %in% c( "harmonics.1", "harmonics.2", "harmonics.3"  ) ) {
          x$cos.w  = cos( 2*pi*x$tiyr )
          x$sin.w  = sin( 2*pi*x$tiyr )
         
          years.with.data = unique( x$yr)
          no.years = which( !( zz$yr %in% years.with.data) )
          zz$yr[ no.years ] = median( years.with.data) 
          zz$cos.w  = cos( zz$tiyr )
          zz$sin.w  = sin( zz$tiyr )
          
          # compute additional harmonics only if required (to try to speed things up a bit)
          if ( p$tsmethod %in% c( "harmonics.2", "harmonics.3"  ) ) {
            x$cos.w2 = cos( 2*x$tiyr )
            x$sin.w2 = sin( 2*x$tiyr )
            zz$cos.w2 = cos( 2*zz$tiyr )
            zz$sin.w2 = sin( 2*zz$tiyr )
          }
          if ( p$tsmethod %in% c( "harmonics.3"  ) ) {
            x$cos.w3 = cos( 3*x$tiyr )
            x$sin.w3 = sin( 3*x$tiyr )
            zz$cos.w3 = cos( 3*zz$tiyr )
            zz$sin.w3 = sin( 3*zz$tiyr )
          }
        }

        # estimate model parameters
        tsmodel = NULL 
        tsmodel = switch( p$gam.optimizer ,
          bam = try( bam( mf, data=x, weights=w ) ) ,
          bfgs = try( gam( mf, data=x, weights=w, optimizer=c("outer","bfgs")  ) ) ,
          perf = try( gam( mf, data=x, weights=w, optimizer=c("perf")  ) ) ,
          newton = try( gam( mf, data=x, weights=w, optimizer=c("outer","newton")  ) ) ,
          nlm = try( gam( mf, data=x, weights=w, optimizer=c("outer","nlm")  ) ) 
        )
        
        if ( ! "try-error" %in% class(tsmodel) ) {
          out = try( predict( tsmodel, newdata=zz, type="response", se.fit=T ) ) 
          if ( ! "try-error" %in% class( out ) ) {
            zz$fit = out$fit
            zz$se = out$se.fit 
            break()  # candidate predictions found exit inner loop (dm)
          }
        }
      }
    } # end for dm loop						
 
    return(zz)
}


