
temperature.db = function ( ip=NULL, year=NULL, p, DS, varnames=NULL, yr=NULL, dyear=NULL ) {
  
  if (DS=="lbm.inputs") {

    B = hydro.db( p=p, DS="bottom.gridded.all"  ) # might be better add all data instead of the gridded ... gridding is just to make things faster
    B$tiyr = lubridate::decimal_date ( B$date )

    # globally remove all unrealistic data
    keep = which( B$t >= -3 & B$t <= 25 ) # hard limits
    if (length(keep) > 0 ) B = B[ keep, ]
    TR = quantile(B$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    keep = which( B$t >=  TR[1] & B$t <=  TR[2] )
    if (length(keep) > 0 ) B = B[ keep, ]
    
    # default output grid
    Bout = bathymetry.db( p, DS="baseline", varnames=p$varnames )
    coords = p$variables$LOCS
    covars = setdiff( p$varname, p$variables$LOCS )
    OUT  = list( 
      LOCS = Bout[,coords],
      COV = as.list( Bout[,covars] ) 
    )          

    return (list(input=B, output=OUT))
  }

  # -----------------


  if ( DS %in% c("lbm.finalize.redo", "lbm.finalize" )) {
    #// temperature( p, DS="lbm.finalize(.redo)" return/create the
    #//   lbm interpolated method formatted and finalised for production use
    fn = file.path(  project.datadirectory("bio.temperature"), "lbm",
      paste( "temperature", "lbm", "finalized", p$spatial.domain, "rdata", sep=".") )
    if (DS =="lbm.finalize" ) {
      B = NULL
      if ( file.exists ( fn) ) load( fn)
      return( B )
    }

    # merge into statistics
    B = bathymetry.db( p=p, DS="baseline" ) # add to the input data

    BS = lbm_db( p=p, DS="stats.to.prediction.grid" )
    colnames(BS) = paste("t", colnames(BS), sep=".")
    B = cbind( B, BS )
    rm (BS); gc()
  
    save( B, file=fn, compress=TRUE)
    return(fn)

    if (0) {
      aoi = which( B$z > 5 & B$z < 3000 & B$z.range < 500)
      levelplot( log(z) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( log(t.ar_1) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      
      levelplot( log(t.range) ~ plon + plat, B[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( Z.rangeSD ~ plon + plat, B[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }
  }


  if (DS %in% c(  "lbm.prediction.mean", "lbm.prediction.sd", "lbm.prediction.redo", "lbm.prediction" )){

		# NOTE: the primary interpolated data were already created by lbm. This routine points to this data and also creates 
    # subsets of the data where required, determined by "subregions" 
		savedir = file.path(p$project.root, "lbm", p$spatial.domain ) # already saved by lbm_db .. this is a synonym

    if (DS %in% c("lbm.prediction", "lbm.prediction.mean")) {
      P = NULL
      fn1 = file.path( savedir, paste("lbm.prediction.mean",  yr, "rdata", sep=".") )
      if (file.exists( fn1) ) load(fn1)
      return ( P )
    }

    if (DS %in% c("lbm.prediction.sd")) {
      V = NULL
      fn2 = file.path( savedir, paste("lbm.prediction.sd",  yr, "rdata", sep=".") )
      V =NULL
      if (file.exists( fn2) ) load(fn2)
      return ( V )
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if ( is.null(ip) ) ip = 1:p$nruns

    #downscale
    for ( r in ip ) {
      yr = p$runs[r, "yrs"]
      # default domain
      PP0 = temperature.db( p=p, DS="lbm.prediction.mean", yr=yr)
      VV0 = temperature.db( p=p, DS="lbm.prediction.sd", yr=yr)
      p0 = spatial_parameters( p=p, type=p$default.spatial.domain ) # from
      p$wght = fields::setup.image.smooth( nrow=p0$nplons, ncol=p0$nplats, dx=p0$pres, dy=p0$pres,
              theta=p$phi, xwidth=p$nsd*p$phi, ywidth=p$nsd*p$phi )
      L0 = bathymetry.db( p=p0, DS="baseline" )
      sreg = setdiff( p$subregions, p$spatial.domain.default ) 
      for ( gr in sreg ) {
        p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
          L1 = bathymetry.db( p=p1, DS="baseline" )
          P = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
          V = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
          for (iw in 1:p$nw) {
            P[,iw] = spatial_warp ( Z0=PP0[,iw], L0, L1, p0=p, p1=p1 )
            V[,iw] = spatial_warp ( Z0=VV0[,iw], L0, L1, p0=p, p1=p1 )
          }
          savedir_sg = file.path(p$project.root, "lbm", p1$spatial.domain ) 
          dir.create( savedir_sg, recursive=T, showWarnings=F )
          fn1_sg = file.path( savedir_sg, paste("lbm.prediction.mean",  yr, "rdata", sep=".") )
          fn2_sg = file.path( savedir_sg, paste("lbm.prediction.sd",  yr, "rdata", sep=".") )
          save( P, file=fn1_sg, compress=T )
          save( V, file=fn2_sg, compress=T )
          print (fn1_sg)
      }
    } 
    return ("Completed")
  }


  # -----------------

  if (DS %in% c(  "bottom.statistics.annual", "bottom.statistics.annual.redo" )){

		tstatdir = project.datadirectory("bio.temperature",  "data", "stats", p$spatial.domain )
    dir.create( tstatdir, showWarnings=F, recursive = TRUE )

		if (DS %in% c("bottom.statistics.annual")) {
      O = NULL
      fn = file.path( tstatdir, paste("bottom.statistics.annual",  yr, "rdata", sep=".") )
      if (file.exists( fn) ) load(fn)
      return ( O )
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if ( is.null(ip)) ip = 1:p$nruns

    for ( r in ip ) {
      y = p$runs[r, "yrs"]
			print ( paste("Year:", y)  )

			O = bathymetry.db( p=p, DS="baseline" )
      P = temperature.db( p=p, DS="lbm.prediction.mean", yr=y  )
   		P[ P < -2 ] = -2
		  P[ P > 30 ] = 30
		  ibaddata = which( !is.finite(P) )
			P[ ibaddata ] = mean(P, na.rm=T )

			V = temperature.db( p=p, DS="lbm.prediction.sd", yr=y  )
			V[ V < 0.1 ] = 100  # shrink weighting of unreasonably small SEs
		  V[ which( !is.finite(V)) ] = 1000 # "
			V[ ibaddata ] = 10000 # " smaller still

      O$yr = y
      O$wmin = NA
      O$wmax = NA
      O$tmin = NA
      O$tmax = NA
      O$tsd = NA
      O$tmean = NA
      O$tamplitude = NA
      O$thalfperiod = NA

      O$wmin = apply( P, 1, which.min )
      O$wmax = apply( P, 1, which.max )
			O$tmin = apply( P, 1, quantile, probs=0.005 )
      O$tmax = apply( P, 1, quantile, probs=0.995 )

			W = 1/V^2   # weights: inverse variance, normalised
			W = W / rowSums(W)
			O$tmean = apply( P*W, 1, sum, na.rm=T)

			SS = (P-O$tmean)^2 # sums of squares
			O$tsd  = apply( SS*W, 1, sum, na.rm=T ) # weighted seasonal mean sums of squares

			O$tamplitude = O$tmax- O$tmin  # approximate as sinusoid can span 2 yrs .. max amplitude

			# half-period .. also approximate as sinusoid can also span 2 yrs
			# sin tranf required to make circular and then take difference and rescale
      O$thalfperiod = abs( sin(O$wmax/p$nw*pi) - sin(O$wmin/p$nw*pi) ) * p$nw/pi

      fn =  file.path( tstatdir, paste("bottom.statistics.annual", y, "rdata", sep=".") )
      save( O, file=fn, compress=T )

      rm (O, P) ; gc()
    }
    return ("Completed")
  }

  # -----------------

  if (DS %in% c("climatology", "climatology.redo") ) {

    # form a basic prediction surface in planar coords for SS habitat for
    # factors that do not "change" rapidly and
    # e.g. climatological means of temperature characteristics, substrate, bathymetry

    outdir = file.path( project.datadirectory("bio.temperature"), "data", "interpolated" )
    dir.create(outdir, recursive=T, showWarnings=F)

    outfile =  file.path( outdir, paste("PS.climatology", p$spatial.domain, "rdata", sep="." ) )
    if ( p$spatial.domain =="snowcrab" ) outfile = file.path( outdir, paste("PS.climatology", "ESS", "rdata", sep="." ) )

    if ( DS=="climatology" ) {
      if (file.exists(outfile)) load( outfile )
      if ( p$spatial.domain =="snowcrab" ) {
        id = bathymetry.db( DS="lookuptable.sse.snowcrab" )
        PS = PS[ id, ]
      }
      return (PS)
    }

    PS = bathymetry.db( p=p, DS="baseline" )  # SS to a depth of 500 m  the default used for all planar SS grids
    PS$id = 1:nrow(PS)
 
    TC = temperature.db( p=p, DS="lbm.finalize")
    PS = merge( PS, TC, by="?") # or cbind ... Not finished

    B = matrix( NA, nrow=nrow(PS), ncol=length(p$tyears.climatology ) )
    climatology = matrix ( NA, nrow=nrow(PS), ncol=length(p$bstats)  )

		for ( iv in 1:length(p$bstats) )  {
			vn = p$bstats[iv]
      print (vn)
      B[] = NA
      for ( iy in 1:length(p$tyears.climatology) ) {
				y = p$tyears.climatology[iy]
        H = temperature.db( p=p, DS="bottom.statistics.annual", yr=y )
				if (! is.null( H ) ) B[,iy] = H[,vn]
			}
			climatology[,iv] = rowMeans(B, na.rm=T)
		}
    colnames (climatology) = p$bstats
    PS = cbind( PS, climatology )
    save (PS, file=outfile, compress=T )
    return( outfile )
  }


  # ----------------------------


  if (DS %in% c("timeslice", "timeslice.redo") ) {

    # form a timeslice of temperatures at a given time of year (dyear)

    outdir = file.path( project.datadirectory("bio.temperature"), "data", "timeslice" )
    dir.create(outdir, recursive=T, showWarnings=F)

    outfile =  file.path( outdir, paste("PS.timeslice", p$spatial.domain, "rdata", sep="." ) )
    if ( p$spatial.domain =="snowcrab" ) outfile = file.path( outdir, paste("PS.timeslice", "ESS", "rdata", sep="." ) )

    if ( DS=="timeslice" ) {
      if (file.exists(outfile)) load( outfile )
      if ( p$spatial.domain =="snowcrab" ) {
        id = bathymetry.db( DS="lookuptable.sse.snowcrab" )
        PS = PS[ id, ]
      }
      return (PS)
    }
  
    B0 = bathymetry.db( p=p, DS="baseline" ) 
    id = which.min( abs( dyear - p$dyear_centre))
    PS = matrix( NA, ncol=p$ny, nrow=nrow(B0) )
    for (iy in 1:p$ny){
      u = NULL
      u = temperature.db(p=p, DS="lbm.prediction.mean", yr=p$tyears[iy] )
      if (!is.null(u)) PS[,iy] = u
    }

    save( PS, file=outfile, compress=TRUE )
    return(outfile)
  }

  # ----------------------------

  if (DS %in% c("complete", "complete.redo" )) {
    ### a conveniance data table to reduce number of merges occuring during modelling steps
    ### annual stats and climatology are merged together
    ### essentially the base level data set for indicators.db but needed at a lower level as it is used for the other indicators

    if (DS=="complete") {
      PS = NULL
      outdir =  file.path( project.datadirectory("bio.temperature"), "data", "interpolated", "complete", p$spatial.domain )
      outfile =  file.path( outdir, paste( "PS", year, "rdata", sep= ".") )
      if ( file.exists( outfile ) ) load( outfile )
      return(PS)
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if (is.null(ip)) ip = 1:p$nruns

    print ( "Completing and downscaling data where necessary ..." )

    # default domain climatology
    p0 = spatial_parameters( type=p$spatial.domain.default )
    Z0 = matrix( NA, nrow=p0$nplons, ncol=p0$nplats)
    PS0 = bathymetry.db ( p=p0, DS="baseline" )
    PS0$id =1:nrow(PS0)
    CL = temperature.db(p=p0, DS="climatology")
    CL$plon = CL$plat = CL$z = NULL
    names(CL)[ which(names(CL)=="tmean") ] = "tmean.cl"
    names(CL)[ which(names(CL)=="tamplitude") ] = "tamp.cl"
    names(CL)[ which(names(CL)=="wmin") ] = "wmin.cl"
    names(CL)[ which(names(CL)=="thalfperiod") ] = "thp.cl"
    names(CL)[ which(names(CL)=="tsd") ] = "tsd.cl"
    PS0 = merge( PS0, CL,  by =c("id"), all.x=T, all.y=F, sort=F )
    PS0_m = cbind( (PS0$plon-p0$plons[1])/p0$pres + 1, (PS0$plat-p0$plats[1])/p0$pres + 1) # row, col indices in matrix form
    rm (CL); gc()

    for (iy in ip) {
      yr = p$runs[iy, "yrs"]
      print (paste( yr))
           # default domain annual stats
      E = temperature.db( DS="bottom.statistics.annual", p=p0, yr=yr  )
      E$z = NULL
      if (is.null(E)) print( paste( "bottom.statistics.annual not found for:" , yr ) )
      names(E)[ which(names(E)=="tamplitude") ] = "tamp"  # fix this at the level of "bottom statistics"
      names(E)[ which(names(E)=="thalfperiod") ] = "thp"
      PS0y = merge( PS0, E,  by =c("plon", "plat"), all.x=T, all.y=F, sort=F)
      PS0y = PS0y[ order( PS0$id), ]

      if ( p$spatial.domain == p$spatial.domain.default ) {
        PS = PS0y
      } else {
        # down scale data to alternate grids
        PS = bathymetry.db ( p=p, DS="baseline" )

        PS = planar2lonlat( PS, proj.type=p$internal.projection )  # convert new locations to lon lat
        PS$yr = yr
        PS$plon0 = PS$plon
        PS$plat0 = PS$plat
        PS = lonlat2planar( PS, proj.type=p0$internal.projection )  # convert lon lat to coord system of p0
        locsout = PS[, c("plon", "plat")]
        p0$wgts = fields::setup.image.smooth( nrow=p0$nplons, ncol=p0$nplats, dx=p0$pres, dy=p0$pres,
              theta=p$phi, xwidth=p$nsd*p$phi, ywidth=p$nsd*p$phi )
        vn = setdiff( names(PS0y), c("plon", "plat", "z" , "yr" ) )
        for ( ww in vn ) {
          Z = Z0
          Z[PS0_m] = PS0y[,ww]
          # simple linear interpolations
          is = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=Z), loc=locsout )
          ii = which( is.na( is) )
          if ( length( ii)> 0 ) {
            # smoothed surface ..but fast!  mostly complex edges such as coastlines ..
            kd = try( fields::image.smooth( Z, dx=p0$pres, dy=p0$pres, wght=p0$wgts )$z  )
            if ( ! (class(kd) %in% "try-error") ) is[ii] = kd[ii]
          }
          jj = which( is.na( is) )
          if ( length( jj)> 0 ) is[jj] = median( PS0y[,ww], na.rm=TRUE )
          PS[,ww] = is
        }
        # return to coordinate system of original projection
        PS$plon = PS$plon0
        PS$plat = PS$plat0
        PS = PS[ , names(PS0y) ]
      }
      PS$id = NULL
      outdir =  file.path( project.datadirectory("bio.temperature"), "data", "interpolated", "complete", p$spatial.domain )
      outfile =  file.path( outdir, paste( "PS", yr, "rdata", sep= ".") )
      dir.create(outdir, recursive=T, showWarnings=F)
      save (PS, file=outfile, compress=T )
    }
    return( outdir )
  }



}
