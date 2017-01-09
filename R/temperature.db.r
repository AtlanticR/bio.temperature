
temperature.db = function ( ip=NULL, year=NULL, p, DS, varnames=NULL, yr=NULL, dyear=NULL, ret="NULL" ) {
  
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


  if ( DS %in% c("predictions", "predictions.redo" ) ) {
    # NOTE: the primary interpolated data were already created by lbm. 
    # This routine points to this data and also creates 
    # subsets of the data where required, determined by "subregions" 
 
    outdir = file.path(project.datadirectory("bio.temperature"), "lbm", p$spatial.domain)
   
    if (DS %in% c("predictions")) {
      P = V = NULL
      fn = file.path( outdir, paste("lbm.prediction", ret,  yr, "rdata", sep=".") )
      if (is.null(ret)) ret="mean"
      if (file.exists(fn) ) load(fn) 
      if (ret=="mean") return (P)
      if (ret=="sd") return( V)
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if ( is.null(ip) ) ip = 1:p$nruns

    # downscale and warp from p(0) -> p1

    for ( r in ip ) {
      print (r)
      yr = p$runs[r, "yrs"]
      # default domain
      PP0 = lbm_db( p=p, DS="lbm.prediction", yr=yr, ret="mean")
      VV0 = lbm_db( p=p, DS="lbm.prediction", yr=yr, ret="sd")
      p0 = spatial_parameters( p=p, type=p$spatial.domain ) # from
      L0 = bathymetry.db( p=p0, DS="baseline" )
      L0i = array_map( "xy->2", L0[, c("plon", "plat")], 
        corner=c(p0$plons[1], p0$plats[1]), res=c(p0$pres, p0$pres) )
       
      sreg = setdiff( p$subregions, p$spatial.domain ) 

      for ( gr in sreg ) {
        p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], 
          corner=c(p1$plons[1], p1$plats[1]), res=c(p1$pres, p1$pres) )

        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )

        P = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
        V = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
        for (iw in 1:p$nw) {
          P[,iw] = spatial_warp( PP0[,iw], L0, L1, p0, p1, L0i, L1i )
          V[,iw] = spatial_warp( VV0[,iw], L0, L1, p0, p1, L0i, L1i )
        }

        outdir_p1 = file.path(project.datadirectory("bio.temperature"), "lbm", p1$spatial.domain)
        dir.create( outdir_p1, recursive=T, showWarnings=F )
        fn1_sg = file.path( outdir_p1, paste("lbm.prediction.mean",  yr, "rdata", sep=".") )
        fn2_sg = file.path( outdir_p1, paste("lbm.prediction.sd",  yr, "rdata", sep=".") )
        save( P, file=fn1_sg, compress=T )
        save( V, file=fn2_sg, compress=T )
        print (fn1_sg)
      }
    } 
    return ("Completed")

    if (0) {
      aoi = which( PS$z > 5 & PS$z < 3000 & PS$z.range < 500)
      levelplot( log(z) ~ plon + plat, PS[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( log(t.ar_1) ~ plon + plat, PS[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      
      levelplot( log(t.range) ~ plon + plat, PS[ aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( Z.rangeSD ~ plon + plat, PS[aoi,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }
  
  }


  # -----------------


  if (DS %in% c(  "bottom.statistics.annual", "bottom.statistics.annual.redo" )){

		tstatdir = file.path( project.datadirectory("bio.temperature"),  "stats", p$spatial.domain )
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
      P = temperature.db( p=p, DS="predictions", yr=y, ret="mean"  )
#      P = filter.bathymetry( DS=domain, Z=P )

   	  # ibaddata = which( !is.finite(P) )
			# P[ ibaddata ] = mean(P, na.rm=T )

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



  # ----------------------------


  if (DS %in% c("complete", "complete.redo" )) {
    ### a conveniance data table to reduce number of merges occuring during modelling steps
    ### annual stats and climatology are merged together
    ### essentially the base level data set for indicators.db but needed at a lower level as it is used for the other indicators

    #// temperature( p, DS="lbm.finalize(.redo)" return/create the
    #//   lbm interpolated method formatted and finalised for production use
    # NOTE: the primary interpolated data were already created by lbm. This routine points to this data and also creates 
    # subsets of the data where required, determined by "subregions" 

    if (DS=="complete") {
      TM = NULL
      outdir =  file.path( project.datadirectory("bio.temperature"), "lbm", p$spatial.domain )
      outfile =  file.path( outdir, paste( "temperature", "complete", p$spatial.domain, "rdata", sep= ".") )
      if ( file.exists( outfile ) ) load( outfile )
      return(TM)
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if (is.null(ip)) ip = 1:p$nruns

    print ( "Completing and downscaling data where necessary ..." )


    # default domain 
    p0 = p # rename for clarity
    L0 = bathymetry.db( p=p0, DS="baseline" )
    L0i = array_map( "xy->2", L0[, c("plon", "plat")], 
      corner=c(p0$plons[1], p0$plats[1]), res=c(p0$pres, p0$pres) )
 
    # merge statistics from lbm
    BS = lbm_db( p=p, DS="stats.to.prediction.grid" )
    colnames(BS) = paste("t", colnames(BS), sep=".")
    TM = cbind( TM0, BS )
  
    # merge climatology
    B = matrix( NA, nrow=nrow(TM), ncol=length(p$tyears.climatology ) )
    climatology = matrix ( NA, nrow=nrow(TM), ncol=length(p$bstats)  )

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
    TM = cbind( TM, climatology )

    TM$plon = TM$plat = TM$z = NULL
    names(TM)[ which(names(TM)=="tmean") ] = "tmean.cl"
    names(TM)[ which(names(TM)=="tamplitude") ] = "tamp.cl"
    names(TM)[ which(names(TM)=="wmin") ] = "wmin.cl"
    names(TM)[ which(names(TM)=="thalfperiod") ] = "thp.cl"
    names(TM)[ which(names(TM)=="tsd") ] = "tsd.cl"


    # bring in last stats
 

    for (iy in ip) {
      yr = p$runs[iy, "yrs"]
      print (paste( yr))
           # default domain annual stats
      E = temperature.db( DS="bottom.statistics.annual", p=p0, yr=yr  )
      E$z = NULL
      if (is.null(E)) print( paste( "bottom.statistics.annual not found for:" , yr ) )
      names(E)[ which(names(E)=="tamplitude") ] = "tamp"  # fix this at the level of "bottom statistics"
      names(E)[ which(names(E)=="thalfperiod") ] = "thp"
      TM0y = merge( TM0, E,  by =c("plon", "plat"), all.x=T, all.y=F, sort=F)
      TM0y = TM0y[ order( TM0$id), ]

      sreg = setdiff( p$subregions, p$spatial.domain ) 

      for ( gr in sreg ) {
        p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], 
          corner=c(p1$plons[1], p1$plats[1]), res=c(p1$pres, p1$pres) )

        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )

        vn = setdiff( names(TM0y), c("plon", "plat", "z" , "yr" ) )
        
        for ( ww in vn ) {
          L1[,ww] = spatial_warp( TM0y[,ww], L0, L1, p0, p1, L0i, L1i )
        }
        # return to coordinate system of original projection
        L1$plon = L1$plon0
        L1$plat = L1$plat0
        TM = L1[ , names(TM0y) ]

      }

      TM$id = NULL
      outdir_p1 = file.path(project.datadirectory("bio.temperature"), "lbm", p1$spatial.domain)
      dir.create( outdir_p1, recursive=T, showWarnings=F )
      outfile_p1 =  file.path( outdir_p1, paste( "temperature", "complete", p1$spatial.domain, "rdata", sep= ".") )
      save( TM, file=outfile_p1, compress=T )
      print( outfile_p1 )

    }

    return( outdir )
  }



}



