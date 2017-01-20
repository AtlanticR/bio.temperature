
temperature.db = function ( ip=NULL, year=NULL, p, DS, varnames=NULL, yr=NULL, ret="NULL", voi="t" ) {

  # over-ride default dependent variable name if it exists
  if (exists("variables",p)) if(exists("Y", p$variables)) voi=p$variables$Y
  
  if (DS=="lbm.inputs") {

    B = hydro.db( p=p, DS="bottom.gridded.all"  ) # might be better add all data instead of the gridded ... gridding is just to make things faster
    B$tiyr = lubridate::decimal_date ( B$date )

    # globally remove all unrealistic data
    keep = which( B$t >= -3 & B$t <= 25 ) # hard limits
    if (length(keep) > 0 ) B = B[ keep, ]
    TR = quantile(B$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    keep = which( B$t >=  TR[1] & B$t <=  TR[2] )
    if (length(keep) > 0 ) B = B[ keep, ]
    
    keep = which( B$z >=  1 ) # ignore very shallow areas ..
    if (length(keep) > 0 ) B = B[ keep, ]
    
    # default output grid
    Bout = bathymetry.db( p, DS="baseline", varnames=p$varnames )
    coords = p$variables$LOCS
    covars = setdiff( p$varnames, p$variables$LOCS )
    if (length(covars)==1) {
      covs = list( Bout[,covars] )
      names(covs) = covars
      OUT  = list( LOCS = Bout[,coords], COV=covs ) 
    } else {
      OUT  = list( LOCS = Bout[,coords], COV=as.list( Bout[,covars] ) ) 
    }

    return (list(input=B, output=OUT))
  }

  # -----------------


  if ( DS %in% c("predictions", "predictions.redo" ) ) {
    # NOTE: the primary interpolated data were already created by lbm. 
    # This routine points to this data and also creates 
    # subsets of the data where required, determined by "spatial.domain.subareas" 
 
    outdir = file.path(project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain )
    
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
      L0i = array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
      sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain ) 

      for ( gr in sreg ) {
        p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )

        P = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
        V = matrix( NA, ncol=p$nw, nrow=nrow(L1) )
        for (iw in 1:p$nw) {
          P[,iw] = spatial_warp( PP0[,iw], L0, L1, p0, p1, "fast", L0i, L1i )
          V[,iw] = spatial_warp( VV0[,iw], L0, L1, p0, p1, "fast", L0i, L1i )
        }
        outdir_p1 = file.path(project.datadirectory("bio.temperature"), "modelled", voi, p1$spatial.domain)
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


  if (DS %in% c(  "bottom.statistics.annual", "bottom.statistics.annual.redo", "bottom.statistics.climatology" )){

		tstatdir = file.path( project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain )
    dir.create( tstatdir, showWarnings=F, recursive = TRUE )

    if (DS == "bottom.statistics.climatology" ) {
      clim = NULL
      fn.climatology = file.path( tstatdir, paste("bottom.statistics.climatology", p$spatial.domain, "rdata", sep=".") )
      if (file.exists( fn.climatology ) ) load(fn.climatology)
      return ( clim )
    }

		if (DS %in% c("bottom.statistics.annual")) {
      BS = NULL
      fn = file.path( tstatdir, paste("bottom.statistics.annual", p$spatial.domain, ret, "rdata", sep=".") )
      if (file.exists( fn) ) load(fn)
      return ( BS )
    }

    if (exists( "libs", p)) RLibrary( p$libs )
    if ( is.null(ip)) ip = 1:p$nruns
     
    p0 = p # to make it explicit

    grids = unique( c(p0$spatial.domain.subareas , p0$spatial.domain ) ) # operate upon every domain
 
    for (gr in grids ) {
      print(gr)

      p1 = spatial_parameters( type=gr ) #target projection    
      L1 = bathymetry.db(p=p1, DS="baseline")

      tstatdir_p1 = file.path( project.datadirectory("bio.temperature"), "modelled", voi, p1$spatial.domain )

      O = array( NA, dim=c( nrow(L1), p0$ny, length(p$bstats)) )

      for ( iy in 1:p0$ny ) {

        y = p0$yrs[iy]

  			print ( paste("Year:", y)  )
        
  			V = temperature.db( p=p1, DS="predictions", yr=y, ret="sd"  )
  			V[ V < 0.1 ] = 100  # shrink weighting of unreasonably small SEs
  		  V[ which( !is.finite(V)) ] = 1000 # "
        W = 1/V^2   # weights: inverse variance, normalised
        W = W / rowSums(W)
        V = NULL
        
        P = temperature.db( p=p1, DS="predictions", yr=y, ret="mean"  )
        O[,iy,1] = c(apply( P*W, 1, sum, na.rm=T))
        O[,iy,2]  = c(apply( (P-rowMeans(P))^2*W, 1, sum, na.rm=T )) # weighted seasonal mean sums of squares
        W = NULL
  			O[,iy,3] = c(apply( P, 1, quantile, probs=0.005, na.rm=TRUE ))
        O[,iy,4] = c(apply( P, 1, quantile, probs=0.995, na.rm=TRUE ))
  			P = NULL
        O[,iy,5] = O[,iy,4] - O[,iy,3]  # approximate as sinusoid can span 2 yrs .. max amplitude
  			# half-period .. also approximate as sinusoid can also span 2 yrs
  			# sin tranf required to make circular and then take difference and rescale
      }

      # save sp-time matrix for each stat .. easier to load into lbm this way   
      for ( st in 1:length(p$bstats) ){
        BS = O[,,st]
        fn = file.path( tstatdir_p1, paste("bottom.statistics.annual", p1$spatial.domain, p$bstats[st], "rdata", sep=".") )
        save( BS, file=fn, compress=T )
        BS = NULL
        gc()

      }

      # climatology
      clim = matrix( NA, nrow=nrow(L1), ncol=length(p$bstats) )
      for (si in 1:length(p$bstats)) {
        clim[,si] = rowMeans(O[,,si], na.rm=T)
      }
      fn.climatology = file.path( tstatdir_p1, paste("bottom.statistics.climatology", p1$spatial.domain, "rdata", sep=".") )
      save( clim, file=fn.climatology, compress=T )
      clim = NULL
      gc()

    }

    return ("Completed")
  }


  # -----------------


  if (DS %in% c(  "timeslice", "timeslice.redo" )){

    tslicedir = file.path( project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain )
    dir.create( tslicedir, showWarnings=F, recursive = TRUE )

    dyear_index = 1
    if (exists("dyears", p) & exists("prediction.dyear", p))  dyear_index = which.min( abs( p$prediction.dyear - p$dyears))

    if (DS %in% c("timeslice")) {
      O = NULL
      if (is.null(ret)) ret="mean"
      outfile =  file.path( tslicedir, paste("bottom.timeslice", p$prediction.dyear, ret, "rdata", sep=".") )
      if (file.exists( outfile ) ) load(outfile)
      return ( O )
    }

    p0 = p  # the originating parameters
    L0 = bathymetry.db( p=p0, DS="baseline" )
    L0i = array_map( "xy->2", L0, gridparams=p0$gridparams )

    Op = matrix(NA, ncol=p$ny, nrow=nrow(L0) )
    Ov = matrix(NA, ncol=p$ny, nrow=nrow(L0) )
    
    for ( r in 1:p$ny ) {
      print ( paste("Year:", p$tyears[r])  )
      V = temperature.db( p=p, DS="predictions", yr=p$tyears[r], ret="sd"  )
      if (!is.null(V)) Ov[,r] = V[,dyear_index]
      P = temperature.db( p=p, DS="predictions", yr=p$tyears[r], ret="mean"  )
      if (!is.null(P))  Op[,r] = P[,dyear_index]
      V = P = NULL
    }

    outfileP =  file.path( tslicedir, paste("bottom.timeslice", p$prediction.dyear, "mean", "rdata", sep=".") )
    save( Op, file=outfileP, compress=T )

    outfileV =  file.path( tslicedir, paste("bottom.timeslice", p$prediction.dyear, "sd", "rdata", sep=".") )
    save( Ov, file=outfileV, compress=T )

    # warp the other grids .. copy original as *0
    Op0 = Op
    Ov0 = Ov
  
    #using fields
    grids = setdiff( unique( p0$spatial.domain.subareas ), p0$spatial.domain )
    for (gr in grids ) {
      print(gr)
      p1 = spatial_parameters( type=gr ) #target projection
      L1 = bathymetry.db( p=p1, DS="baseline" )
      L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
      L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
      L1$plon_1 = L1$plon # store original coords
      L1$plat_1 = L1$plat
      L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
      p1$wght = fields::setup.image.smooth( 
        nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
        theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
      Op = Ov = matrix(NA, ncol=p$ny, nrow=nrow(L1) )
      for (iy in 1:p$ny) {
        Op[,iy] = spatial_warp( Op0[,iy], L0, L1, p0, p1, "fast", L0i, L1i )
        Ov[,iy] = spatial_warp( Ov0[,iy], L0, L1, p0, p1, "fast", L0i, L1i )
      }
  
      tslicedirp1 = file.path( project.datadirectory("bio.temperature"),  "modelled", voi, p1$spatial.domain )
      outfileP =  file.path( tslicedirp1, paste("bottom.timeslice", p$prediction.dyear, "mean", "rdata", sep=".") )
      O = Op
      save( O, file=outfileP, compress=T )
     
      outfileV =  file.path( tslicedirp1, paste("bottom.timeslice", p$prediction.dyear, "sd", "rdata", sep=".") )
      O = Ov
      save( O, file=outfileV, compress=T )
    }

    return ("Completed")
  }

  # ----------------------------


  if (DS %in% c(  "lbm.stats", "lbm.stats.redo" )){

    outdir = file.path(project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain )
    
    if (DS %in% c("lbm.stats")) {
      stats = NULL
      fn = file.path( outdir, paste( "lbm.statistics", "rdata", sep=".") )
      if (file.exists(fn) ) load(fn) 
      return( stats )
    }

    # downscale and warp from p(0) -> p1
    # default domain
    S0 = lbm_db( p=p, DS="stats.to.prediction.grid" )
    Snames = colnames(S0)
    p0 = spatial_parameters( p=p, type=p$spatial.domain ) # from
    L0 = bathymetry.db( p=p0, DS="baseline" )
    L0i = array_map( "xy->2", L0[, c("plon", "plat")], gridparams=p0$gridparams )
    sreg = setdiff( p$spatial.domain.subareas, p$spatial.domain ) 

    for ( gr in sreg ) {
      p1 = spatial_parameters( p=p, type=gr ) # 'warping' from p -> p1
      L1 = bathymetry.db( p=p1, DS="baseline" )
      L1i = array_map( "xy->2", L1[, c("plon", "plat")], gridparams=p1$gridparams )
      L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
      L1$plon_1 = L1$plon # store original coords
      L1$plat_1 = L1$plat
      L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
      p1$wght = fields::setup.image.smooth( nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
      stats = matrix( NA, ncol=ncol(S0), nrow=nrow(L1) )
      for ( i in 1:ncol(S0) ) {
        stats[,i] = spatial_warp( S0[,i], L0, L1, p0, p1, "fast", L0i, L1i )
      }
      colnames(stats) = Snames
      outdir_p1 = file.path(project.datadirectory("bio.temperature"), "modelled", voi, p1$spatial.domain)
      dir.create( outdir_p1, recursive=T, showWarnings=F )
      fn1_sg = file.path( outdir_p1, paste("lbm.statistics", "rdata", sep=".") )
      save( stats, file=fn1_sg, compress=T )
      print (fn1_sg)
    }
    return ("Completed")
  }



  # ----------------------------



  if (DS %in% c("complete", "complete.redo" )) {
    # static summaries
    
    if (DS=="complete") {
      TM = NULL
      outdir =  file.path( project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain )
      outfile =  file.path( outdir, paste( "temperature", "complete", p$spatial.domain, "rdata", sep= ".") )
      if ( file.exists( outfile ) ) load( outfile )
      Tnames = names(TM)
      if (is.null(varnames)) varnames=Tnames
      varnames = intersect( Tnames, varnames )
      if (length(varnames) == 0) varnames=Tnames  # no match .. send all
      TM = TM[ , varnames]
      return(TM)
    }


    grids = unique( c(p$spatial.domain.subareas , p$spatial.domain ) ) # operate upon every domain
 
    for (gr in grids ) {
      print(gr)

      p1 = spatial_parameters( type=gr ) #target projection    
      L1 = bathymetry.db(p=p1, DS="baseline")

      BS = temperature.db( p=p1, DS="lbm.stats" )
      colnames(BS) = paste("t", colnames(BS), sep=".")
      TM = cbind( L1, BS )

      CL = temperature.db( p=p1, DS="bottom.statistics.climatology" )
      colnames(CL) = paste(colnames(CL), "climatology", sep=".")
      TM = cbind( TM, CL )

      # bring in last stats
      outdir = file.path(project.datadirectory("bio.temperature"), "modelled", voi, p$spatial.domain)
      dir.create( outdir, recursive=T, showWarnings=F )
      outfile =  file.path( outdir, paste( "temperature", "complete", p1$spatial.domain, "rdata", sep= ".") )
      save( TM, file=outfile, compress=T )

      print( outfile )

    }

    return( outdir )
  }

}



