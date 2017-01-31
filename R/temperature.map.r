
temperature.map = function( ip=NULL, p=NULL, type="all", vname=NULL ) {

  # ip is the first parameter passed in the parallel mode
  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) ip = 1:p$nruns

  loc = bathymetry.db(p=p, DS="baseline" )

	require( lattice )


  if ( type=="all" ) {

    allgrids = unique(c( p$spatial.domain.subareas, p$spatial.domain) )
    for ( gr in allgrids ) {
      print (gr)
      p1 = spatial_parameters(  p=p, type= gr )
      p1 = make.list( list( yrs=p1$yrs), Y=p1 )
      temperature.map( p=p1, type="lbm.stats" ) # no parallel option .. just a few
      temperature.map( p=p1, type="climatology" ) # no parallel option .. just a few
      parallel.run( temperature.map, p=p1, type="seasonal" ) # all seasonal predicted means
      for ( bs in p$bstats ) {
        parallel.run( temperature.map, p=p1, type="annual", vname=bs ) # bottom.statistics.annual
      }
    }

  }

  # -----------------------


  if ( type=="seasonal" ) {

    bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain, "bottom.predictions", "seasonal" )
    dir.create( bottomdir.maps, recursive=T, showWarnings=F )
    datarange = seq(-1, 12, length.out=100)
    cols = color.code( "blue.black", datarange )
    for (iy in ip ) {
      y = p$runs[iy, "yrs"]
      H = temperature.db( p=p, DS="predictions", yr=y, ret="mean"  )
      if (is.null(H)) next ()
      for (w in 1:p$nw ) {
        wchar = paste( "0", w, sep="" )
        wchar = substr( wchar, nchar(wchar)-1, nchar(wchar) )
        outfn = paste( "temperatures.bottom", y, wchar, sep=".")
        annot = paste("Temperature\n", y, " - ", w, sep="")
        bio.spacetime::map( xyz=cbind(loc, H[,w]), cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
          fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
          corners=p$corners, spatial.domain=p$spatial.domain )
      }
    }
    return( "Completed maps")
  }

  # ---------------------------

  if ( type %in% c("annual" ) ) {
    bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain , "bottom.predictions", "annual" )
    dir.create( bottomdir.maps, recursive=T, showWarnings=F )
    for (iy in ip ) {
      y = p$runs[iy, "yrs"]
      print(y)
      H = temperature.db( p=p, DS="bottom.statistics.annual", yr=y, ret=vname )

      if (is.null(H)) next ()

      if (vname %in% c("tmean") ) {
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom", y, sep=".")
        annot = y
      } 

      if (vname %in% c("tsd") ) {
        datarange = seq(0.001, 6, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.sd", y, sep=".")
        annot = y
      } 
   
      if (vname %in% c("tmin") ) {
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.min", y, sep=".")
        annot = y
      }
   
      if (vname %in% c("tmax") ) {
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.max", y, sep=".")
        annot = y
      }

      if (vname %in% c("amplitude") ) {
        datarange = seq(0,5, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.amplitude", y, sep=".")
        annot = y
      }

      bio.spacetime::map( xyz=cbind(loc, H[,iy]), cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 

    } 
  }

  # ------------------------------

  if ( type %in% c("climatology" ) ) {

    bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain, "bottom.predictions", "climatology" )
    dir.create( bottomdir.maps, recursive=T, showWarnings=F )

    H = temperature.db( p=p, DS="bottom.statistics.climatology" )

    for (vname in p$bstats) {
      if (vname %in% c("tmean") ) {
        # datacols = c("plon", "plat", "tmean")
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom", sep=".")
        annot = paste("Temperature bottom climatology\n", sep="")
      } 

      if (vname %in% c("tsd") ) {
        # datacols = c("plon", "plat", "tmean")
        datarange = seq(0.001, 6, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.sd", sep=".")
        annot = paste("Temperature bottom SD climatology\n", sep="")
      } 
   
      if (vname %in% c("tmin") ) {
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.min", sep=".")
        annot = paste("Temperature bottom minimum climatology\n",  sep="")
      }
   
      if (vname %in% c("tmax") ) {
        datarange = seq(-0.5, 10, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.max", sep=".")
        annot = paste("Temperature bottom maximum climatology\n",  sep="")
     }

      if (vname %in% c("amplitude") ) {
        datarange = seq(0,5, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "temperatures.bottom.amplitude", sep=".")
        annot = paste("Temperature bottom amplitude climatology\n",  sep="")
      }

      bio.spacetime::map( xyz=cbind(loc, H[,which( p$bstats==vname)]), cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 
      print( file.path( bottomdir.maps,outfn))
    }

  }


  # ------------------------------


  if ( type %in% c("lbm.stats" ) ) {

    bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain, "bottom.predictions", "lbm.stats" )
    dir.create( bottomdir.maps, recursive=T, showWarnings=F )

    H = temperature.db( p=p, DS="lbm.stats" )

    for (vname in colnames(H)) {
      
      if (vname %in% c("sdTotal") ) {
        datarange = seq(0.1, 6, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "sdTotal", sep=".")
        annot = paste("Temperature bottom sdTotal\n", sep="")
        xyz=cbind(loc, H[,vname])
      } 

      if (vname %in% c("rsquared") ) {
        datarange = seq(0.1, 0.99, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "rsquared", sep=".")
        annot = paste("Temperature bottom rsquared\n", sep="")
        xyz=cbind(loc, H[,vname])
      } 
   
      if (vname %in% c("ndata") ) {
        datarange = seq(1, 8000, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "ndata", sep=".")
        annot = paste("Temperature bottom ndata\n",  sep="")
        xyz=cbind(loc, H[,vname])
      }
   
      if (vname %in% c("sdSpatial") ) {
        datarange = seq(0.1, 7, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "sdSpatial", sep=".")
        annot = paste("Temperature bottom sdSpatial\n",  sep="")
        xyz=cbind(loc, H[,vname])
     }

      if (vname %in% c("sdObs") ) {
        datarange = seq(0, 4, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "sdObs", sep=".")
        annot = paste("Temperature bottom sdObs\n",  sep="")
        xyz=cbind(loc, (H[,vname]))
      }

      if (vname %in% c("range") ) {
        datarange = seq( log(0.5), log(500), length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "range", sep=".")
        annot = paste("Temperature bottom range\n",  sep="")
        xyz=cbind(loc, log(H[,vname]))
        xyz = xyz[which( is.finite(rowSums(xyz))),]
      }

      if (vname %in% c("phi") ) {
        datarange = seq(log(0.15), log(100), length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "phi", sep=".")
        annot = paste("Temperature bottom phi\n",  sep="")
        xyz=cbind(loc, log(H[,vname]))
        xyz = xyz[which( is.finite(rowSums(xyz))),]
      }

      if (vname %in% c("nu") ) {
        datarange = seq(0.1, 3.5, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "nu", sep=".")
        annot = paste("Temperature bottom nu\n",  sep="")
        xyz=cbind(loc, (H[,vname]))
      }

      if (vname %in% c("ar_timerange") ) {
        datarange = seq(0, 6, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "ar_timerange", sep=".")
        annot = paste("Temperature bottom ar_timerange\n",  sep="")
        xyz=cbind(loc, (H[,vname]))
      }

      if (vname %in% c("ar_1") ) {
        datarange = seq(0, 1, length.out=100)
        cols = color.code( "blue.black", datarange )
        outfn = paste( "ar_1", sep=".")
        annot = paste("Temperature bottom ar_1\n",  sep="")
        xyz=cbind(loc, H[,vname])
      }

      bio.spacetime::map( xyz=xyz, cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
        corners=p$corners, spatial.domain=p$spatial.domain ) 
      print( file.path( bottomdir.maps,outfn))

    }

  }



  return (NULL)
}



