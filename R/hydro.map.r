
  hydro.map = function( ip=NULL, p=NULL, type="annual", ... ) {

    # ip is the first parameter passed in the parallel mode
    if (exists( "libs", p)) RLibrary( p$libs )
    if (is.null(ip)) ip = 1:p$nruns

		require( lattice )
		datarange = seq(-0.5, p$nw, length.out=50)
    cols = color.code( "blue.black", datarange )

    if ( type=="dyear" ) {
      warning( "This only works for the default spatial domain .. ")
      bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain, "bottom.predictions", "dyear" )
      dir.create( bottomdir.maps, recursive=T, showWarnings=F )
      datarange = seq(-0.5, p$nw, length.out=50)
      cols = color.code( "blue.black", datarange )
      xyz = bathymetry.db( p=p, DS="baseline" )
      xyz = xyz[, c("plon", "plat")]
      for (iy in ip ) {
        y = p$runs[iy, "yrs"]
        H = temperature.db( p=p, DS="temporal.interpolation", yr=y  )
        if (is.null(H)) next ()
        for (w in 1:p$nw ) {
          wchar = paste( "0", w, sep="" )
          wchar = substr( wchar, nchar(wchar)-1, nchar(wchar) )
          outfn = paste( "temperatures.bottom", y, wchar, sep=".")
          annot = paste("Temperature\n", y, " - ", w, sep="")
          map( xyz=cbind(xyz, H[,w]), cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners , spatial.domain=p$spatial.domain)
        }
      }
      return( "Completed maps")
    }

    if ( type %in% c("annual", "amplitudes", "temperatures", "tsd" ) ) {

      bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain , "bottom.predictions", "annual" )
      dir.create( bottomdir.maps, recursive=T, showWarnings=F )

      for (iy in ip ) {
        y = p$runs[iy, "yrs"]
        print(y)

        H = temperature.db( p=p, DS="complete", year=y )

        if ( p$spatial.domain=="snowcrab" ) {
          i = which( H$plon< 990 &  H$plon > 220  &   ## these are in planar coords  ..should fix this hack one day
                     H$plat< 5270 &  H$plat > 4675
          )
          H = H[ i, ]
        }
        if (is.null(H)) next ()

        if (type %in% c("temperatures", "annual") ) {
          datacols = c("plon", "plat", "tmean")
          datarange = seq(-1,11, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom", y, sep=".")
          annot = y
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners, spatial.domain=p$spatial.domain )
        }

        if (type %in% c("amplitudes", "annual") ) {
          datacols = c("plon", "plat", "tamp")
          datarange = seq(0,10, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom.amplitude", y, sep=".")
          annot = y
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners, spatial.domain=p$spatial.domain   )
        }


     
        if (type %in% c("tsd", "annual") ) {
          datacols = c("plon", "plat", "tsd")
          datarange = seq(0, 4, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom.sd", y, sep=".")
          annot = y
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners , spatial.domain=p$spatial.domain )
        }

      } }

      if ( type %in% c("global", "amplitudes", "temperatures", "tsd" ) ) {

      bottomdir.maps = file.path( project.datadirectory("bio.temperature"), "maps", p$spatial.domain, "bottom.predictions", "global" )
      dir.create( bottomdir.maps, recursive=T, showWarnings=F )

      H = temperature.db( p=p, DS="complete", year=p$tyears[1] )

        if (type %in% c("temperatures", "global") ) {
          datacols = c("plon", "plat", "tmean.climatology")
          datarange = seq(-1,11, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom", sep=".")
          annot = paste("Temperature\n", sep="")
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners, spatial.domain=p$spatial.domain  )
        }

        if (type %in% c("amplitudes", "global") ) {
          datacols = c("plon", "plat", "tamplitude.climatology")
          datarange = seq(0,10, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom.amplitude", sep=".")
          annot = paste("Temperature amplitude\n",  sep="")
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners , spatial.domain=p$spatial.domain )
        }

  
       if (type %in% c("tsd", "global") ) {
          datacols = c("plon", "plat", "tsd.climatology")
          datarange = seq(0, 5, length.out=50)
          cols = color.code( "blue.black", datarange )
          outfn = paste( "temperatures.bottom.sd",  sep=".")
          annot = paste("Temperature SD\n", sep="")
          map( xyz=H[,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
            fn=outfn, loc=bottomdir.maps, at=datarange , col.regions=cols,
            corners=p$corners  , spatial.domain=p$spatial.domain)
        }

    }
    return (NULL)
  }



