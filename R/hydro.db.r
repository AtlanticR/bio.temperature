
hydro.db = function( ip=NULL, p=NULL, DS=NULL, yr=NULL, additional.data=c("groundfish", "snowcrab", "USSurvey_NEFSC"), ...) {

  # manipulate temperature databases from osd, groundfish and snow crab and grid them
  # OSD data source is
  # http://www.meds-sdmm.dfo-mpo.gc.ca/zmp/climate/climate_e.htm
  # http://www.mar.dfo-mpo.gc.ca/science/ocean/database/data_query.html
  ## must download manually to this directory and run gzip
  ## use choijae/Jc#00390
  ## depths: 500,500, "complete profile"   .. raw data  for the SS
  # (USER Defined -- region: jc.ss")

  # no time records, just day/mon/year .. assume utc

  basedir = file.path( project.datadirectory("bio.temperature"), "data" )
  loc.archive = file.path( basedir, "archive", "profiles")
  loc.basedata = file.path( basedir, "basedata", "rawdata" )
  dir.create( loc.basedata, recursive=T, showWarnings=F )

  # OSD data series variables of interest


  if ( DS == "osd.rawdata" ) {
    # simple loading of annual data files
    out = NULL
    for ( y in yr ) {
      print (y)
      fn = file.path( loc.basedata, paste( "osd.rawdata", y, "rdata", sep=".") )
      if (file.exists ( fn ) ) {
        load(fn)
        out = rbind( out, X )
      }
    }
    return ( out )
  }


  if ( DS=="osd.rawdata.allfiles.redo" ) {
    fn.all = list.files( path=loc.archive, pattern="osd.clim.*.gz", full.names=T)
    X = NULL
    varlist = c("DEPTH","PRESSURE","CRUISE_DATE","LATITUDE" ,"LONGITUDE" ,"TEMPERATURE","SALINITY" ,"SIGMAT" )
    varstosave = c( "depth", "pressure", "latitude" ,"longitude" ,"temperature" ,"salinity" ,"sigmat", "date" )
    for (fn in fn.all) {
      f = read.csv( gzfile(fn), header=T, as.is=T, sep=",", na.strings="9999")
      f = f[,varlist]
      fyears = as.numeric( matrix( unlist( strsplit( f$CRUISE_DATE, "/" ) ), ncol=3, byrow=T) [,3] )
      years = sort( unique( fyears ))
      for (yrs in years) {
        fn.out = file.path( loc.basedata,  paste( "osd.rawdata", yrs, "rdata", sep=".") )
        print( paste(yrs, ":", fn.out) )
        X = f[ which( fyears == yrs) ,]
        names(X) = tolower( names(X) )
        X$date = lubridate::dmy( X$cruise_date )   # default is UTC ... need to confirm-- no time records .. assume utc

        X = X[ , varstosave ]
        save( X, file=fn.out, compress=T)
      }
    }

  }

  if (DS=="osd.rawdata.singleyear.redo" ) {
    varlist = c("DEPTH","PRESSURE","CRUISE_DATE","LATITUDE" ,"LONGITUDE" ,"TEMPERATURE","SALINITY" ,"SIGMAT" )
    varstosave = c( "depth", "pressure", "latitude" ,"longitude" ,"temperature" ,"salinity" ,"sigmat", "date" )
    for ( y in yr) {
      X = NULL
      fn.all = list.files( path=loc.archive, pattern="osd.clim.*.gz", full.names=T)
      fn = fn.all[ grep (as.character(y), fn.all) ]
      f = read.csv( gzfile(fn), header=T, as.is=T, sep=",", na.strings="9999")
      X = f[,varlist]
      fn.out = file.path( loc.basedata, paste( "osd.rawdata", y, "rdata", sep=".") )
      names(X) = tolower( names(X) )
      X$date = lubridate::ymd( X$cruise_date, tz= )  # default is UTC ... need to confirm
      X= X[, varstosave ]
      save( X, file=fn.out, compress=T)
    }
  }

  # ----------------------

  if (DS=="osd.initial" ) {
    ## this is a data dump directly from Roger Pettipas for 2008 to 2015
    varstosave = c( "depth", "pressure", "latitude" ,"longitude" ,"temperature" ,"salinity" ,"sigmat", "date" )
    fndata = file.path( loc.archive, "Data_2008-2015.csv.xz" )
    XX = read.csv( file=xzfile(fndata), header=FALSE, skip=2 , stringsAsFactors=FALSE, na.strings="9999" )
    header = c("MissionID", "Latitude", "Longitude", "Year", "Month", "Day", "Hour", "Minute", "Pressure", "Temperature", "Salinity", "SigmaT" ,"StationID" )
    names(XX) = tolower( header )
    XX$depth = decibar2depth ( P=XX$pressure, lat=XX$latitude )
    if (!exists( "sigmat", XX))  XX$sigmat = XX$sigma.t  # naming is variable
    XX$date_string = paste( XX$year, XX$month, XX$day, sep="-" )
    XX$date = lubridate::ymd( XX$date_string )   # default is UTC ... need to confirm

    yrs = sort( unique( XX$year) )
    for ( y in yrs ) {
      print (y)
      fn.out = file.path( loc.basedata, paste( "osd.rawdata", y, "rdata", sep=".") )
      ii = which ( XX$year == y )
      if (length(ii) > 1) {
        X= XX[ ii, varstosave ]
        save( X, file=fn.out, compress=T)
      }
    }
  }


  # ----------------------


  if (DS=="osd.current" ) {
    ## this is a data dump directly from Roger Pettipas for 2015 and on
    varstosave = c( "depth", "pressure", "latitude" ,"longitude" ,"temperature" ,"salinity" ,"sigmat", "date" )
    for ( y in yr ) {
      print (y)
      fndata = file.path( loc.archive, paste( "Data_", y, ".csv.xz", sep="" ) )
      fn.out = file.path( loc.basedata, paste( "osd.rawdata", y, "rdata", sep=".") )
      X = read.csv( file=xzfile(fndata), skip=2, stringsAsFactors=FALSE, na.strings="9999" )
      # insert Header :
      header = c("MissionID", "Latitude", "Longitude", "Year", "Month", "Day", "Hour", "Minute", "Pressure", "Temperature", "Salinity", "SigmaT" ,"StationID" )
      names(X) = tolower( header )
      X$depth = decibar2depth ( P=X$pressure, lat=X$latitude )
      if (!exists( "sigmat", X))  X$sigmat = X$sigma.t  # naming is variable
      X$date_string = paste( X$year, X$month, X$day, sep="-" )
      X$date = lubridate::ymd( X$date_string )   # default is UTC ... need to confirm

      X= X[, varstosave ]
      save( X, file=fn.out, compress=T)
    }
  }

  ##

  if (DS=="USSurvey_NEFSC") {
    # data dump supplied by Adam Cook .. assumed to tbe bottom temperatures from their surveys in Gulf of Maine area?
    ne = NULL
    fn = file.path( project.datadirectory("bio.temperature"), "archive", "NEFSCTemps.rdata" )
    if (file.exists(fn)) load(fn)
    ne$id = paste(ne$plon, ne$plat, lubridate::date( ne$timestamp), sep="~" )
    ne$salinity = NA
    ne$oxyml = NA
    ne$sigmat = NA
    ne$date = ne$timestamp
    ne$yr = lubridate::year( ne$timestamp )
    ne$dyear = lubridate::decimal_date( ne$timestamp ) - ne$yr
    ne = planar2lonlat( ne, proj.type=p$internal.projection )  # convert lon lat to coord system of p0
    if (is.null(yr)) return(ne) # everything
    i = which( lubridate::year( ne$timestamp) %in% yr )
    if (length(i) > 0) ne = ne[i,]
    return (ne)
  }


  # ----------------

  if ( DS %in% c("ODF_ARCHIVE", "ODF_ARCHIVE.redo") ) {
    # PTRAN/CHOIJ
    loc = file.path( project.datadirectory("bio.temperature"), "data" )

    DataDumpFromWindows = F
    if ( DataDumpFromWindows ) {
      loc = file.path("C:", "datadump")
    }
    dir.create( path=loc, recursive=T, showWarnings=F )

    fn.root =  file.path( loc, "ODF_ARCHIVE" )
    dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

    out = NULL
    if ( is.null(DS) | DS=="ODF_ARCHIVE" ) {
      fl = list.files( path=fn.root, pattern="*.rdata", full.names=T )
        for ( fny in fl ) {
        load (fny)
        out = rbind( out, odfdat )
      }
      return (out)
    }

    con = ROracle::dbConnect( DBI::dbDriver("Oracle"), username=oracle.personal.user, password=oracle.personal.password, dbname="PTRAN" )
    cruises   <- ROracle::dbGetQuery(con, "select * from ODF_ARCHIVE.ODF_CRUISE_EVENT" )

    for ( y in yr ) {
      fny = file.path( fn.root, paste( y, "rdata", sep="."))
      odfdat = ROracle::dbGetQuery( con,  paste(
      " select * " ,
      " from ODF_ARCHIVE.ODF_CRUISE_EVENT i, ODF_ARCHIVE.ODF_DATA j " ,
      " where i.CRUISE_EVENT_ID(+)=j.DATA_VAL_ID ",
      " and EXTRACT(YEAR from start_date_time) =", y, ";"
      ) )

      names(odfdat) =  tolower( names(odfdat) )
      print(fny)
      save(odfdat, file=fny, compress=T)
      gc()  # garbage collection
      print(y)
    }

    ROracle::dbDisconnect(connect)
    
    return (fn.root)

  }

  # ----------------

  if (DS %in% c( "profiles.annual.redo", "profiles.annual" ) ) {
    # read in annual depth profiles then extract bottom temperatures

    basedir = project.datadirectory("bio.temperature", "data" )
    loc.profile = file.path( basedir, "basedata", "profiles" )
    dir.create( loc.profile, recursive=T, showWarnings=F )

    if (DS=="profiles.annual") {
      fn = file.path(  loc.profile, paste("depthprofiles", yr, "rdata", sep="."))
      Y = NULL
			if (file.exists( fn) ) load (fn )
      return(Y)
    }

    ####### "ip" is the first parameter expected when run in parallel mode .. do not move this one

      if ( is.null(ip)) {
        if( exists( "nruns", p ) ) {
          ip = 1:p$nruns
        } else {
          if ( !is.null(yr)) {
            # if only selected years being re-run
            ip = 1:length(yr)
            p$runs = data.frame(yrs = yr)
          } else {
            ip = 1:length(p$tyears)
            p$runs = data.frame(yrs = p$tyears)
          }
        }
      }

    if (exists( "libs", p)) RLibrary( p$libs )

    # bring in snow crab, groundfish and OSD data ...
    set = bio.snowcrab::snowcrab.db( DS="setInitial" )
    mlu = bio.snowcrab::minilog.db( DS="set.minilog.lookuptable" )
    slu = bio.snowcrab::seabird.db( DS="set.seabird.lookuptable" )
    set = merge( set, mlu, by= c("trip", "set"), all.x=TRUE, all.y=FALSE )
    set = merge( set, slu, by= c("trip", "set"), all.x=TRUE, all.y=FALSE )
    set$longitude =set$lon
    set$latitude = set$lat
    set$oxyml = NA
    set$salinity = NA
    set$sigmat = NA

    set = set[ ,c("minilog_uid", "seabird_uid", "longitude", "latitude", "oxyml", "salinity", "sigmat" ) ]

    grdfish = bio.groundfish::groundfish.db( "gshyd.georef" )

    Ydummy = bio.temperature::hydro.db( DS="osd.rawdata", yr=2000, p=p ) [1,]  # dummy entry using year=2000
    Ydummy$yr = NA
    Ydummy$dyear = 0.5
    Ydummy$id =  "dummy"
    Ydummy$depth = -1
    Ydummy$oxyml = NA

    dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks

    for (iy in ip) {
      yt = p$runs[iy, "yrs"]

      Y =  bio.temperature::hydro.db( DS="osd.rawdata", yr=yt, p=p )
        if ( is.null(Y) ) {
          Y = Ydummy
          Y$yr = yt
        } else {
          Y$yr = yt
          Y$dyear = lubridate::decimal_date( Y$date ) - Y$yr

          Yid = cut( Y$dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE )
          Y$id =  paste( round(Y$longitude,2), round(Y$latitude,2), Yid , sep="~" )
          Y$depth = bio.utilities::decibar2depth ( P=Y$pressure, lat=Y$latitude )
          Y$oxyml = NA
          # next should not be necessary .. but just in case the osd data types get altered
          Y$temperature = as.numeric(Y$temperature )
          Y$salinity= as.numeric(Y$salinity)
          Y$sigmat = as.numeric(Y$sigmat)
        }

      Y$pressure = NULL

      if ("groundfish" %in% additional.data ) {
        gfkeep = c( "id", "sdepth", "temp", "sal", "oxyml", "lon", "lat", "yr", "timestamp")
        gf = grdfish[ which( grdfish$yr == yt ) , gfkeep ]
				if (nrow(gf) > 0) {
					gf$sigmat = NA
					gf$date = gf$timestamp
          # gf$date = as.POSIXct(gf$date, origin=lubridate::origin)
          gf$dyear = lubridate::decimal_date( gf$date ) - gf$yr
          names(gf) = c( "id", "depth", "temperature", "salinity", "oxyml", "longitude", "latitude", "yr", "date", "dyear", "sigmat"  )
          Y = rbind( Y, gf[, names(Y)] )
				}
			}

      if ("snowcrab" %in% additional.data ) {

        minilog = bio.snowcrab::minilog.db( DS="basedata", Y=yt )

        if (! is.null( nrow( minilog ) ) ) {
          minilog = merge( minilog, set, by="minilog_uid", all.x=TRUE, all.y=FALSE )
          minilog$id = minilog$minilog_uid
          minilog$date = minilog$timestamp
          # minilog$date = as.POSIXct(minilog$chron, origin=lubridate::origin)
          minilog$yr = yt
          minilog$dyear = lubridate::decimal_date( minilog$date ) - minilog$yr
          Y = rbind( Y, minilog[, names(Y) ] )
        }

        seabird = bio.snowcrab::seabird.db( DS="basedata", Y=yt )
				if ( !is.null( nrow( seabird ) ) ) {
          seabird = merge( seabird, set, by="seabird_uid", all.x=TRUE, all.y=FALSE )
          seabird$id = seabird$seabird_uid
          seabird$yr = yt
          seabird$date = seabird$timestamp
          # seabird$date = as.POSIXct(seabird$chron, origin=lubridate::origin)
          seabird$dyear = lubridate::decimal_date( seabird$date ) - seabird$yr
          seabird$oxyml = NA
          Y = rbind( Y, seabird[, names(Y) ] )
        }

      }

      oo = which( Y$id == "dummy" )
      if (length(oo) > 0 ) Y = Y[ -oo, ]

      if ( is.null( nrow(Y) ) ) next()
      if ( nrow(Y) < 5 ) next()

      if ( is.null(Y) ) next()

      iiY = which(duplicated(Y))
      if (length(iiY)>0) Y = Y [ -iiY, ]

      bad = which( Y$temperature < -5 | Y$temperature > 30 )
      if (length(bad)>0) Y=Y[-bad,]

      fn = file.path( loc.profile, paste("depthprofiles", yt, "rdata", sep="."))
      print( fn )
      save( Y, file=fn, compress=T )
    }

    return ("Completed")
  }


  # ----------------

  if (DS %in% c( "bottom.annual", "bottom.annual.redo", "bottom.all" ) ) {
    # extract bottom temperatures

    basedir = project.datadirectory("bio.temperature", "data" )
    loc.bottom = file.path( basedir, "basedata", "bottom"  )
    dir.create( loc.bottom, recursive=T, showWarnings=F )
    
    fbAll = file.path( loc.bottom, "bottom.all.rdata" ) 
    if (DS=="bottom.all") {
      O = NULL
      if (file.exists(fbAll) ) load (fbAll)
      return(O)
    }

    if (DS=="bottom.annual") {
      fn = file.path( loc.bottom, paste("bottom", yr, "rdata", sep="."))
      Z = NULL
			if (file.exists(fn) ) load (fn )
      return(Z)
    }

      if ( is.null(ip)) {
        if( exists( "nruns", p ) ) {
          ip = 1:p$nruns
        } else {
          if ( !is.null(yr)) {
            # if only selected years being re-run
            ip = 1:length(yr)
            p$runs = data.frame(yrs = yr)
          } else {
            ip = 1:length(p$tyears)
            p$runs = data.frame(yrs = p$tyears)
          }
        }
      }

    if (exists( "libs", p)) RLibrary( p$libs )

    tne = hydro.db( p=p, DS="USSurvey_NEFSC" )

    dyears = (c(1:(p$nw+1))-1)  / p$nw # intervals of decimal years... fractional year breaks

    for (iy in ip) {

      yt = p$runs[iy, "yrs"]
      Y = bio.temperature::hydro.db( DS="profiles.annual", yr=yt, p=p )
      if (is.null(Y)) next()
      igood = which( Y$temperature >= -3 & Y$temperature <= 25 )  ## 25 is a bit high but in case some shallow data
      Y = Y[igood, ]

      # Bottom temps
      Yid = cut( Y$dyear, breaks=dyears, include.lowest=T, ordered_result=TRUE )
      Y$id =  paste( round(Y$longitude,2), round(Y$latitude,2), Yid, sep="~" )
      ids =  sort( unique( Y$id ) )
      Z = copy.data.structure( Y)

      for (i in ids ) {
        W = Y[ which( Y$id == i ), ]
        jj = which( is.finite( W$depth ) )
        if ( length(jj) < 3 ) next()
        Wmax = max( W$depth, na.rm=T ) - 10  # accept any depth within 10 m of the maximum depth
        kk =  which( W$depth >= Wmax )
        R = W[ which.max( W$depth ) , ]
        R$temperature = median( W$temperature[kk] , na.rm=T )
        R$salinity = median( W$salinity[kk] , na.rm=T )
        R$sigmat = median( W$sigmat[kk] , na.rm=T )
        R$oxyml = median( W$oxyml[kk] , na.rm=T )
        Z = rbind( Z, R )
      }

      Z = rename.df( Z, "longitude", "lon")
      Z = rename.df( Z, "latitude", "lat")
      Z = rename.df( Z, "temperature", "t")
      Z = rename.df( Z, "depth", "z")

      utne = which( tne$yr== yt)
      if ( length(utne) > 0) Z = rbind( Z, tne[,names(Z)] )

      Z$date = as.Date( Z$date ) # strip out time of day information
      Z$ddate = lubridate::decimal_date( Z$date )
      Z$dyear = Z$ddate - Z$yr
      Z = lonlat2planar( Z, proj.type=p$internal.projection )

      igood = which( Z$t >= -3 & Z$t <= 25 )  ## 25 is a bit high but in case some shallow data
      Z = Z[igood, ]

      igood = which( Z$lon >= p$corners$lon[1] & Z$lon <= p$corners$lon[2]
          &  Z$lat >= p$corners$lat[1] & Z$lat <= p$corners$lat[2] )
      Z = Z[igood, ]


      Z = Z[ which( is.finite( Z$lon + Z$lat + Z$plon + Z$plat ) ) , ]
      ## ensure that inside each grid/time point
      ## that there is only one point estimate .. taking medians
      vars = c("z", "t", "salinity", "sigmat", "oxyml")
      Z$st = paste( Z$ddate, Z$plon, Z$plat )

      o = which( ( duplicated( Z$st )) )
      if (length(o)>0) {
        dupids = unique( Z$st[o] )
        for ( dd in dupids ) {
          e = which( Z$st == dd )
          keep = e[1]
          drop = e[-1]
          for (v in vars) Z[keep, v] = median( Z[e,v], na.rm=TRUE )
          Z$st[drop] = NA  # flag for deletion
        }
        Z = Z[ -which( is.na( Z$st)) ,]
      }
      Z$st = NULL
      fn = file.path( loc.bottom, paste("bottom", yt, "rdata", sep="."))
			print (fn)
      save( Z, file=fn, compress=T)
    }

    O = NULL
    for ( yr in p$tyears ) {
      o = hydro.db( p=p, DS="bottom.annual", yr=yr )
      if (!is.null(o)) O = rbind(O, o)
    }
    save(O, file=fbAll, compress=TRUE)

    return ("Completed")

  }


}


