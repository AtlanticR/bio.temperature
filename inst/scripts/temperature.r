
  # ----------------
  # Prep OSD, snow crab and groundfish temperature profiles
  # this one has to be done manually .. no longer mainted by anyone ..

  p = list( project.name = "bio.temperature" )
  p$project.root = project.datadirectory( p$project.name )

  p$libs = RLibrary( c( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "fields" ) )
  p$init.files = bioLibrary(  "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" )

  # p$tyears = c(1910:2013)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
  p$tyears = c(1950:2015)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
  p$ny = length(p$tyears)
  p$nw = 10 # number of intervals in time within a year
  p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
  p$gam.optimizer = "nlm" ## other optimizers:: "bam" (slow), "perf"(ok), "nlm" (good), "bfgs" (ok), "newton" (default)
  p$nMin.tbot = p$ny*3 # min number of data points req before attempting to model timeseries in a localized space
  p$dist.km = c( 2.5, 5, 7.5, 10, 12.5, 15 ) # "manhattan" distances to extend search for data
  p$maxdist = 20 # if using gstat  max dist to interpolate in space
  # choose: temporal interpolation method ... harmonic analysis seems most reasonable
  # .. do not use more than 2 as it chases noise too much .. 1 harmonic seems the best in terms of not chasing after noise
  # possible methods: "annual", "seasonal.basic", "seasonal.smoothed", "harmonics.1", "harmonics.2", "harmonics.3", "inla.ts.simple"
  p$tsmethod = "harmonics.1"
  # p$spmethod = "inverse.distance"  ## too slow
  # p$spmethod = "gam" ## too smooth
  p$spmethod = "kernel.density" ## best
  p$theta = 5 # FFT kernel bandwidth (SD of kernel) for method p$spmethod = "kernel.density"
  p$nsd = 6 # number of SD distances to pad boundaries with 0 for FFT  in method  p$spmethod = "kernel.density

  p$newyear = 2015

  p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution
  p$spatial.domain.default = "canada.east"
  p = spatial.parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution

  # ------------------------------

  if ( create.baseline.database ) {

    if (historical.data.redo) {
      hydro.db( DS="osd.rawdata.allfiles.redo", p=p )   # redo whole data set (historical) from 1910 to 2010
      hydro.db( DS="osd.initial", p=p ) # 2008:2014
    }
    # Roger Petipas has been maintaining a database, the following loads this data
    hydro.db( DS="osd.current", p=p, yr=2014:p$newyear ) # specify range or specific year

    # Merge depth profiles from all data streams: OSD, groundfish, snowcrab
    p = make.list( list( yrs=c(2008:p$newyear)), Y=p )   # specify range or specific year
    p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce ep = make.list( list( yrs=c(2008:p$newyear), Y=p ))   # specify range or specific year
    hydro.db( DS="profiles.annual.redo", p=p  )  # specify range or specific year
    # parallel.run( hydro.db, p=p, yr=p$tyears, DS="profiles.annual.redo", init.files=p$init.files )

    # Extract bottom data from each profile
    p = make.list( list( yrs=2008:p$newyear), Y=p )  # specify range or specific year
    hydro.db( DS="bottom.annual.redo", yr=2008:p$newyear, p=p ) # yr argument overrides p$tyears .. e.g. for a new year of data
    # hydro.db( DS="bottom.annual.redo", p=p )
    # parallel.run( hydro.db, p=p, yr=p$tyears, DS="bottom.annual.redo", init.files=p$init.files )
  }

  # ------------------------------
  # Basic data uptake now complete  .. move to interpolations
  # ------------------------------

  if (create.interpolated.results ) {

    # 1. grid bottom data to internal spatial resolution ; <1 min
    p = make.list( list( yrs=p$tyears), Y=p )
    # parallel.run( hydro.db, p=p, DS="bottom.gridded.redo" )
    hydro.db( p=p, DS="bottom.gridded.redo" )  # all p$tyears, for a single year use with yr argument: yr=p$newyear
    hydro.db( p=p, DS="bottom.gridded.all.redo" )  # all p$tyears, for a single year use with yr argument: yr=p$newyear


    # 2. temporal interpolations assuming some seasonal pattern
    # 1950-2013, SSE took ~ 35 hrs on laptop (shared RAM, 24 CPU; 1950-2013 run April 2014 ) ... 17 GB req of shared memory
    # 1950-2015, SSE 22 hrs, 42 GB RAM, 8 CPU on hyperion (10 Jan 2015), using NLM .. not much longer for "canada.east"
    # define output mattrix
    # predictions are made upon the locations defined by bathymetry "baseline"
    p$nP = nrow( bathymetry.db( p=p, DS="baseline" ) )
    p = temperature.db(p=p, DS="bigmemory.initiate" ) # return pointers to bigmemory objects

      # test to make sure attach is possible ..
      # sometimes requires restart of R or even computer .. bigmemory memory leak or file caching
      tbot <- bigmemory::attach.big.matrix( p$descriptorfile.tbot  )
      tbot.se <- bigmemory::attach.big.matrix( p$descriptorfile.tbotse  )

    p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce efficiency using the beowulf network
    # p = make.list( list( loc=sample.int( p$nP) ), Y=p ) # random order helps use all cpus
    p = make.list( list( loc=sample.int( p$nP ) ), Y=p ) # random order helps use all cpus
    parallel.run( temperature.timeseries.interpolate, p=p)
    # temperature.timeseries.interpolate ( p=p )
    temperature.db( p=p, DS="temporal.interpolation.redo" ) # save interpolation as time slices to disk


    # 3. simple spatial interpolation
    # NOTE this one does a regridding as selected by p$subregions ...
    # ... it is required for the habitat lookup .. no way around it
    # (complex/kriging takes too much time/cpu) ==> 3-4 hr/run
    # using localhost in 2014 6+ hr for each run but with multiple cycles ~ 10 hr total
    # use all clusters if available
    p$clusters = c( rep("kaos",23), rep("nyx",24), rep("tartarus",24) )
    p$clusters = rep("localhost", detectCores() )
    p = make.list( list( yrs=p$tyears), Y=p )
    parallel.run( temperature.db, p=p, DS="spatial.interpolation.redo" )
    #  temperature.db( p=p, DS="spatial.interpolation.redo" ) # 2hr in serial mode


    # 4. extract relevant statistics:: only for default grid . TODO might as well do for each subregion/subgrid
    # temperature.db(  p=p, DS="bottom.statistics.annual.redo" )
    # or parallel runs: ~ 1 to 2 GB / process
    # 4 cpu's ~ 10 min
    p$clusters = c( rep("kaos",23), rep("nyx",24), rep("tartarus",24) )
    p = make.list( list( yrs=p$tyears), Y=p )
    parallel.run( temperature.db, p=p, DS="bottom.statistics.annual.redo" )
    #  temperature.db( p=p, DS="bottom.statistics.annual.redo" )



    # 5. climatology database ... ~ 2 min :: only for  default grid . TODO might as well do for each subregion/subgrid
    p$bstats = c("tmean", "tamplitude", "wmin", "thalfperiod", "tsd" )
    p$tyears.climatology = p$tyears  # or redefine it with : p$tyears.climatology = 1950:2015
    temperature.db ( p=p, DS="climatology.redo")

  }


  ### to this point everything is run on p$spatial.domain.default domain (except spatial.interpolation), now take subsets:

    # 7. downscale to appropriate domain: simple interpolations and maps
    # p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce efficiency using th
    # p$clusters = c( rep("kaos",23), rep("nyx",24), rep("tartarus",24) )
    p = make.list( list( yrs=p$tyears), Y=p )

    for ( gr in p$subregions ) {
      print (gr)
      p = spatial.parameters(  p=p, type= gr )
      if ( length(p$clusters) > 1 ) {
        parallel.run( temperature.db, p=p, DS="complete.redo")
        parallel.run( hydro.map, p=p, type="annual"  )
        parallel.run( hydro.map, p=p, type="global")
      } else {
        temperature.db( p=p, DS="complete.redo")
        hydro.map( p=p, type="annual" )
        hydro.map( p=p, type="global" )
      }
    }



  # finished interpolations


  #-----------------------------------------
  # the following is just for testing / notes

  if (0) {
    # to access data:
    p = spatial.parameters( p=p, type="SSE" )
    tp = hydro.db(  p=p, DS="profiles.annual", yr=2007 )
    tp = hydro.db(  p=p, DS="bottom.annual", yr=2007 )


    #  test analysis spatial variations in temperature
    require( gstat )
    p = spatial.parameters( p=p, type="SSE" )
    testyear = 2006
    O = hydro.db( p=p, DS="bottom.gridded", yr=testyear )
    O = O[, c("plon", "plat", "t", "yr", "dyr")]
    O = O[ which( is.finite(O$t)) ,]
    ee.g = gstat( id = "t", formula=t~plon+plat, locations=~plon+plat, data=O )
    vg <- variogram( ee.g, boundaries=c(2, 4, 8, 16, 20, 32, 40, 64, 80, 128, 150 ) )
    vm0 = vgm( psill=max(vg$gamma) * 0.75, model="Mat", range=max(vg$dist)*0.5, nugget=max(vg$gamma) * 0.25, kappa = 2 ) #initial settings
    vg.fit = fit.variogram( vg, vm0 )
    plot( vg, model=vg.fit, plot.numbers=T )
    # running various models suggest on dates with good data :
    #   range: about 40 km
    #   nugget: 1.5
    #   partial sill: 3.7


    # ----------------
    #  test analysis GAMM ... localised


    # ----------------
    #  test analysis GAMM ... global

    require( mgcv )
    p = spatial.parameters( p=p, type="SSE" )

    O = hydro.db( p=p, DS="bottom.gridded.all"  )
    O = O[, c("plon", "plat", "t", "yr", "dyear")]
    O = O[ which( is.finite(O$t)) ,]


    testyears = c(1950:2007)
    testyears = sort( unique(O$yr) )

    i = which( O$yr %in% testyears )
    # i = sample( i, floor(length(i) * 0.1 ))

    O = O[ i, ]

    e = gam( t ~ s(yr, dyear) + s(plon, plat) , data=O[i,] )
    e = gam( t ~ s(yr, ttime) + s(plon, plat) , data=O )
    e = lme( t ~ yr + ttime  + plon+plat, data=O, random=list(yr=~1 ) )
    e = lme( t ~ yr + ttime  + plon+plat, data=O, random=list(yr=~1,ttime=~1 ) )
    e = gamm( t~ s(yr) + ttime + s(plon,plat), data=O, random=list(plonplat=~1 ) )

    e = gamm( t~ s(yr, dyear) + s(plon,plat), data=O, random=list(plon=~1, plat=~1), correlation=corSpher( c(range, nugg), form = ~ plon+plat, nugget=T) )
    summary(e)
    AIC(e)
    plot(e, all.terms=T, pers=T,theta=60)
  }


  if (0) {
    # data extraction for Katja Fennel
    p = list()
    p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.temperature" )
    p$tyears = c(1950:2012)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p = spatial.parameters( p=p, type= "SSE" )
    out = NULL
    out = temperature.db( p=p, DS="bottom.statistics.annual", yr=p$tyears[1] )
    out = planar2lonlat( out, proj.type=p$internal.projection )
    out = out[ , c("lon", "lat", "tmean" ) ]
    names(out) = c( "lon", "lat", paste("tmean", p$tyears[1],sep="_") )
    for ( y in p$tyears[-1] ) {
      r = NULL
      r = temperature.db( p=p, DS="bottom.statistics.annual", yr=y )
      names(r) = paste( names(r), y, sep="_")
      iname = grep( "tmean", names(r) )
      out = cbind( out, r[,iname] )
    }
    names(out) = c( "lon", "lat", paste( "tmean", p$tyears, sep="_") )
    library(R.matlab)
    writeMat("tdata.mat", tdata=out )
  }





