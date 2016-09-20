
temperature.parameters = function( DS, p=NULL, current.year=NULL ) {

  if (is.null( current.year )) current.year=lubridate::year(lubridate::now())
  if ( is.null(p) ) p=list()

  if ( !exists("project.name", p) ) p$project.name="bio.temperature"
  p$project.root = project.datadirectory( p$project.name )

  if (DS=="inla") {

    p$libs = RLibrary( "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "INLA", "lattice" )
    p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" ) )

    p$tyears = c(1990:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$nw = 10  # number of intervals within a year
    p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required
    return(p)
  }

  # ----------

  if (DS=="gam"){
      p$libs = RLibrary( c( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "fields" ) )
      p$libs = unique( c( p$libs, bioLibrary(  "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" ) ) )

      # p$tyears = c(1910:2013)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
      p$tyears = c(1950:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
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
      p$newyear = current.year
      p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution
      p$spatial.domain.default = "canada.east"
      p = spatial.parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution
    return(p)
  }


  if (DS=="bio.bathymetry.spacetime") {
   
    p$libs = RLibrary( c( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "ff", "ffbase", "fields" ) )
    p$libs = unique( c( p$libs, bioLibrary(  "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" ) ) )
  
    p$rootdir = file.path( p$project.root, "spacetime" )

    p$fn.P =  file.path( p$rootdir, paste( "spacetime", "predictions", p$spatial.domain, "rdata", sep=".") )
    p$fn.S =  file.path( p$rootdir, paste( "spacetime", "statistics", p$spatial.domain, "rdata", sep=".") )
    p$fn.results.covar =  file.path( p$rootdir, paste( "spatial", "covariance", p$spatial.domain, "rdata", sep=".") )
    p$variogram.engine = "gstat"  # "geoR" seg faults frequently ..
    p$dist.mwin = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$upsampling = c( 1.1, 1.2, 1.5, 2 )  # local block search fractions
    p$downsampling = c( 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2 ) # local block search fractions  -- need to adjust based upon data density
    p$mesh.boundary.resolution = 150
    p$mesh.boundary.convex = -0.025

    p$variables = list( Y="t", LOCS=c("plon", "plat") )
    p$spacetime.link = function( X ) { log(X + 1000) }  ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
    p$spacetime.invlink = function( X ) { exp(X) - 1000 }
    p$dist.max = 100 # length scale (km) of local analysis .. for acceptance into the local analysis/model
    p$dist.min = 75 # lower than this .. subsampling occurs
    p$dist.pred = 0.95 # % of dist.max where **predictions** are retained (to remove edge effects)
    p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time
    p$expected.range = 50 #+units=km km , with dependent var on log scale
    p$expected.sigma = 1e-1  # spatial standard deviation (partial sill) .. on log scale
    p$spatial.field.name = "spatial.field"  # name used in formula to index the spatal random field
    p$modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covariance model .. defined in spacetime.interpolate.inla.local (below)
    p$spacetime.family = "gaussian"
    p$spacetime.outputs = c( "predictions.projected", "statistics" ) # "random.field", etc.
    
    # if not in one go, then the value must be reconstructed from the correct elements:
    p$sbbox = spacetime.db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$spacetime.stats.boundary.redo = FALSE ## estimate boundart of data to speed up stats collection? Do not need to redo if bounds have already been determined
    
    p$nPreds = p$nplons * p$nplats
    # p$tyears = c(1910:2013)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$tyears = c(1950:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$ny = length(p$tyears)
    p$nw = 10 # number of intervals in time within a year
    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
    p$gam.optimizer = "nlm" ## other optimizers:: "bam" (slow), "perf"(ok), "nlm" (good), "bfgs" (ok), "newton" (default)
    p$nMin.tbot = p$ny*3 # min number of data points req before attempting to model timeseries in a localized space
    p$dist.km = c( 2.5, 5, 7.5, 10, 12.5, 15 ) # "manhattan" distances to extend search for data
    p$maxdist = 20 # if using gstat  max dist to interpolate in space
    # possible methods: "annual", "seasonal.basic", "seasonal.smoothed", "harmonics.1", "harmonics.2", "harmonics.3", "inla.ts.simple"
    p$tsmethod = "harmonics.1"
    # p$spmethod = "inverse.distance"  ## too slow
    # p$spmethod = "gam" ## too smooth
    p$spmethod = "kernel.density" ## best
    p$theta = 5 # FFT kernel bandwidth (SD of kernel) for method p$spmethod = "kernel.density"
    p$nsd = 6 # number of SD distances to pad boundaries with 0 for FFT  in method  p$spmethod = "kernel.density
    p$newyear = current.year
    p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution
    p$spatial.domain.default = "canada.east"
    p = spatial.parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution
    return(p)
  }


}


