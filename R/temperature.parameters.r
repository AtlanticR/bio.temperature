
temperature.parameters = function( DS, p=NULL, sp="inla", current.year=NULL ) {

  if (is.null( current.year )) current.year=lubridate::year(lubridate::now())
  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  if (sp=="inla") {
    p$project.root = project.datadirectory( p$project.name )

    p$libs = RLibrary( "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "INLA", "lattice" )
    p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" ) )

    p$tyears = c(1990:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$nw = 10  # numbe rof intervals within a year
    p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required
    return(p)
  }

  # ----------

  if (sp=="gam"){
      p$project.root = project.datadirectory( p$project.name )
      p$libs = RLibrary( c( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "fields" ) )
      p$init.files = bioLibrary(  "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" )

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

}


