
temperature.parameters = function( p=NULL, current.year=NULL, DS="default" ) {

  if ( DS=="default") {
    if ( is.null( current.year )) current.year=lubridate::year(lubridate::now())
    if ( is.null( p ) ) p=list()

    if ( !exists("project.name", p) ) p$project.name="bio.temperature"
    p$project.root = project.datadirectory( p$project.name )

    p$libs = RLibrary( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "raster", "fields",
      "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons", "bio.temperature" )

    p$spatial.domain.default = "canada.east"
    p = spatial_parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution
    p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution for additional data subsets .. add here your are of interest

    p$newyear = current.year
    p$tyears = c(1950:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.

    p$yrs = p$tyears  # yr labels for output
    p$ny = length(p$yrs)
    p$nw = 10 # number of intervals in time within a year
    p$nt = p$nw*p$ny # must specify, else assumed = 1
    p$tres = 1/ p$nw # time resolution

    tout = expand.grid( yr=p$tyears, dyear=1:p$nw, KEEP.OUT.ATTRS=FALSE )
    tout$tiyr = tout$yr + tout$dyear/p$nw - p$tres/2 # mid-points
    tout = tout[ order(tout$tiyr), ]
    p$ts = tout$tiyr

    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
    p$dyear_centre = p$dyears[ round(p$nw/2) ] + p$tres/2

    p$spatial_distance_max = 25 # obsolete? .. used by old inverse distance method
    p$phi # FFT -based methods range parameter .. == covariance range parameter

    return(p)
  }
  

  if (DS=="hivemod") {

    p$libs = RLibrary( c( p$libs, "hivemod" ) ) # required for parallel processing
    p$storage.backend="bigmemory.ram"
    p$boundary = TRUE 
    p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land
 
    p$hivemod_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$hivemod_lowpass_phi = p$pres / 5 # FFT-baed methods cov range parameter
    p$hivemod_lowpass_nu = 0.5
    p$hivemod_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$hivemod_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$hivemod_rsquared_threshold = 0.25 # lower threshold
    p$hivemod_distance_prediction = 7.5 # this is a half window km
    p$hivemod_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$hivemod_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$hivemod_distance_min = 1
    p$hivemod_distance_max = 50 

    # other options might work depending upon data density but GP are esp slow .. too slow for bathymetry .. here?
    p$hivemod_variogram_method = "fast"
  
    p$n.min = 150 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    p$n.max = 1000 # numerical time/memory constraint -- anything larger takes too much time
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # 


    p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" )
 
    # using covariates as a first pass essentially makes it kriging with external drift
    p$hivemod_global_modelengine = NULL #"gam"
    p$hivemod_global_modelformula = NULL # formula( t ~ s(z, bs="ts") ) # marginally useful .. consider removing it.
    p$hivemod_global_family = gaussian()

    if (!exists("hivemod_local_modelengine", p)) p$hivemod_local_modelengine = "gam" # "twostep" 
  
  
    if (p$hivemod_local_modelengine == "gaussianprocess2Dt") {
 
      message( "NOTE:: The gaussianprocess2Dt method is really slow .. " )
 
    } else if (p$hivemod_local_modelengine =="gam") {
      # 32 hours on nyx all cpus
      p$hivemod_local_family = gaussian()
      p$hivemod_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s( log(z), k=3, bs="ts")
          + s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts")
          + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )  
      # more than 100 knots and it takes a very long time
      # other possibilities:
        #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ',
        #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
        #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ',
        p$hivemod_local_model_distanceweighted = TRUE
 
    } else if (p$hivemod_local_modelengine =="twostep") {
      # 34 hr with 8 CPU RAM on thoth, using 48 GB RAM .. about 1/3 faster than 24 cpus systems
      # 42 hrs on tartarus all cpus 
      # 18 GB RAM for 24 CPU .. 
      p$hivemod_local_family = gaussian()
      p$hivemod_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s(log(z), k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts") )
        # similar to GAM model but no spatial component .. space is handled via FFT
      p$hivemod_local_model_distanceweighted = TRUE
      p$hivemod_fft_filter = "spatial.process"

    } else if (p$hivemod_local_modelengine =="spate") {
 
      # similar to the two-step but use "spate" (spde, bayesian, mcmc) instead of "fields" (GMRF, ML)
      p$hivemod_local_family = gaussian()
      p$hivemod_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s( log(z), k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts") )
        # similar to GAM model but no spatial component , space and time are handled via FFT but time is seeded by the averge local TS signal (to avoid missing data isses in time.)
      p$hivemod_local_model_distanceweighted = TRUE
 
    } else if (p$hivemod_local_modelengine == "bayesx") {
 
      # bayesx families are specified as characters, this forces it to pass as is and 
      # then the next does the transformation internal to the "hivemod__bayesx"
      p$hivemod_local_family = gaussian() 
      p$hivemod_local_family_bayesx = "gaussian" 

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$hivemod_local_modelformula = formula(
        t ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps")
          + sx(plon, bs="ps") + sx(plat,  bs="ps")
          + sx(plon, plat, cos.w, sin.w, yr, bs="te")  # te is tensor spline
      )
      p$hivemod_local_model_bayesxmethod="MCMC"
      p$hivemod_local_model_distanceweighted = FALSE
    
    } else {
    
      message( "The specified hivemod_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }
    

    return(p)
  }

}

