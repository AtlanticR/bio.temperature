
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
    p$theta # FFT kernel bandwidth (SD of kernel if gaussian) required for method "harmonic.1/kernel.density, etc .."

    return(p)
  }
  

  if (DS=="conker") {

    p$libs = RLibrary( c( p$libs, "conker" ) ) # required for parallel processing
    p$storage.backend="bigmemory.ram"
 
    if (!exists("conker_local_modelengine", p)) {
      message( "'conker_local_modelengine' was not specified, using gam as default")
      p$conker_local_modelengine = "gam" # default is gam ..
      p$conker_local_modelengine = "twostep" # testing twostep ~ hybrid of gam + kernel density ..
    }

 
    if (p$conker_local_modelengine == "gaussianprocess2Dt") {
      message( "NOTE:: The gaussianprocess2Dt method is really slow .. " )
    } else if (p$conker_local_modelengine =="gam") {
      p$conker_local_family = gaussian()
      p$conker_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s(z, k=3, bs="ts")
          + s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts")
          + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )  
      # more than 100 knots and it takes a very long time
      # other possibilities:
        #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ',
        #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
        #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ',
        p$conker_local_model_distanceweighted = TRUE
    } else if (p$conker_local_modelengine =="twostep") {
      # 18 GB RAM for 24 CPU .. 
      # 34 hr with 8 CPU RAM on thoth, using 48 GB RAM
      p$conker_local_family = gaussian()
      p$conker_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s(z, k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts") )
        # similar to GAM model but no spatial component .. space is handled via FFT
      p$conker_local_model_distanceweighted = TRUE
    } else if (p$conker_local_modelengine =="spate") {
      # similar to the two-step but use "spate" (spde, bayesian, mcmc) instead of "fields" (GMRF, ML)
      p$conker_local_family = gaussian()
      p$conker_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s(z, k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts") )
        # similar to GAM model but no spatial component , space and time are handled via FFT but time is seeded by the averge local TS signal (to avoid missing data isses in time.)
      p$conker_local_model_distanceweighted = TRUE
    } else if (p$conker_local_modelengine == "bayesx") {
      # bayesx families are specified as characters, this forces it to pass as is and 
      # then the next does the transformation internal to the "conker__bayesx"
      p$conker_local_family = gaussian() 
      p$conker_local_family_bayesx = "gaussian" 

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$conker_local_modelformula = formula(
        t ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps")
          + sx(plon, bs="ps") + sx(plat,  bs="ps")
          + sx(plon, plat, cos.w, sin.w, yr, bs="te")  # te is tensor spline
      )
      p$conker_local_model_bayesxmethod="MCMC"
      p$conker_local_model_distanceweighted = FALSE
    
    } else {
    
      message( "The specified conker_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }
    

    # using covariates as a first pass essentially makes it kriging with external drift
    p$conker_global_modelengine = NULL #"gam"
    p$conker_global_modelformula = NULL # formula( t ~ s(z, bs="ts") ) # marginally useful .. consider removing it.
    p$conker_global_family = gaussian()

    p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" )
  
    p$conker_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$conker_distance_prediction = 7.5 # this is a half window km
    p$conker_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$conker_distance_min = 1
    p$conker_distance_max = 100 

    # other options might work depending upon data density but GP are esp slow .. too slow for bathymetry .. here?
    p$conker_variogram_method = "fast"
  
    p$n.min = 200 # n.min/n.max changes with resolution
    # min number of data points req before attempting to model timeseries in a localized space
    p$n.max = 1000 # numerical time/memory constraint -- anything larger takes too much time

    p$conker_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$conker_theta = p$theta # FFT kernel bandwidth (SD of kernel if gaussian) required for method "harmonic.1/kernel.density, etc .."

    p$conker_rsquared_threshold = 0.25 # lower threshold
    p$conker_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$conker_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations

    p$boundary = TRUE 
    p$depth.filter = log(0.5) # depth is given as log(depth) so, choose andy stats locations with elevation > 0.5 m as being on land

    return(p)
  }

}

