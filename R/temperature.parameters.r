

temperature.parameters = function( p=NULL, year.assessment=NULL, DS="default" ) {

  if ( DS=="default") {
    if ( is.null( p ) ) p=list()

    p$libs = RLibrary( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "raster", "fields",
      "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons", "bio.temperature" )

    if ( !exists("project.name", p) ) p$project.name="bio.temperature"
    if ( !exists("project.root", p) ) p$project.root = project.datadirectory( p$project.name )

    if ( !exists("spatial.domain", p) ) p$spatial.domain = "canada.east" # canada.east.highres and canada.east.superhighres result in memory overflow
    p = spatial_parameters( p=p, type=p$spatial.domain )  # default grid and resolution

    if ( !exists("spatial.domain.subareas", p) )  p$spatial.domain.subareas = c( "SSE.mpa", "SSE", "snowcrab" ) # target domains and resolution for additional data subsets .. add here your are of interest

    if ( is.null( year.assessment )) year.assessment=lubridate::year(lubridate::now())
    p$newyear = year.assessment
    p$tyears = c(1950:year.assessment)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$tyears.climatology = p$tyears

    p$bstats = c("tmean", "tsd", "tmin", "tmax", "amplitude", "degreedays" )

    if ( !exists("yrs", p) )  p$yrs = p$tyears  # yr labels for output

    p$ny = length(p$yrs)
    p$nw = 10 # number of intervals in time within a year
    p$nt = p$nw*p$ny # must specify, else assumed = 1
    p$tres = 1/ p$nw # time resolution

    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
    p$dyear_centre = p$dyears[ round(p$nw/2) ] + p$tres/2

    p$prediction.dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01")) # used for creating timeslices .. needs to match the values in indicators.parameters()

    # output timeslices for predictions
    tout = expand.grid( yr=p$yrs, dyear=1:p$nw, KEEP.OUT.ATTRS=FALSE )
    tout$tiyr = tout$yr + tout$dyear/p$nw - p$tres/2 # mid-points
    tout = tout[ order(tout$tiyr), ]
    p$prediction.ts = tout$tiyr   # predictions at these time values (decimal-year)

    return(p)
  }


  if (DS=="lbm") {

    p$libs = RLibrary( c( p$libs, "lbm" ) ) # required for parallel processing
    if (!exists("storage.backend", p))  p$storage.backend="bigmemory.ram"
    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )

    if (!exists("boundary", p)) p$boundary = FALSE
    if (!exists("depth.filter", p)) p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)
    if (!exists("lbm_eps", p)) p$lbm_eps = 0.1  # distance units for eps noise to permit mesh gen for boundaries
    if (!exists("lbm_quantile_bounds", p)) p$lbm_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations

    if (!exists("lbm_rsquared_threshold", p)) p$lbm_rsquared_threshold = 0.25 # lower threshold

    if (!exists("lbm_distance_prediction", p)) p$lbm_distance_prediction = 4 # this is a half window km
    if (!exists("lbm_distance_statsgrid", p)) p$lbm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    if (!exists("lbm_distance_scale", p)) p$lbm_distance_scale = 25 # km ... approx guess of 95% AC range
    if (!exists("lbm_distance_min", p)) p$lbm_distance_min = p$lbm_distance_statsgrid
    if (!exists("lbm_distance_max", p)) p$lbm_distance_max = 65

    if (!exists("n.min", p)) p$n.min = 500 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("n.max", p)) p$n.max = 5000 # numerical time/memory constraint -- anything larger takes too much time
    if (!exists("sampling", p)) p$sampling = c( 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25 )  # mostly used to down sample when there is too much data (depth, substrate)

    # fast variogram options
    if (!exists("lbm_variogram_method", p)) p$lbm_variogram_method = "fast"

    # global model options
    # using covariates as a first pass essentially makes it ~ kriging with external drift
    p$lbm_global_modelengine = NULL # no covars .. alternatives .. "gam"
    p$lbm_global_modelformula = NULL
    # p$lbm_global_modelformula = formula( t ~ s(z, bs="ts" + s(s.range, bs="ts") + s(dZ, bs="ts") + s(ddZ, bs="ts") + s(log.substrate.grainsize, bs="ts")  ) ) # marginally useful .. consider removing it.
    if (!exists("lbm_global_family", p))  p$lbm_global_family = gaussian()




    # local model options    
    # if (!exists("lbm_local_modelengine", p)) p$lbm_local_modelengine = "spate" # "twostep" might be interesting to follow up
    if (!exists("lbm_local_modelengine", p)) p$lbm_local_modelengine = "gam" # "twostep" might be interesting to follow up
    if (!exists("lbm_local_modelformula", p)) p$lbm_local_modelformula = formula(
      t ~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")
        + s(plon, k=10, bs="ts") + s(plat, k=10, bs="ts")
        + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )
    # more than 100 knots and it takes a very long time, 50 seems sufficient, given the large-scaled pattern outside of the prediction box
    # other possibilities:
      #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ',
      #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
      #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ',
    if (!exists("lbm_local_family", p)) p$lbm_local_family = gaussian()
    if (!exists("lbm_local_model_distanceweighted", p)) p$lbm_local_model_distanceweighted = TRUE


    # local engine-specific options
    if (p$lbm_local_modelengine =="gam") {
      # 32 hours on nyx all cpus;
      # XX hrs on thoth all cpus
       # p$lbm_gam_optimizer="perf"
        p$lbm_gam_optimizer=c("outer", "bfgs")

    # } else if (p$lbm_local_modelengine =="fft") {
    #   # lowpass seems a bit too noisy
    #   # spatial process and lowpass_spatial.process are over-smooth
    #   # p$lbm_fft_filter = "lowpass" # only act as a low pass filter .. depth has enough data for this. Otherwise, use:
    #   p$lbm_fft_filter = "spatial.process" # to ~ unviersal krige with external drift

    # } else if (p$lbm_local_modelengine == "gaussianprocess2Dt") {
    #   message( "NOTE:: The gaussianprocess2Dt method is really slow .. " )
    # }
    } else if (p$lbm_local_modelengine =="twostep") {
      # 34 hr with 8 CPU RAM on thoth, using 48 GB RAM .. about 1/3 faster than 24 cpus systems
      # 42 hrs on tartarus all cpus
      # 18 GB RAM for 24 CPU .. 
      p$lbm_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")
          + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts")
          + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )
        # similar to GAM model but no spatial component .. space is handled via FFT
  
      # p$lbm_twostep_space = "spatial.process"
      p$lbm_twostep_space = "krige"

    } else if (p$lbm_local_modelengine =="spate") {
      p$lbm_spate_method = "mcmc_fast"  # "mcmc" is too slow for temperature
      p$lbm_distance_prediction = 4 # this is a half window km
      p$lbm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      p$lbm_spate_boost_timeseries = TRUE  # use simple GAM spectral contraint to structure timeseries as spate's fft in time seems to cause overfitting ?
      p$lbm_spate_nburnin = 1000
      p$lbm_spate_nposteriors = 1000
      p$lbm_spate_nCovUpdates = 20 # no of times to update cov matrix during simulations

    } else if (p$lbm_local_modelengine == "bayesx") {

      # bayesx families are specified as characters, this forces it to pass as is and
      # then the next does the transformation internal to the "lbm__bayesx"
      p$lbm_local_family = "gaussian"

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$lbm_local_modelformula = formula(
        t ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps")
          + sx(plon, bs="ps") + sx(plat,  bs="ps")
          + sx(plon, plat, cos.w, sin.w, yr, bs="te")  # te is tensor spline
      )
      p$lbm_local_model_bayesxmethod="MCMC"
      p$lbm_local_model_distanceweighted = FALSE

    } else {

      message( "The specified lbm_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }

    # for fft-based methods that require lowpass:


    p$lbm_lowpass_phi = p$pres / 5 # FFT-baed methods cov range parameter .. not required for "spatial.process" ..
    p$lbm_lowpass_nu = 0.5

    p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" )

    p$varnames = c( p$variables$LOCS, p$variables$COV ) # to extract for prediction

    return(p)
  }

}

