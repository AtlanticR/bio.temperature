
temperature.parameters = function( p=NULL, current.year=NULL, DS="default" ) {

  if ( DS=="default") {
    if ( is.null( current.year )) current.year=lubridate::year(lubridate::now())
    if ( is.null( p ) ) p=list()

    p$libs = RLibrary( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "raster", "fields",
      "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons", "bio.temperature" )

    if ( !exists("project.name", p) ) p$project.name="bio.temperature"
    if ( !exists("project.root", p) ) p$project.root = project.datadirectory( p$project.name )

    if ( !exists("spatial.domain.default", p) ) p$spatial.domain.default = "canada.east"
    p = spatial_parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution
    
    if ( !exists("subregions", p) )  p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution for additional data subsets .. add here your are of interest

    p$newyear = current.year
    p$tyears = c(1950:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.

    if ( !exists("yrs", p) )  p$yrs = p$tyears  # yr labels for output
    
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

    p$spatial_distance_max = 25 # obsolete .. used by old inverse distance method .. to be retired shortly

    return(p)
  }
  

  if (DS=="lbm") {

    p$libs = RLibrary( c( p$libs, "lbm" ) ) # required for parallel processing
    p$clusters = rep("localhost", detectCores() )    
    p$storage.backend="bigmemory.ram"

    p$boundary = TRUE 
    p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)
    p$lbm_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$lbm_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$lbm_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$lbm_rsquared_threshold = 0.25 # lower threshold
    p$lbm_distance_prediction = 7.5 # this is a half window km
    p$lbm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$lbm_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$lbm_distance_min = p$lbm_distance_statsgrid 
    p$lbm_distance_max = 50 

  
    p$n.min = 200 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    p$n.max = 2000 # numerical time/memory constraint -- anything larger takes too much time
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # 


    p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" )
    
    if (!exists("lbm_variogram_method", p)) p$lbm_variogram_method = "fast"
    if (!exists("lbm_local_modelengine", p)) p$lbm_local_modelengine = "gam" # "twostep" might be interesting to follow up

    # using covariates as a first pass essentially makes it ~ kriging with external drift
    p$lbm_global_modelengine = NULL #"gam"
    p$lbm_global_modelformula = NULL # formula( t ~ s(z, bs="ts") ) # marginally useful .. consider removing it.
    p$lbm_global_family = gaussian()
  
    p$lbm_local_family = gaussian()

    if (p$lbm_local_modelengine == "gaussianprocess2Dt") {
 
      message( "NOTE:: The gaussianprocess2Dt method is really slow .. " )
 
    } else if (p$lbm_local_modelengine =="gam") {
      # 32 hours on nyx all cpus; 
      # XX hrs on thoth all cpus
      
      p$lbm_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") + s( log(z), k=3, bs="ts")
          + s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts")
          + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )  
      # more than 100 knots and it takes a very long time, 50 seems sufficient, given the large-scaled pattern outside of the prediction box
      # other possibilities:
        #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ',
        #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
        #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ',
        p$lbm_local_model_distanceweighted = TRUE
        p$lbm_gam_optimizer="perf"
        # p$lbm_gam_optimizer=c("outer", "bfgs") 
    } else if (p$lbm_local_modelengine =="fft") {

      # p$lbm_fft_filter = "lowpass" # only act as a low pass filter .. depth has enough data for this. Otherwise, use: 
      p$lbm_fft_filter = "spatial.process" # to ~ unviersal krige with external drift
      p$lbm_lowpass_phi = p$pres / 5 # FFT-baed methods cov range parameter
      p$lbm_lowpass_nu = 0.5

    } else if (p$lbm_local_modelengine =="twostep") {
      # 34 hr with 8 CPU RAM on thoth, using 48 GB RAM .. about 1/3 faster than 24 cpus systems
      # 42 hrs on tartarus all cpus 
      # 18 GB RAM for 24 CPU .. 
      p$lbm_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") + s(log(z), k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts", k=9) )
        # similar to GAM model but no spatial component .. space is handled via FFT
      p$lbm_local_model_distanceweighted = TRUE

      p$lbm_fft_filter = "spatial.process"
      p$lbm_lowpass_phi = p$pres / 5 # FFT-baed methods cov range parameter .. not required for "spatial.process" ..
      p$lbm_lowpass_nu = 0.5

    } else if (p$lbm_local_modelengine =="spate") {
 
      # similar to the two-step but use "spate" (spde, bayesian, mcmc) instead of "fields" (GMRF, ML)
      p$lbm_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") + s( log(z), k=3, bs="ts")
          + s(cos.w, sin.w, yr, bs="ts") )
        # similar to GAM model but no spatial component , space and time are handled via FFT but time is seeded by the averge local TS signal (to avoid missing data isses in time.)
      p$lbm_local_model_distanceweighted = TRUE
 
    } else if (p$lbm_local_modelengine == "bayesx") {
 
      # bayesx families are specified as characters, this forces it to pass as is and 
      # then the next does the transformation internal to the "lbm__bayesx"
      p$lbm_local_family_bayesx = "gaussian" 

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$lbm_local_modelformula = formula(
        t ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps")
          + sx(plon, bs="ps") + sx(plat,  bs="ps")
          + sx(plon, plat, cos.w, sin.w, yr, bs="te")  # te is tensor spline
      )
      p$lbm_local_model_bayesxmethod="MCMC"
      p$lbm_local_model_distanceweighted = FALSE
    
    } else {
    
      message( "The specified lbm_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }
    

    return(p)
  }

}

