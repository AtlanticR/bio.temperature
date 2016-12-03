
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

    return(p)
  }

  if (DS=="conker") {

    p$libs = RLibrary( c( p$libs, "conker" ) ) # required for parallel processing

    if (!exists("conker_local_modelengine", p)) {
      message( "'conker_local_modelengine' was not specified, using gam as default")
      p$conker_local_modelengine = "gam" 
    }

    if (p$conker_local_modelengine == "spate" ){

    } else if (p$conker_local_modelengine == "gaussianprocess2Dt") {
      
      message( "NOTE:: The gaussianprocess2Dt method is really slow .. " )

    } else if (p$conker_local_modelengine =="gam") {
      p$conker_local_family = gaussian()
      p$conker_local_modelformula = formula(
        t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") +s(z, k=3, bs="ts")
          + s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts")
          + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )  
      # more than 100 knots and it takes a very long time
      # other possibilities:
        #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ',
        #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
        #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ',
        #     seasonal.smoothed.depth.lonlat.complex = ' s(yr, dyear, bs="ts") + s(yr, k=3, bs="ts") + s(dyear, bs="cc") +s(z, bs="ts") +s(plon, bs="ts") +s(plat, bs="ts") + s(plon, plat, by=tiyr, k=10, bs="ts" ) ',
        #     harmonics.1 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ',
        #     harmonics.2 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' ,
        #     harmonics.3 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) ',
        #     harmonics.1.depth = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) +s(z)  ',
        #     harmonics.1.depth.lonlat = 's(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") +s(z, k=3, bs="ts") +s(plon,k=3, bs="ts") +s(plat, k=3, bs="ts") + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") ',
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

    p$conker_rsquared_threshold = 0.3 # lower threshold
    p$conker_distance_prediction = 7.5 # this is a half window km
    p$conker_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    # using covariates as a first pass essentially makes it kriging with external drift
    p$conker_global_modelengine = "gam"
    p$conker_global_modelformula = formula( t ~ s(z, bs="ts") )
    p$conker_global_family = gaussian()

    p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" )
  
    # other options might work depending upon data density but GP are esp slow .. too slow for bathymetry .. here?
    p$conker_variogram_method = "fast"
  
    p$n.min = p$ny*3 # n.min/n.max changes with resolution
    # min number of data points req before attempting to model timeseries in a localized space
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time

    # if not in one go, then the value must be reconstructed from the correct elements:
    p$conker_sbox = conker_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$conker_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$conker_theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"

    p$conker_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

    return(p)
  }

}

