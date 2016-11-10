
temperature.parameters = function( p=NULL, current.year=NULL ) {

  if ( is.null( current.year )) current.year=lubridate::year(lubridate::now())
  if ( is.null( p ) ) p=list()

  if ( !exists("project.name", p) ) p$project.name="bio.temperature"
  p$project.root = project.datadirectory( p$project.name )

  p$libs = RLibrary( "lubridate", "gstat", "sp", "rgdal", "parallel", "mgcv", "raster", "fields",
    "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons", "bio.temperature" )

  p$spatial.domain.default = "canada.east"
  p = spacetime_parameters( p=p, type=p$spatial.domain.default )  # default grid and resolution
  p$subregions = c("canada.east", "SSE", "SSE.mpa", "snowcrab" ) # target domains and resolution for additional data subsets .. add here your are of interest

  p$newyear = current.year
  p$tyears = c(1950:current.year)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
  
  p$spacetime_yrs = p$tyears  # yr labels for output

  p$ny = length(p$spacetime_yrs)
  p$nw = 10 # number of intervals in time within a year
  p$nt = p$nw*p$ny # must specify, else assumed = 1

  p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
  p$dyear_centre = p$dyears[ round(p$nw/2) ] + (1/p$nw/2)

  p$spacetime_variogram_engine = "gstat"  # "geoR" seg faults frequently ..
  p$spacetime_rsquared_threshold = 0.3 # lower threshold
  p$spacetime_distance_prediction = 7.5 # this is a half window km
  p$spacetime_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale (dist.max, dist.min)/2 to try in local block search
  
  p$spacetime_engine = "gam" # see model form in spacetime.r (method="xyts")
  p$spacetime_engine_modelformula = formula( 
    t ~ s(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") +s(z, k=3, bs="ts")
      + s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") 
      + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") )  

    if (0) {
      # alternative models .. testing
      p$spacetime_engine = "bayesx"
      p$spacetime_engine_modelformula = formula( 
        t ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps")
          + sx(plon, bs="ps") + sx(plat,  bs="ps") 
          + sx(plon, plat, cos.w, sin.w, yr, bs="te")
      )
      p$bayesx.method="MCMC"

      
    }


  p$spacetime_model_distance_weighted = TRUE

  p$spacetime_covariate_modeltype="gam"
  p$spacetime_covariate_modelformula = formula( t ~ s(z, bs="ts") )

  p$variables = list( Y="t", LOCS=c("plon", "plat"), TIME="tiyr", COV="z" ) 

  p$spacetime_family = gaussian()
  # or to make your own
  # p$spacetime_family = function(offset=0) {
  #   structure(list(
  #     linkfun = function(mu) mu + offset, 
  #     linkinv = function(eta) mu - offset,
  #     mu.eta = function(eta) NA, 
  #     valideta = function(eta) TRUE, 
  #     name = paste0("logexp(", offset, ")") ),
  #     class = "link-glm" )
  # }

  p$dist.max = 75 # length scale (km) of local analysis .. for acceptance into the local analysis/model
  p$dist.min = 2 # lower than this .. subsampling occurs

  p$n.min = p$ny*3 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
  # min number of data points req before attempting to model timeseries in a localized space
  p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time

  # if not in one go, then the value must be reconstructed from the correct elements:
  p$sbbox = spacetime_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
  p$non_convex_hull_alpha = 20  # radius in distance units (km) to use for determining boundaries
  p$theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"

  p$spacetime.noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

  return(p)
}



