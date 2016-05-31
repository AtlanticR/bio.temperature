
  p = list( project.name = "temperature" )
  p$project.root = project.datadirectory( p$project.name )

  p$libs = RLibrary( "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "INLA", "lattice" )
  p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons" , "bio.temperature" ) )


  p$tyears = c(1990:1995)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
  p$nw = 10  # numbe rof intervals within a year
  p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required

  # p = temperature.db( p=p, DS="bigmemory.filenames" )


  # wrap up
  p = temperature.db( p=p, DS="bigmemory.cleanup" )


  # test timeseries methods and spatial, spatio-temporal methods on ocean temperature data
  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors


    # -------------------
    # 1. example model fit of a simple AR1 data series with inla
    ar.coef = 0.8
    n = 200
    ti =  1:n

    y0 = arima.sim(n=n, model=list( ar=ar.coef ))
    missing.values = sample( ti, trunc(n * 0.3))
    y1 = y0
    y1[ missing.values ]  = NA

    y1  = y1 + rnorm( n, mean=0, sd=sd( y0 ) )


    dat = data.frame( ti=ti, ti2=ti, y0=y0, y1=y1 )
    r = inla( y1 ~ 0 + f( ti, model="ar1" )  ,
      data = dat
    )
    dat$predictions =  r$summary.random$ti[["mean"]]
    ar.pred =  r$summary.hyperpar["Rho for ti", "mean" ]
    mm = glm( predictions~y0, data=dat )

    graphics.off()
    par(mfrow = c(3, 1) )
    plot(y0 ~ ti, dat, type="l", col="grey" )  # original data
    points(y1 ~ ti, dat, type="p", col="grey" )  # tainted data
    legend( "topright", legend=( c( paste( "Ar=", ar.coef ), paste( "Ar estimated=", round(ar.pred, 3) ) ) ) )

    plot( y0 ~ ti, dat, type="n" )  # original data scales
    points( y1 ~ ti, dat, type="p", col="grey" )  # tainted data
    lines( y0 ~ ti, dat, col="grey" )
    lines( predictions ~ ti, dat, col="green" )

    plot(predictions ~ y0, dat)
    abline( reg=mm, col="green" )
    mmc = round(coef(mm),3)
    legend( "topleft", legend= c(
      paste("b =", mmc[1]),
      paste("m =", mmc[2]),
      paste("R2 = ",  round( Rsquared(mm),3) ),
      paste("AIC = ", round( AIC(mm),3) )
    ) )


    u0 = acf( dat$y0 )
    up = acf( dat$predictions )
    u = merge(  cbind(acfy0 =u0$acf, lag=u0$lag), cbind( acfpred= up$acf, lag=up$lag), by="lag" )


    x11()
    par(mfrow = c(3, 1) )
    plot(u0)
    plot(up)

    plot( acfpred ~ lag, u, type="l" , col="red" )
    lines( acfy0 ~ lag, u, type="l" , col="gray" )



    # -------------------
    # 2. now temperature data timeseries

    # resource for kernel smoothers, spectral analysis, etc
    # https://onlinecourses.science.psu.edu/stat510/?q=book/export/html/57

    p = list()
    p$libs = RLibrary( "lubridate", "gstat", "sp", "parallel", "mgcv", "bigmemory", "INLA", "lattice" )
    p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.temperature", "bio.polygons" ) )

    p$tyears = c(1970:2013)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$nw = 10
    p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    t0 = hydro.db( p=p, DS="bottom.gridded.all"  )
    t0 = t0[ which( t0$yr %in% p$tyears ), ]
  # t0 = t0[ filter.region.polygon( t0, region="isobath1000m" ) , ]

   # scotianshelf = locator()
    scotianshelf = polygon.db( id="scotia.fundy.with.buffer"  )
    plot( t0$lat ~ t0$lon, pch="." )
    lines( scotianshelf )

    a = which( point.in.polygon( t0$lon, t0$lat, scotianshelf$lon, scotianshelf$lat ) != 0 )

    t0 = t0[a,]
 --- add prediction locations /times / depths / etc ....

    fp = data.frame( lon= -63, lat= 45 )
    dl = 0.5
    subarea = which( t0$lon > (fp$lon-dl) & t0$lon < (fp$lon+dl) & t0$lat < (fp$lat+dl) &  t0$lat > (fp$lat-dl) )
    t0 = t0[ subarea ,]

    t0$yrindex = t0$yr - min(t0$yr) + 1
    t0$pryr = t0$dyear
    t0$ti = t0$yr + t0$pryr
    t0$logZ = log(t0$z)
    t0$tC = t0$t
    t0$b0 = 1  # intercepts

    R <- inla(
        tC ~ 0 + b0
          + f( z, model='rw2')
          + f( yr, model='ar1', param=c(1,0.01))
          + f( pryr, model='ar1', cyclic=TRUE, param=c(1,0.01) ),
        family='gaussian',
        data=t0,
        control.compute=list(dic=TRUE),
        control.predictor=list( compute=TRUE),
        verbose=TRUE
    )
    summary(R)
    plot(R)

    plot( t ~ ti, t0  )
    lines( R$summary.fitted.values$mean ~ t0$ti, col="green" )

    round( R$summary.fixed)

Fixed effects:
     mean     sd 0.025quant 0.5quant 0.975quant   mode   kld
b0 3.6505 2.9344       -2.1   3.6705     9.2672 3.6834 1e-04

Random effects:
Name	  Model
 z   IID model
yr   AR1 model
pryr   AR1 model

Model hyperparameters:
                                        mean   sd     0.025quant 0.5quant 0.975quant mode
Precision for the Gaussian observations 0.3308 0.0363 0.2617     0.3305   0.4037     0.3314
Precision for z                         0.5211 0.1061 0.3567     0.5054   0.7699     0.4730
Precision for yr                        0.5285 0.1594 0.2723     0.5107   0.8919     0.4774
Rho for yr                              0.5352 0.1511 0.2112     0.5457   0.7941     0.5649
Precision for pryr                      0.1365 0.0669 0.0376     0.1271   0.2914     0.0981
Rho for pryr                            0.9212 0.0393 0.8309     0.9266   0.9795     0.9440

Expected number of effective parameters(std dev): 410.20(53.28)
Number of equivalent replicates : 2.606

Deviance Information Criterion: 4617.06
Effective number of parameters: 412.15

Marginal Likelihood:  -2485.71
Posterior marginals for linear predictor and fitted values computed


    p$gam.optimizer = "nlm"
    p$tsmethod = "seasonal.basic"
    p$tsmethod = "seasonal.smoothed"
    p$tsmethod = "harmonics.1"


    # choose model formula for GAM-based models
    mf = switch( p$tsmethod ,
      annual = ' t ~ s(yr) ',
      seasonal.basic = ' t ~ s(yr) + s(dyear, bs="cc") ',
      seasonal.smoothed = ' t ~ s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
      harmonics.1 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ',
      harmonics.2 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' ,
      harmonics.3 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) '
    )
    mf = formula(mf)

   t0$w = 1 / (( t0$lon - fp$lon)**2 + (t0$lat - fp$lat)**2 )# weight data in space: inverse distance squared
   t0$w[ which( t0$w < 1e-3 ) ] = 1e-3
   # t0 = t0[,c("t", "w", "yr", "dyear" )]
   # data transformations and creation of new variables where required for raw data

   x = t0

   # default output grid
   dyears = (c(1:(p$nw))-1)  / p$nw # intervals
   OP = expand.grid( dyear=dyears, yr=p$tyears )
   OP$fit = NA  # these will be filled in with predicted fits and se's
   OP$se  = NA


    OP$tiyr = OP$yr + OP$dyears
    x$tiyr =  x$yr + x$dyears

   if ( p$tsmethod %in% c( "harmonics.1", "harmonics.2", "harmonics.3"  ) ) {
      x$cos.w  = cos( x$tiyr )
      x$sin.w  = sin( x$tiyr )
      OP$cos.w  = cos( OP$tiyr )
      OP$sin.w  = sin( OP$tiyr )
      # compute additional harmonics only if required (to try to speed things up a bit)
      if ( p$tsmethod %in% c( "harmonics.2", "harmonics.3"  ) ) {
        x$cos.w2 = cos( 2*x$tiyr )
        x$sin.w2 = sin( 2*x$tiyr )
        OP$cos.w2 = cos( 2*OP$tiyr )
        OP$sin.w2 = sin( 2*OP$tiyr )
      }
      if ( p$tsmethod %in% c( "harmonics.3"  ) ) {
        x$cos.w3 = cos( 3*x$tiyr )
        x$sin.w3 = sin( 3*x$tiyr )
        OP$cos.w3 = cos( 3*OP$tiyr )
        OP$sin.w3 = sin( 3*OP$tiyr )
      }
    }

    # estimate model parameters
    tsmodel = NULL
    tsmodel = switch( p$gam.optimizer ,
      bam = try( bam( mf, data=x, weights=w ) ) ,
      bfgs = try( gam( mf, data=x, weights=w, optimizer=c("outer","bfgs")  ) ) ,
      perf = try( gam( mf, data=x, weights=w, optimizer=c("perf")  ) ) ,
      newton = try( gam( mf, data=x, weights=w, optimizer=c("outer","newton")  ) ) ,
      nlm = try( gam( mf, data=x, weights=w, optimizer=c("outer","nlm")  ) )
    )

    summary(tsmodel)
    AIC(tsmodel)

    if ( ! "try-error" %in% class(tsmodel) ) {
      out = try( predict( tsmodel, newdata=OP, type="response", se.fit=T ) )
      OP$fit = out$fit
      OP$se = out$se
      x11()
      plot( OP$fit ~ OP$tiyr,  type="l", main=p$tsmethod, sub=paste( "AIC", round(AIC(tsmodel))) )
    }

    # other methods
    OP$time = OP$yr + OP$dyears
    OP = OP[ order( OP$time ) ,]

    # STL method
    vn = "fit"
    OPts = ts( OP[,vn] , start=c( min(p$tyears), 1), frequency=12 )
    plot.default(  OPts, type="l" )

    OPstl = stl( OPts, s.window=4) # seasonal decomposition using loess
    plot(OPstl)

    # StructTS method
    OPstr  = StructTS( OPts, type = "BSM") ### much slower ... used only as a daignostic tool checking stability
    OPstrsm = tsSmooth( OPstr )
    plot( OPstrsm )

    spectrum( OPts )
    kn = kernel( "modified.daniell", c(2,1))

    vn = "fit"
    opop = kernapply( as.vector(OP[,vn] ), kn)
    plot(OP[,vn], type="l")
    lines(opop, type="l", col="orange")

    ppp = kernapply( OPts, kn)
    plot.default(  OPts, type="l", col="green" )
    lines( ppp, col = "red")

    ooo = spectrum( OPts, kn, plot=FALSE )
    plot( ooo, plot.type = "marginal") # the default type
    plot( ooo, plot.type = "coherency")
    plot( ooo, plot.type = "phase")


    require(forecast)
    Arima( OPts )



  # -------------------
  # 3. spatial model only on temperature

    p = list()
    p$libs = RLibrary( c( "chron", "gstat", "sp", "parallel", "mgcv", "bigmemory", "INLA", "lattice" ) )
    p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.temperature", "bio.polygons" ) )

    p$tyears = c(2012)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
    p$nw = 10
    p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    t0 = hydro.db( p=p, DS="bottom.gridded.all"  )

    # scotianshelf = locator()
    scotianshelf = read.table( find.bio.gis( "scotia.fundy.with.buffer.dat"  ) )
    names( scotianshelf) = c("lon", "lat")


    plot( t0$lat ~ t0$lon, pch="." )
    lines( scotianshelf )

    a = which( point.in.polygon( t0$lon, t0$lat, scotianshelf$lon, scotianshelf$lat ) != 0 )
    t0 = t0[a,]
  # boundary domain

    locs0  = as.matrix( t0[,c("plon", "plat")] )

  M0.domain = inla.nonconvex.hull( locs0, convex=10, resolution=120 )


  M0 = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary=M0.domain,
      max.edge=c( 20, 120 ),  # max size of a triange (in, out)
      cutoff=15 # min distance allowed
  )




  M0 = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary=M0.domain,
      offset=c( 10, 100 ),  # how much to extend in the c(inner, outer) domains
      max.edge=c( 10, 100 ),  # max size of a triange (in, out)
      min.angle=c(21),   # min angle (in, out)
      cutoff=10 # min distance allowed
  )


# plot(M0, asp=1 ) # visualise mesh
    t0$yrindex = t0$yr - min(t0$yr) + 1
    t0$pryr = t0$dyear
    t0$ti = t0$yr + t0$pryr
    t0$logZ = log(t0$z)
    t0$b0 = 1  # intercepts

  tvar = var(t0$t)

  sdCL = tvar * c(0.1, 0.9 ) # 95% CL of SD
  rangeCL = c( 20, 200 )

  spatsd = estimParamsFromCI( sdCL, pm0=c(0.1, 0.1) )
  spatrange =  estimParamsFromCI( rangeCL, pm0=c(0.1, 0.1) )
  hyper.space = list(
    range=list( param=spatrange, prior="loggamma" ),
    prec =list( param= spatsd, prior="loggamma" ) )


  # SPDE components
  # matern representation using mesh M
  S0 = inla.spde2.matern( M0, alpha=2 ) # alpha=2 is exponential correlation function


  sigma0 = 1  # SD of the latent field
  maxdist = sqrt( diff(range( locs0[,1]))^2 +  diff(range( locs0[,2]))^2  )
  range0 = maxdist / 5  # spatial range
  kappa0 = sqrt(8) / range0
  tau0 = 1/ ( sqrt(4*pi) * kappa0 *sigma0 )

  S0 = inla.spde2.matern( M0, alpha=2,
    B.tau =cbind(log(tau0), -1,1 ),     # parameter basis functions
    B.kappa = cbind( log(kappa0), 0, 1 ), # parameter basis functions
    theta.prior.mean = c(0,0),    # theta1 controls variance .. vague; theta2 controls range   .. means 0
    theta.prior.prec = c(0.1, 1)  #  precisions are vague for theta1;  for range .. theta2 prec 1 ==> 95% prior prob that range is smaller than domain size

  )

  # indices of SPDE
  # only time and mesh dependence ,, not location!
  i <- inla.spde.make.index('i', n.spde=S0$n.spde )

  # projection matrix A to translate from mesh nodes to data nodes
  A = inla.spde.make.A( mesh=M0, loc=locs0 )

  # data stack for occurence (PA)
  Z = inla.stack(
      tag="tdata",
      data=list( tC=t0$t ) ,
      A=list(A,1),
      effects=list( i=i, t0 )
  )

  theta.ar1 = list( theta1=list(param=c(1, 0.1) ), rho=list( param=c( 0.9, 1/0.1 ) ) )  # N( mean=0.9, var=0.1)

  R <- inla(
      tC ~ 0 + b0 + f(i, model=S0),
      family='gaussian',  # log transf by default .. (?)
      data=inla.stack.data(Z),
   #   control.compute=list(dic=TRUE),
      control.predictor=list(A=inla.stack.A(Z), compute=TRUE),
      verbose=TRUE
  )

    # Prediction onto a grid using projections (method A/B) is faster/computationaly more efficient .. method of choice
    pG = inla.mesh.projector( M0, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats) )
    # pG = inla.mesh.projector( M, dims=c(p$nplons, p$nplats) )

    Pmean = inla.mesh.project( pG, R$summary.random$i$mean )  # posterior mean
    Psd = inla.mesh.project( pG, R$summary.random$i$s )       # posterior SD

    image(log(Pmean))

    graphics.off()
    # field parameters on user scale
    oo = inla.spde2.result(R, 'i', S0, do.transf=TRUE)
    exp(oo$summary.log.range.nominal)
    exp(oo$summary.log.kappa)
    exp(oo$summary.log.tau)


    plot(oo$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[x]), ylab='Density')
    abline(v=exp(oo$summary.log.variance.nominal$mean) )

    plot(oo$marginals.kappa[[1]], type='l', xlab=expression(kappa), ylab='Density')
    abline(v=exp(oo$summary.log.kappa$mean) )

    plot(oo$marginals.range.nominal[[1]], type='l', xlab='range nominal', ylab='Density')
    abline(v=exp(oo$summary.log.range$mean) )
    # or
    abline(v=sqrt(8)/exp(oo$summary.log.kappa$mean) )


    # indices for random field at data locations
    idat <- inla.stack.index( Z, 'data')$data

    # correlation between the the posterior mean and the response by
    cor( z, R$summary.linear.predictor$mean[idat])



    plot( R$marginals.fixed[[1]], pch="o", xlab="Beta" )


    str(R$marginals.hyperpar)
    plot( R$marginals.hyper[[1]], type="l", xlab="" )
    plot.default( inla.tmarginal( function(x) {1/exp(x)}, R$marginals.hyperpar[[1]]), xlab="", type="l")



  # -------------------
  # 4. spatial-temporal model on temperature

  p = list()
  p$libs = RLibrary( c( "gstat", "sp", "parallel", "mgcv", "bigmemory", "INLA", "lattice" ) )
  p$libs = c( p$libs, bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.temperature", "bio.polygons" ) )

  p$tyears = c(1990:1995)  # 1945 gets sketchy -- mostly interpolated data ... earlier is even more sparse.
  p$nw = 10
  p = spatial.parameters( p=p, type="SSE" ) #  type="canada.east"  can be completed later (after assessment) when time permits if required

  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

  t0 = hydro.db( p=p, DS="bottom.gridded.all"  )

  # scotianshelf = locator()
  scotianshelf = read.table( find.bio.gis( "scotia.fundy.with.buffer.dat"  ) )
  names( scotianshelf) = c("lon", "lat")

    plot( t0$lat ~ t0$lon, pch="." )
    lines( scotianshelf )

  a = which( point.in.polygon( t0$lon, t0$lat, scotianshelf$lon, scotianshelf$lat ) != 0 )
  t0 = t0[a,]

  t0$yrindex = t0$yr - min(t0$yr) + 1
  t0$pryr = t0$dyear
  t0$ti = t0$yr + t0$pryr
  t0$logZ = log(t0$z)
  t0$b0 = 1  # intercepts

  if (anyDuplicated( paste( t0$plon, t0$plat, t0$yr, t0$dyear, sep="_"))) stop("Dups!")


  # boundary domain
  locs0  = as.matrix( t0[,c("plon", "plat")] )
  M0.domain = inla.nonconvex.hull( locs0, convex=10, resolution=200 )

  cutoff = 8
  max.edge = c(8, 80)
  offset = c(8, 8)

  debug=T
  if (debug) {
    cutoff=15
    max.edge = c(15, 120)
    offset = c(15, 120)
  }


  M0 = inla.mesh.2d (
    loc=locs0, # locations of data points
    boundary=M0.domain,
    offset=offset,  # how much to extend in the c(inner, outer) domains
    max.edge=max.edge,  # max size of a triange (in, out)
    min.angle=c(22),   # min angle (in, out)
    cutoff=cutoff # min distance allowed  /.... use 8 or less for production
  )

  plot(M0, asp=1 ) # visualise mesh

  # SPDE components
  # matern representation using mesh M
  #   spatial scale parameter kappa(u)
  #     variance rescaling parameter tau(u)
  #     (kappa^2(u)-Delta)^(alpha/2) (tau(u) x(u)) = W(u)
  #
  # hyperparamters for matern2d:
  #   theta1[log prec], prior=loggamma, param=c(1, 5e-5), initial=4
  #   theta2[log range], prior=loggamma, param=c(1, 0.01), initial=2
  #   hyper.spde2 = list( theta1=list(param=c(1,0.001), theta2=list(param=c(1,0.001))) )
  # BUT... hyper expects only one as follows! why?
  hyper.spde2 = list( theta=list(param=c(1,0.001)))

  S0 = inla.spde2.matern( M0, alpha=2 ) # alpha=2 is exponential correlation function

  # indices of SPDE
  # only time and mesh dependence ,, not location!
  ny = length(unique(t0$yr))
  i <- inla.spde.make.index('i', n.spde=S0$n.spde, n.group=ny )

  # projection matrix A to translate from mesh nodes to data nodes
  A = inla.spde.make.A( mesh=M0, loc=locs0, n.group=ny, group=t0$yrindex )

  tC=t0$t
#  varstokeep = c("b0", "pryr", "yrindex" )
  varstokeep = c("plon", "plat", "yrindex", "pryr", "b0", "z" )
  t0 = t0[,varstokeep]

  gc()

  # data stack for occurence (PA)
  Z = inla.stack(
      tag="tdata",
      data=list( tC=tC ) ,
      A=list(A,1),
      effects=list( i=i, t0 )
  )

  # hyperparmeters for ar1:
  #   theta1[prec] prior=loggamma, param=c(1, 5e-5) [mean, precision], initial=4
  #   theta2[rho]  prior=normal, param=c(0, 0.15) [mean, precision?], initial=2 [[ logit lag 1 correlation]]

Max.post.marg(theta): log(dens) = -54155.754292 fn = 774

theta =  -0.817457 1.335154 -4.208128 4.837452 0.095769 3.823501

List of hyperparameters:
		theta[0] = [Log precision for the Gaussian observations]
		theta[1] = [Theta1 for i]
		theta[2] = [Theta2 for i]
		theta[3] = [Group rho_intern for i]
		theta[4] = [Log precision for pryr]
		theta[5] = [Rho_intern for pryr]


  # hyper.year   = list( theta1=list(param=c(1, 0.01), initial=1 ), rho=list( param=c( 0.9, 1 ), initial=0.5 ) )
  hyper.year   = list( rho=list( param=c( 4, 0.1 ), initial=4 ) )
  # hyper.season = list( theta1=list(param=c(1, 0.01), initial=1 ), rho=list( param=c( 0.8, 1 ), initial=0.5 ) )
  hyper.season = list( rho=list( param=c( 4, 0.1 ), initial=4 ) )

  gc()

  # 1990:2000 -- 11 years takes ~ 1 hr
  # to complete 2000:2013
  #   mesh resolutionat 5 .. 60GB required and 10hrs
  #   mesh resoltuion at 15 .. 50 GB and 7hrs
  R <- inla(
      tC ~ 0 + b0
             + f( z, model='rw2' )
             + f( i, model=S0, hyper=hyper.spde2,
                 group=i.group, control.group=list(model='ar1', hyper=hyper.year ))
             + f( pryr, model='ar1', cyclic=TRUE, hyper=hyper.season ) ,
#      family='gaussian',  # log transf by default .. (?)
      data=inla.stack.data(Z),
#      control.compute=list(dic=TRUE, mlik=TRUE, openmp.strategy='huge'),
#      control.compute=list(dic=TRUE, mlik=TRUE),
             # control.compute = list(cpo=TRUE, pit=TRUE ),  # cpo=conditional predictive ordinate .. leave one out measures of fit to id extreme values (p(y_i|y_{-i}) .. ie. posterior value; # PIT=probability Integral Transforms Pr( y_i {new} <= y_i | y_{-i} ) .. ie on pr scale

      #   quantiles=NULL,
   #      control.results=list(return.marginals.random=FALSE, return.marginals.predictor=FALSE ),
      control.predictor=list(A=inla.stack.A(Z), compute=TRUE),
   #   control.inla=list(strategy='gaussian'),
   #   control.inla=list( h=0.002, restart=3, stupid.search=FALSE, stencil=7),
    control.inla = list( h=0.002, strategy="laplace", npoints=21, stencil=7 , strategy='gaussian'),  # more points for tails (default is 9)
    #  control.inla=list( stencil=7),
   #   num.threads=ncpu,
   #   working.directory = working.directory,
      verbose=TRUE
  )

  fn = "~/tmp/R.spatio.temporal.1990.2000.rdata"
  save(R, file=fn, compress=TRUE )
  # load(fn)

summary(R)
R$summary.hyperpar

~ 1.5 days to complete
Max.post.marg(theta): log(dens) = -126269.867546 fn = 1003
theta =
-0.922535
1.191479
-3.791064
3.987786
0.096397
3.965233
       List of hyperparameters:
                theta[0] = [Log precision for the Gaussian observations]
                theta[1] = [Theta1 for i]
                theta[2] = [Theta2 for i]
                theta[3] = [Group rho_intern for i]
                theta[4] = [Log precision for pryr]
                theta[5] = [Rho_intern for pryr]



# seasonal effect
graphics.off()

plot(R$summary.random$pryr$mean)
lines(R$summary.random$pryr$mean + R$summary.random$pryr$sd*2)
lines(R$summary.random$pryr$mean - R$summary.random$pryr$sd*2)
round(res$summary.hy[3,], 5)  # temporal autocorrelation

# annual effect



#correlation between data and predictions
idat <- inla.stack.index( Z, 'tdata')$data
cor( tC, R$summary.linear.predictor$mean[idat], use="pairwise.complete.obs" ) # 89%
plot( tC, R$summary.linear.predictor$mean[idat] )

cor( R$summary.linear.predictor$mean[idat], R$summary.fitted.values, use="pairwise.complete.obs" ) # 89%

# prediction of the random field: two ways .. with analysis or after analysis
    # (aside: mesh points are already predicted )

    newlocations = rbind(  c( 330, 4790), c(788, 5080) )
    Anew = inla.spde.make.A ( mesh=M, loc=newlocations )  # sparse matrix of the projection coefficients


    # A) simple projection, post-analysis
      Anew %*% R$summary.random$i$mean

      # or more directly:
      R$summary.fix[1,1] + Anew %*% R$summary.random$i$mean  # for the intercept
      sqrt( 1^2 + R$summary.fix[1,2]^2 + drop( Anew %*% R$summary.random$i$sd^2 ) ) # for SE



    # B) or using inla functionality
      # inla.mesh.project( inla.mesh.projector( M, loc=newlocations), R$summary.random$i$mean )

      pG = inla.mesh.projector( M0, xlim=(p$corners$plon), ylim=(p$corners$plat), dims=c(p$nplons, p$nplats) )
      x11()
      out = inla.mesh.project( pG, R$summary.random$i$mean[ i$i.group==8 ] ) # first
      levelplot( out, aspect="iso", at=seq(-5.5, 8, 0.2 )  )

      # Prediction onto a grid using projections (method A/B) is faster/computationaly more efficient .. method of choice
      pG = inla.mesh.projector( M0, xlim=p$corners$plon, ylim=p$corners$plat, dims=c(p$nplons, p$nplats) )
      # pG = inla.mesh.projector( M0, dims=c(p$nplons, p$nplats) )

      # to get climatology .. must find overall mean ... todo
      Pmean = inla.mesh.project( pG, R$summary.random$i$mean )  # posterior mean
      Psd = inla.mesh.project( pG, R$summary.random$i$s )       # posterior SD
      image(log(Pmean))



    # C) or using inla in a global analysis

    # simple model so this is not needed but to be consistent:
      i <- inla.spde.make.index('i', n.spde=S$n.spde )
      # i = 1:S$n.spde  # alternatively this
      Znew = inla.stack(
        data = list( B=NA),
        A = list( Anew ),
        effects = list( i=i, m=rep(1, nrow(newlocations) ) ),
        tag="spatialpredictions"
      )
      Z = inla.stack( Z, Znew )
      # refit the same model
      ipredictions = inla.stack.index( Z, tag="spatialpredictions")
      R$summary.linear.predctor[ ipredictions$data ,]

      # Prediction of the response variable (B) .. ie. imputation
      Zimp = inla.stack(
        data=list(B=NA),
        A=list( Anew, 1 ), # 1 is for the intercept
        effects=list( i=1:S$n.spde, m=rep(1, nrow(newlocations) ) ), # 1 is for the intercept
        tag="prediction.response"
      )

      Z = inla.stack( Z, Zimp )
      R = inla(
        B ~ 0 + m + f(i, model=S),
        data=inla.stack.data(Z),
        control.predictor=list( A=inla.stack.A(Z), compute = TRUE ) ## new addition here to compute linear predictions but this is CPU expensive
      )

      iimputations = inla.stack.index( Z, tag="prediction.response")
      R$summary.fitted.values[ iimputations$data ,]
      oo = R$marginals.fitted.values[ iimputations$data ]
      inla.hpdmarginal( 0.95, oo[[2]] )
      inla.zmarginal( oo[[2]])




      --------

  shape = 2

  formula = ~ + f(space, model="matern2d", nu=shape, hyper=hyper.space )


