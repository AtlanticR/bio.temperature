
  # ----------------
  # Prep OSD, snow crab and groundfish temperature profiles
  # this one has to be done manually .. no longer mainted by anyone ..

  current.year=2016

  p = bio.temperature::temperature.parameters( current.year=current.year )

  # ------------------------------

  if ( create.baseline.database ) {
    # 0 data assimilation
    
    if (historical.data.redo) {
      hydro.db( DS="osd.rawdata.allfiles.redo", p=p )   # redo whole data set (historical) from 1910 to 2010
      hydro.db( DS="osd.initial", p=p ) # 2008:2015
      hydro.db( DS="ODF_ARCHIVE", p=p, yr=1969:2015 ) # specify range or specific year
    }
    # Roger Petipas has been maintaining a database, the following loads this data
    hydro.db( DS="osd.current", p=p, yr=2014:p$newyear ) # specify range or specific year
    hydro.db( DS="ODF_ARCHIVE", p=p, yr=p$newyear ) # specify range or specific year

    # Merge depth profiles from all data streams: OSD, groundfish, snowcrab, USSurvey_NEFSC
    p = make.list( list( yrs=c(2008:p$newyear)), Y=p )   # specify range or specific year
    p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce ep = make.list( list( yrs=c(2008:p$newyear), Y=p ))   # specify range or specific year
    hydro.db( DS="profiles.annual.redo", p=p  )  # specify range or specific year
    # parallel.run( hydro.db, p=p, yr=p$tyears, DS="profiles.annual.redo" )

    # Extract bottom data from each profile
    p = make.list( list( yrs=2008:p$newyear), Y=p )  # specify range or specific year
    hydro.db( DS="bottom.annual.redo", yr=2008:p$newyear, p=p ) # yr argument overrides p$tyears .. e.g. for a new year of data
    # hydro.db( DS="bottom.annual.redo", p=p )
    # parallel.run( hydro.db, p=p, yr=p$tyears, DS="bottom.annual.redo" )
  }

  # ------------------------------
  # Basic data uptake now complete  .. move to interpolations
  # ------------------------------

  if (create.interpolated.results.lbm ) {
    
    # 1. lbm interpolations assuming some seasonal pattern
    # 1950-2013, SSE took ~ 35 hrs on laptop (shared RAM, 24 CPU; 1950-2013 run April 2014 ) ... 17 GB req of shared memory
    # 1950-2015, SSE 22 hrs, 42 GB RAM, 8 CPU on hyperion (10 Jan 2015), using NLM .. not much longer for "canada.east"

    # p$lbm_local_modelengine = "twostep" -- with krige would take months ... ignore for now
    # p$lbm_local_modelengine = "gam"

    p$clusters = rep("localhost", 2 ) 
    
    p$lbm_spate_method="mcmc_fast" 

    p = bio.temperature::temperature.parameters( DS="lbm", p=p )
   

    DATA='temperature.db( p=p, DS="lbm.inputs" )' 
    lbm( p=p, tasks=c("initiate"), DATA=DATA ) # no global model, 5 min
    lbm( p=p, tasks=c( "stage0" ) ) #  24 hrs for gam
    lbm( p=p, tasks=c( "stage1" ) ) #  24 hrs for gam
    lbm( p=p, tasks=c( "stage2" ) ) #   3.5 hrs for gam
    lbm( p=p, tasks=c( "stage3" ) )
    lbm( p=p, tasks=c( "save" ) )

    # to view progress in terminal:
    # watch -n 120 cat /home/jae/bio.data/bio.temperature/modelled/t/canada.east/lbm_current_status


    # 2.  collect predictions from lbm and warp/break into sub-areas defined by 
    #     p$spatial.domain.subareas = c( "SSE", "SSE.mpa", "snowcrab" ) 
    p = make.list( list( yrs=p$tyears), Y=p )
    parallel.run( temperature.db, p=p, DS="predictions.redo" ) # 10 min
    temperature.db( p=p, DS="lbm.stats.redo" ) # warp to sub grids

    # 3. extract relevant statistics 
    # or parallel runs: ~ 1 to 2 GB / process .. ~ 4+ hr
    parallel.run( temperature.db, p=p, DS="bottom.statistics.annual.redo" ) 

    # 4. all time slices in array format
    temperature.db( p=p,  DS="spatial.annual.seasonal.redo" )

    # 5. time slice at prediction time of year
    temperature.db( p=p,  DS="timeslice.redo" )

    # 6. complete statistics and warp/regrid database ... ~ 2 min :: only for  default grid . TODO might as well do for each subregion/subgrid
    temperature.db( p=p, DS="complete.redo")


  # 7. maps 
    current.year=2016
    p = bio.temperature::temperature.parameters( current.year=current.year )
    p = bio.temperature::temperature.parameters( DS="lbm", p=p )
    # p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce efficiency using th
    # p$clusters = c( rep("kaos",23), rep("nyx",24), rep("tartarus",24) )
    p = make.list( list( yrs=p$tyears), Y=p )
 
    temperature.map( p=p )

    # just redo a couple maps for ResDoc
    p$spatial.domain = "SSE"
    #p$bstats = "tmean"
    p = spatial_parameters( p=p, type=p$spatial.domain )  # default grid and resolution
    p$corners = data.frame(plon=c(150, 1022), plat=c(4600, 5320) )
    temperature.map( p=p, DS='climatology' )
    
    temperature.map( p=p, DS='annual' )
    
    temperature.map( p=p, DS='climatology' )


  # finished 



