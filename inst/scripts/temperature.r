
  # ----------------
  # Prep OSD, snow crab and groundfish temperature profiles
  # this one has to be done manually .. no longer mainted by anyone ..

  current.year=2016


  p = bio.temperature::temperature.parameters( current.year=current.year )

  # ------------------------------

  if ( create.baseline.database ) {

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

    p = bio.temperature::temperature.parameters( current.year=current.year )
    p = bio.temperature::temperature.parameters( DS="lbm", p=p )
    
    # 1. grid bottom data to a reasonable internal spatial resolution ; <1 min
    p = make.list( list( yrs=p$tyears), Y=p )
    # parallel.run( hydro.db, p=p, DS="bottom.gridded.redo" )
    hydro.db( p=p, DS="bottom.gridded.redo" )  # all p$tyears, for a single year use with yr argument:
    # hydro.db( p=p, DS="bottom.gridded.redo", yr=p$newyear ) 
    hydro.db( p=p, DS="bottom.gridded.all.redo" )  # all p$tyears, for a single year use with yr argument: 
    # hydro.db( p=p, DS="bottom.gridded.all.redo", yr=p$newyear ) 

    # 2. lbm interpolations assuming some seasonal pattern
    # 1950-2013, SSE took ~ 35 hrs on laptop (shared RAM, 24 CPU; 1950-2013 run April 2014 ) ... 17 GB req of shared memory
    # 1950-2015, SSE 22 hrs, 42 GB RAM, 8 CPU on hyperion (10 Jan 2015), using NLM .. not much longer for "canada.east"

    # p$clusters = c( rep("kaos",16), rep("nyx",16), rep("tartarus",16), rep("hyperion", 4), rep("io", 6) ) # with no clusters defined, use local cpu's only
    
    p = bio.temperature::temperature.parameters( current.year=current.year )
    p$lbm_local_modelengine = "twostep"
    p = bio.temperature::temperature.parameters( DS="lbm", p=p )
   
    DATA='temperature.db( p=p, DS="lbm.inputs" )' 
    p = lbm( p=p, tasks=c("initiate"), DATA=DATA ) # no global model, 5 min
    p = lbm( p=p, tasks=c( "stage1" ) ) #  ~24 hrs 
    p = lbm( p=p, tasks=c( "stage2" ) )
    p = lbm( p=p, tasks=c( "stage3" ) )
    p = lbm( p=p, tasks=c( "save" ) )


    # 3.  collect predictions from lbm and warp/break into sub-areas defined by p$spatial.domain.subareas = c( "SSE", "SSE.mpa", "snowcrab" ) (in paramters() )
    
    p$clusters = rep("localhost", detectCores() )
    p = make.list( list( yrs=p$tyears), Y=p )
    parallel.run( temperature.db, p=p, DS="predictions.redo" )
    


    # 4. extract relevant statistics:: only for default grid 
    # or parallel runs: ~ 1 to 2 GB / process
    # 4 cpu's ~ 10 min
    p$clusters = rep("localhost", detectCores() )
    p = make.list( list( yrs=p$tyears), Y=p )
    parallel.run( temperature.db, p=p, DS="bottom.statistics.annual.redo" )


    # 5. some time slices at specific weeks for fast data lookup for modelling
    temperature.db( p=p,  DS="timeslice.redo" )


    # 6. complete statistics and warp/regrid database ... ~ 2 min :: only for  default grid . TODO might as well do for each subregion/subgrid
    temperature.db( p=p, DS="complete.redo")


  # 7. maps 
 
     # p$clusters = rep("localhost", detectCores() )  # run only on local cores ... file swapping seem to reduce efficiency using th
    # p$clusters = c( rep("kaos",23), rep("nyx",24), rep("tartarus",24) )
    p = make.list( list( yrs=p$tyears), Y=p )
    for ( gr in p$spatial.domain.subareas ) {
      print (gr)
      p = spatial_parameters(  p=p, type= gr )
      parallel.run( hydro.map, p=p, type="annual"  )
      parallel.run( hydro.map, p=p, type="global")
    }
    

  # finished 



