  p=list()
  p$init.files =  loadfunctions( c("temperature", "utility", "parallel" ) ) 
  p$clusters = c( rep( "kaos", 16), rep("nyx", 16), rep("tartarus", 16) )
  p = make.list( list( yrs=yrs), Y=p )
  parallel.run( moord.db, DS="redo.all" )
 


