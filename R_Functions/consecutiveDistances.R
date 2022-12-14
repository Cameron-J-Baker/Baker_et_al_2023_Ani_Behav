# Function that calculates the river distance between conseuctive points
# Function to calculate river distances for each animal
consecDistTravelled <- function(sdata, tr)
{
  # sdata : A dataframe output from the COA and SnapCOAtoRiver functions
  # A transition layer of the area of interest
  
  xy <- sdata
  coordinates(xy) <- c("Longitude.coa.snap", "Latitude.coa.snap") # convert the data to a spatial points format
  proj4string(xy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #give it a proj4string to WGS84
  xy <- spTransform(xy,UTM) # project to UTM
  
  SP <- SpatialPoints(xy,proj4string=UTM) # retrive only the spatial points
  
  # Initialise variables
  iCount <- 1
  rDISTANCE <- rep(0,nrow(xy)-1)
  sDISTANCE <- rep(0,nrow(xy)-1)
  
  # For each line in track file
  while (iCount <= nrow(xy))
    #while (iCount <= 20)
    
  {
    if(iCount == 1){
      
      sDISTANCE[iCount] <- 0     
      rDISTANCE[iCount] <- 0 
      
      iCount <- iCount + 1
    } else {
      
      SP1 <- SP[iCount-1,]
      SP2 <- SP[iCount,]
      #plot(SP1,add=T,col=3)      
      #plot(SP2,add=T,col=4)
      
      crowdist <- spDists(SP1,SP2)
      sPath2 <- shortestPath(x = tr, origin = SP1, goal = SP2, output="SpatialLines")
      riverdist <- SpatialLinesLengths(sPath2)
      
      sDISTANCE[iCount] <- crowdist      
      rDISTANCE[iCount] <- riverdist 
      
      iCount <- iCount + 1
    }
  }
  sdata$rDISTANCE <- rDISTANCE/1000 #divide by 1000 to convert from m to km
  sdata$sDISTANCE <- sDISTANCE/1000 #divide by 1000 to convert from m to km
  return(sdata)
}