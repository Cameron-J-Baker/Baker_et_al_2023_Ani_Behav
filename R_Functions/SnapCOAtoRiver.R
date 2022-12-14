# For studies within riverine or other linear movement system, this function snaps the points produced along the river
# system

SnapCOAtoRiver <- function(coa, river, max_dist = 1000, UTM = CRS("+init=epsg:32754")){
  
  ########################################################################################
  # coa : Output from the COA (VTrack) or COAV2 function
  # river : A line file of the river system to snap the points to
  # max_dist : The maximum distance in m that points should be snapped to the river from
  # UTM : The UTM epsg code of your relevant location. Deafualt is for the Wenlock and Ducie river systems
  ########################################################################################
  
  UTM <- UTM #a value to convert data to the UTM projection
  
  coasp <- coa # Create a version of COA to convert to a spatial points dataframe
  coordinates(coasp) <- c( "Longitude.coa", "Latitude.coa") # convert the duplicate to a spatialpoints dataframe
  proj4string(coasp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # project it to WGS84
  coasp <- spTransform(coasp,UTM) # reproject to UTM
  coasp <- coasp %>% st_as_sf(coords = c("Longitude.coa", "Latitude.coa")) # convert to a sf object
  
  riverUTM <- spTransform(river,UTM) # ensure that the line file is projected to UTM
  riverUTM <- river %>% st_as_sf(coords = c("Lon", "Lat")) # convert to a sf object
  
  # Create the snap points to line function. Function created by TimSalabim
  st_snap_points <- function(x, y, max_dist) {
    
    if (inherits(x, "sf")) n = nrow(x)
    if (inherits(x, "sfc")) n = length(x)
    
    out = do.call(c,
                  lapply(seq(n), function(i) {
                    nrst = sf::st_nearest_points(st_geometry(x)[i], y)
                    nrst_len = st_length(nrst)
                    nrst_mn = which.min(nrst_len)
                    if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                    return(st_cast(nrst[nrst_mn], "POINT")[2])
                  })
    )
    return(out)
  }
  
  # Snap the COA points to the river lines
  COAsnap <- st_snap_points(coasp, riverUTM, max_dist) # GPS Points snapped to the wenlock lines
  
  COAsnap <- do.call(rbind.data.frame, COAsnap) # convert the list output by st_snap_points to a dataframe
  names(COAsnap) <- c("x", "y") # rename the columns to x and y
  
  # convert back from UTM to WGS84
  coordinates(COAsnap) <- c( "x", "y") # make the data frame a spatial points object
  proj4string(COAsnap) <- UTM # set the proj4string to UTM
  # convert from UTM back to WGS84
  COAsnap2 <- spTransform(COAsnap, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  COAsnap3 <- data.frame(COAsnap2)
  
  
  # Extract the new Latitude and longitude and add it to the COA data
  coa$Longitude.coa.snap <- COAsnap3$x
  coa$Latitude.coa.snap <- COAsnap3$y
  
  return(coa)
}
