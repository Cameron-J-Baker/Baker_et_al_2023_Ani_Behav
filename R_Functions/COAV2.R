# The below function is a modification of the VTrack function "COA" to better allow centers of activity for either ATF
# or VTrack format data. The original code was developed by Vinay Udawyer and Colin Simpledorfer. Modifications to the 
# code were completed by Cameron Baker, with this version of the function to be incorporated into the R package VTrack
# in the future
COAV2 <- function (tagdata, pointsdata, VTrack = TRUE,  id, timestep, ...){
  
  ########################################################################
  # tagdata : The acoustic telemetry data in either ATF or VTrack format that you would like to determine COA
  # pointsdata : A dataframe containing the name and location (lat/long) of acoustic hydrophone stations
  # VTrack : True false statement. True if data is VTrack format or False if it is ATF
  # id : Name of the column in tagdata that contains TRANSMITTERID, optional
  # timestep : The time step in minutes that you wish to create the COA
  ########################################################################
  
  if (VTrack == TRUE){
    data <- as.data.frame(tagdata)
    data <- plyr::join(data, pointsdata, by = "STATIONNAME")
    data$dt <- ymd_hms(data$DATETIME)
    data$TRANSMITTERID <- droplevels(as.factor(data$TRANSMITTERID))
    step_sec <- timestep * 60
    ex <- seq(from = trunc(min(data$dt, na.rm = TRUE), "day"), 
              to = trunc(max(data$dt, na.rm = TRUE), "day") + 86400, 
              by = step_sec)
    data$DATETIME <- cut(data$dt, breaks = ex)
    Latitude <- Longitude <- UNITS1 <- SENSOR1 <- TRANSMITTERID  <- NULL # creates a series of empty vectors
    cenac <- plyr::ddply(data, c("DATETIME", id), summarize, TRANSMITTERID = TRANSMITTERID[1], 
                         SENSOR1.coa = mean(SENSOR1), UNITS1 = UNITS1[1], 
                         Latitude.coa = mean(LATITUDE, na.rm = T), Longitude.coa = mean(LONGITUDE, 
                                                                                        na.rm = T), ...)
    cenac <- cenac[!is.na(cenac$Latitude.coa), ]
    if (length(levels(cenac[, id])) > 1) {
      cenac <- dlply(cenac, id)
    }
    return(cenac)
  } else {
    data <- as.data.frame(tagdata)
    data$dt <- ymd_hms(data[, grep(c("Date"), colnames(data))])
    data[, id] <- droplevels(as.factor(data[, id]))
    step_sec <- timestep * 60
    ex <- seq(from = trunc(min(data$dt, na.rm = TRUE), "day"), 
              to = trunc(max(data$dt, na.rm = TRUE), "day") + 86400, 
              by = step_sec)
    data$DateTime <- cut(data$dt, breaks = ex)
    Latitude <- Longitude <- Sensor.Unit <- Sensor.Value <- Transmitter <- Transmitter.Name <- Transmitter.Serial <- NULL
    cenac <- plyr::ddply(data, c("DateTime", id), summarize, Transmitter = Transmitter[1], 
                         Transmitter.Name = Transmitter.Name[1], Transmitter.Serial = Transmitter.Serial[1], 
                         Sensor.Value.coa = mean(Sensor.Value), Sensor.Unit = Sensor.Unit[1], 
                         Latitude.coa = mean(Latitude, na.rm = T), Longitude.coa = mean(Longitude, 
                                                                                        na.rm = T), ...)
    cenac <- cenac[!is.na(cenac$Latitude.coa), ]
    if (length(levels(cenac[, id])) > 1) {
      cenac <- dlply(cenac, id)
    } 
    return(cenac)
  }
}
