# Author: Cameron J Baker
# 2022-12-12
# Email: cameron.baker@uqconnect.edu.au/cameron.baker@cdu.edu.au
# Code for the analyses contained within the manuscript "Long-term tracking reveals a dynamic crocodylian social system" published in
#  the journal Animal Behaviour
# Aim: This script examines whether the social structure of crocodiles is different from random chance following the methods described in the 
# "Is the social structure of estuarine crocodiles different from random?" methods section

library(data.table)
library(lubridate)
library(foreach)
library(doParallel)
library(sf)
library(raster)
library(mapview)
library(data.table)
library(rgdal)
library(rgeos)
library(wicket)
library(ggmap)
library(lubridate)
library(gdistance)
library(tibble)
library(VTrack)
library(spacetime)
library(plotKML)
library(data.table)
library(sp)

######################################################################################################################
#### Import and prepare the data for analysis ####

# Import the raw data
croc_data <- fread("Baker_et_al_data.csv")
croc_data$Month <- quarter(croc_data$DATETIME, with_year = T) # Create a new column that splits each year a crocodile is tracked into quarters to allow the 
                                                      # the calculation of the required metrics for the spatially-explict null model.

####===============================================================================================================####
# Import a file with the distance from the mouth of Port Musgrave to each receiver
receiver_loc <- fread("Data/Baker_et_al_acoustic_hydrophone_location.csv")

# Import the distance matrix between each receiver
distmatrix <- fread("Data/Baker_et_al_hydrophone_distance_matrix.csv")
# Convert the first column to row names and then remove it
distmatrix <- tibble::column_to_rownames(distmatrix, var = "V1") # convert stationame back into row.names

######################################################################################################################
# Calculate the monthly movement metrics (mean speed, sd speed, probability of being observed at receivers) for each
# individual. To prevent the gross overestimation of the speed crocodiles can travel, all speed values will be calculated
# using COAs
#=======================================================================================================================#
# import the functions that will be required
#=======================================================================================================================#
source("R/Functions/COAV2.R")
source("R/Functions/SnapCOAtoRiver.R")
source("R/Functions/consecutiveDistances.R")
#====================================================================================================================#
# Import the required spatial objects and then run the COA functions below
#====================================================================================================================#
# Define UTM to project data into for calculation of home ranges
UTM <- CRS("+init=epsg:32754") #a value to convert data to the UTM projection

# These shapefiles are available from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.mw6m905vw 

# firstly clip the Cape york lines file to the Wenlock river
# Shapefile name: Baker_et_al_line_file
Wen_Riverclip <- readOGR("GIS/Wenlock-Ducie","Wenlock_Ducie_extent_Dissolv") # Load the Wenlock spatiallines  data
Wen_Riverclip <- spTransform(Wen_Riverclip,UTM) # project to UTM
plot(Wen_Riverclip)

# Import the raster of the Wenlock river to calculate distances along
# Baker_et_al_study_site_raster.asc
Wen_Raster <- raster("GIS/Rasters/WenlockRaster2020.asc")
plot(Wen_Raster)
# Convert the raster to a transition object - required for the river distance algorithm
tr <- transition(Wen_Raster, transitionFunction=mean, directions=8) # Create a Transition object from the raster

#====================================================================================================================#
# Next create the COA file for all individuals and then export it so it can be used again later on
crocCOAs <- COAV2(tagdata = croc_data, pointsdata = receiver_loc, Vtrack = TRUE, id = "TRANSMITTERID", timestep = 360) 
crocCOAs <- rbindlist(crocCOAs) # combine into a single dataframe
crocCOAs <- as.data.table(SnapCOAtoRiver(crocCOAs, Wen_Riverclip, max_dist = 5000)) # snap the COA's to the river ployline

write.csv(crocCOAs, "Data/crocodileCOA_2020.csv", row.names = F) # Write a copy to avoid having to duplicate the above steps

# crocCOAs <- fread("Data/crocCOA_Wenlock_only_2020.csv") # Import the COAs calculated from above
#=====================================================================================================================#
# This next function calcualtes the mean and SD speed of individuals per month, along with the probability of them 
# being observed at each receiver

individuals <- unique(croc_data$TRANSMITTERID)
moveMets <- function(i){
  
  # First subset out the individual
  indiv <- croc_data[TRANSMITTERID == individuals[i]] # subset out the individuals data
  indiv_coa <- crocCOAs[TRANSMITTERID == individuals[i]]
  # next determine what months has the individual been recorded within? 
  Months <- unique(indiv$Month)
  
  #indiv_coa$Month <- substr(indiv_coa$DATETIME, 1,7) # add a column of the months to determine the monthly max and mean speed
  indiv_coa$Month <- quarter(indiv_coa$DATETIME, with_year = T)
  
  # Determine the distance travelled  between consecutive detections
  indiv_coa <- consecDistTravelled(sdata = indiv_coa, tr)
  indiv_coa$Speed <- 0 # create an empty column to add the distance travelled between consecutive detections
  
  # extract the speed of each movement
  if(nrow(indiv_coa) > 1){ # if there is more than one COA
    for(e in 2:nrow(indiv_coa)){
      
      indiv_coa$Speed[e] <- indiv_coa$rDISTANCE[e] / as.numeric(difftime(indiv_coa$DATETIME[e], indiv_coa$DATETIME[e-1], units = "hours"))
    }
    
  }
  
  monthly_metrics <- list()
  
  # For each of the months, run the next part
  for(t in 1: length(Months)){
    
    # Subset out the month of interest
    indiv_month <- indiv[Month == Months[t]]
    indiv_coa_month <- indiv_coa[Month == Months[t]]
    
    # Due to how the snapCOA function works, for some months where individuals are only observed at almost the same spot
    # the function may snap all points to a single location, even though there was some movement between two receivers
    # For these cases, recalculate the speed column using the straight line distances between individuals
    if(max(indiv_coa_month$Speed) == 0 & nrow(indiv_coa_month) > 1){
      for(e in 2:nrow(indiv_coa_month)){
        
        indiv_coa_month$Speed[e] <- indiv_coa_month$sDISTANCE[e] / as.numeric(difftime(indiv_coa_month$DATETIME[e], indiv_coa_month$DATETIME[e-1], units = "hours"))
      }
    }
    
    # Next, determine the probabilty of the individual being observed at each of the receivers in the array
    indiv_month$weight <- 1 # add an extra column so that the number of detections per receiver can be counted
    indivProbs <- indiv_month[,.(Prob = sum(weight)/nrow(indiv_month)), by = STATIONNAME]
    
    # Add in a column for the month, along with the mean and maximum speeds travelled by an individual
    indivProbs$Month <- Months[t]
    indivProbs$Mean_speed = mean(indiv_coa_month$Speed)
    indivProbs$SD_speed = ifelse(is.na(sd(indiv_coa_month$Speed)) == T, 0, sd(indiv_coa_month$Speed))
    indivProbs$Max_speed = max(indiv_coa_month$Speed)
    
    monthly_metrics[[t]] <- indivProbs
  }
  
  #convert the list into a data frame
  output <- rbindlist(monthly_metrics)
  
  # Add in the ID of the indivdual and then reorder the columns and export the results
  
  output$TRANSMITTERID <- individuals[i]
  output <- output[, c(7,3,1,2,4,5,6)] # reorder the dataframe
  
  return(output)
}

#test <- moveMets(4)

registerDoParallel(5) # set the number of cores to use

lmig <- foreach(i = 1:(length(individuals)), 
                .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                              "adehabitatLT", "adehabitatHR", "rgeos", "rgdal","VTrack",
                              "geosphere","raster","gdistance","tidyverse", "sf")) %dopar%
  moveMets(i)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

move_met <- rbindlist(lmig) # Combine all of the movement metrics into a single dataframe
# Save the 
#write.csv(move_met, "Data/Movement_metrics_quarter_Wenlock_only.csv", row.names = F)

######################################################################################################################
#### Create and run the spatially explict null model for the required number of iterations ####
rand_walk <- function(i){
  
  ###################################################################
  # i: A value to indicate and individuals position along the list of individuals
  ###################################################################
  
  indiv <- crocs[TRANSMITTERID == individuals[i]] # subset out the individuals data
  indiv_metrics <- moveMetrics[TRANSMITTERID == individuals[i]]
  #indiv_metrics <- test
  # next determine what months has the individual been recorded within? 
  Months <- unique(indiv$Month)
  
  rand_data <- list()
  
  # For each of the months, run the next part
  for(t in 1: length(Months)){
    
    # Subset out the month of interest
    indiv_month <- indiv[Month == Months[t]]
    indiv_metrics_monthly <- indiv_metrics[Month == Months[t]]
    
    # If an individual has only been seen at one station for that month, assign it's detections to that period
    if(nrow(indiv_metrics_monthly) == 1){
      rand_data[[t]] <- data.table(DATETIME = indiv_month$DATETIME,
                                   STATIONNAME = indiv_metrics_monthly$STATIONNAME,
                                   TRANSMITTERID = individuals[i]) 
      
    } else {
      
      # for each observed point after the first one, determine the time difference and maximum distance that could
      # be travelled in that time period
      indiv_month$Max_distance <- 0
      for(q in 2 : nrow(indiv_month)){
        indiv_month$Max_distance[[q]] <- indiv_metrics_monthly$Max_speed[1] * 
          as.numeric(difftime(indiv_month$DATETIME[q], indiv_month$DATETIME[q-1], units = "hours")) 
        
      }
      
      # If it is the first month, set the starting point to a random value within the dataset
      if(t == 1){
        
        # create a vector to store the STATIONNAME of each assigned point to
        sequence <- NA
        
        # For the very firt point of the first month, assign position 1 randomly based on the overal probability of
        # an individual being observed at a receiver
        #sequence[[1]] <- sample(indiv_metrics_monthly$STATIONNAME, size = 1, replace = T, # n - 1 as the first position will be randomly selected
        #                       prob = indiv_metrics_monthly$Prob)
        # The start of the sequence is fixed as the first receiver was identified at
        sequence[[1]] <- indiv_month$STATIONNAME[1]
        
        # Next, for the remaining points of this month, randomly assign the movements of the individual, based on
        # how far they could possibly travel in the time period provided 
        for(y in 2:nrow(indiv_month)){
          # Subset out the distance values from the distance matrix corresponding to the individuals previous position
          # and then subset again by the maximum distance that the individual could have travelled
          indiv_dist<- as.data.frame(rownames_to_column(distmatrix[sequence[y-1]])) # from the matrix, subset out the column of the previous location
          indiv_dist <- subset(indiv_dist,indiv_dist[2] <= indiv_month$Max_distance[y])
          # subset the probabilities to only the potential receivers that the individual could move to in the time between detections
          next_point <- subset(indiv_metrics_monthly, (STATIONNAME %in% indiv_dist$rowname)) 
          # assign the next position in the sequence based on the probabilities of an individual being observed there
          sequence[[y]] <- sample(next_point$STATIONNAME, size = 1, replace = T, # n - 1 as the first position will be randomly selected
                                  prob = next_point$Prob)
          
        }
        # Add the random data sequence frame into the overal list for exporting
        rand_data[[t]] <- data.table(DATETIME = indiv_month$DATETIME,
                                     STATIONNAME = sequence,
                                     TRANSMITTERID = individuals[i]) 
        
        #testing <- data.table(DATETIME = indiv_month$DATETIME,
        #                     Simulated = sequence,
        #                     Observed = indiv_month$STATIONNAME,
        #                   TRANSMITTERID = individuals[i],
        #                     Max_dist = indiv_month$Max_distance) 
        
      } else{
        
        # create a vector to store the STATIONNAME of each assigned point to
        sequence <- NA
        
        # First step, determine the position of the first location based on the previous location of the last month
        # Extract the random sequence generated from the previous run
        previous <- rand_data[[t-1]]
        previous_metrics <- indiv_metrics[Month == Months[t-1]]
        # Determine the maximum distance that the individual could have travelled between the last detection of the last
        # month and first detection of this month
        timediff <- as.numeric(difftime(indiv_month$DATETIME[1], previous$DATETIME[nrow(previous)],units = "hours"))
        
        # Next, determine the maximum distance that the individual could have travelled using the maximum speed the 
        # crocodile could have travelled 
        maxdist <- ifelse(indiv_metrics_monthly$Max_speed[1] > previous_metrics$Max_speed[1],
                          timediff * indiv_metrics_monthly$Max_speed[1], timediff * previous_metrics$Max_speed[1])
        # Subset out the distance values from the distance matrix corresponding to the individuals previous position
        # and then subset again by the maximum distance that the individual could have travelled
        indiv_dist<- as.data.frame(rownames_to_column(distmatrix[previous$STATIONNAME[nrow(previous)]])) # from the matrix, subset out the column of the previous location
        indiv_dist <- subset(indiv_dist,indiv_dist[2] <= maxdist)
        # subset the probabilities to only the potential receivers that the individual could move to in the time between detections
        next_point <- subset(indiv_metrics_monthly, (STATIONNAME %in% indiv_dist$rowname)) 
        # assign the next position in the sequence based on the probabilities of an individual being observed there
        
        # if there is no point that connects the two years, set the first detection of that year as the start
        ifelse(nrow(next_point) == 0, sequence[[1]] <- indiv_month$STATIONNAME[1],
               sequence[[1]] <- sample(next_point$STATIONNAME, size = 1, replace = T, # n - 1 as the first position will be randomly selected
                                       prob = next_point$Prob))
        
        # Next, for the remaining points of this month, randomly assign the movements of the individual, based on
        # how far they could possibly travel in the time period provided 
        for(y in 2:nrow(indiv_month)){
          # Subset out the distance values from the distance matrix corresponding to the individuals previous position
          # and then subset again by the maximum distance that the individual could have travelled
          indiv_dist<- as.data.frame(rownames_to_column(distmatrix[sequence[y-1]])) # from the matrix, subset out the column of the previous location
          indiv_dist <- subset(indiv_dist,indiv_dist[2] <= indiv_month$Max_distance[y])
          # subset the probabilities to only the potential receivers that the individual could move to in the time between detections
          next_point <- subset(indiv_metrics_monthly, (STATIONNAME %in% indiv_dist$rowname)) 
          # assign the next position in the sequence based on the probabilities of an individual being observed there
          sequence[[y]] <- sample(next_point$STATIONNAME, size = 1, replace = T, # n - 1 as the first position will be randomly selected
                                  prob = next_point$Prob)
          
        }
        # Add the random data sequence frame into the overal list for exporting
        rand_data[[t]] <- data.table(DATETIME = indiv_month$DATETIME,
                                     STATIONNAME = sequence,
                                     TRANSMITTERID = individuals[i]) 
        
      }
    }
  }
  
  #indiv_sims[[i]] <- rbindlist(rand_data)
  testing <- rbindlist(rand_data)
  
  return(testing)
}

individuals <- unique(croc_data$TRANSMITTERID) # Create a list of all individuals present

# For loop to run the random model for the required number of times to generate the null models
for(a in 1:10){
  
  registerDoParallel(10) # set the number of cores to use
  
  lmig <- foreach(i = 1:(length(individuals)), 
                  .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                                "adehabitatLT", "adehabitatHR", "rgeos", "rgdal","VTrack",
                                "geosphere","raster","gdistance","tidyverse", "sf")) %dopar%
    rand_walk(i)
  
  stopImplicitCluster() # stop the cluster to return memory and resources to other systems
  
  sims[[a]] <- rbindlist(lmig)
  
  write.csv(sims[[a]], paste0("Data/Simulations/Datasets/Random_walk/Simulation_", a, ".csv"), row.names = F)
  
}

# Once this is complete, the associations and home ranges of each of the random datasets needs to be calculated following the previous R scripts.

######################################################################################################################################################
######################################################################################################################################################
#### Determine if the social structure of the population is different from chance ####

########################################################################################################################
# Home ranges were constructed following the methods described in Baker et. al (2022) Crocodile social environments dictated by male philopatry. Behavioral Ecology, 33(1), 156-166. doi: https://doi.org/10.1093/beheco/arab120
# The required R code and shape files can be found at:  https://github.com/Cameron-J-Baker/Supporting-code-Baker_et-al 
# Import the observed home ranges overlap matrices produced following the code above
# If a dyad was observed overlapping in space, but where not observed associating they were considered a potential dyad and are assigned an SRI of 0
filesPD <- list.files(path = "Data/Overlap_matrices/Monthly_Wenlock/", full.names = T) # create a vector list with the names of all files

potential_dyads <- lapply(filesPD, fread)
yearlyHR <- list() # a list which contains the HR dataframes that contain all of the individuals present (values not subset to be less than 1% overlap)
# Convert each of the files imported into matrices and move the TRANSMITTERIDs into the rownames
for(i in 1:length(potential_dyads)){
  potential_dyads[[i]] <- as.data.frame(potential_dyads[[i]]) # convert to data.frame
  colnames(potential_dyads[[i]]) = potential_dyads[[i]][1, ] # convert the first row back to column names
  potential_dyads[[i]] = potential_dyads[[i]][-1, ] # remove the first row with the names  
  rownames(potential_dyads[[i]]) <- c() # remove the current row names
  potential_dyads[[i]] <- tibble::column_to_rownames(potential_dyads[[i]], var = "NA") # convert iD back into row.names
  
  potential_dyads[[i]] <- as.matrix(potential_dyads[[i]])
  
  potential_dyads[[i]] <- na.omit(as.data.table(as.table(potential_dyads[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
  names(potential_dyads[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns
  
  
  # Create a dyads column so that the they are in order from lowest ID to highest to match what's in the SRI code
  potential_dyads[[i]]$dyads <- ifelse(as.numeric(potential_dyads[[i]]$ID1) < as.numeric(potential_dyads[[i]]$ID2), 
                                       paste0(potential_dyads[[i]]$ID1, "-", potential_dyads[[i]]$ID2),
                                       paste0(potential_dyads[[i]]$ID2, "-", potential_dyads[[i]]$ID1))
  
  # Remove all rows that have the HRoverlap of the same individual
  potential_dyads[[i]] <- potential_dyads[[i]][ID1 != ID2]
  
  # Remove all of the dyads that do not have at least a 1% HR overlap for that year
  #potential_dyads[[i]] <- potential_dyads[[i]][HRoverlap >= 0.01]
  
  # The code to make the df from the matrix didn't do exactly what I wanted, removing duplicated dyads
  potential_dyads[[i]] <- potential_dyads[[i]][!duplicated(potential_dyads[[i]]$dyads),]
  
  # Remove all dyads that had 0% HR overlap. Individuals who do not spatially overlap can not be potentially
  # associating with each other
  yearlyHR[[i]] <- potential_dyads[[i]]
  potential_dyads[[i]] <- potential_dyads[[i]][HRoverlap > 0.01]
  potential_dyads[[i]]$Year <- substr(filesPD[i], 39,45)
  
}

# potential_dyads contains only the dyads that overlapped at least 1 % of their home ranges
# yearlyHR contains all individuals that had a home range generated
potential_dyads <- rbindlist(potential_dyads) # combine into a single data frame for later
potential_dyads$Month <- substr(potential_dyads$Year,6,7) # next create a column for the months present
###########################################################################################################################
# Import the raw co-occurrences for the whole study period
raw_dyads <- fread("crocodile_edges_2020_4.csv")

#----------------------------------------------------------------------------#
# For each of the observed months, go through the raw data and remove any dyads where an indivdual did not have the required data to determine a HR overlap
# for that specific month
raw_dyads2 <- list()
for(i in 1:length(years1)){
  year_sri <- raw_dyads[Year == years1[i]] # subset out all of the co-detections for that year
  year_HR <- potential_dyads[Year == years1[i]] # extract the home range overlaps for the year of interest
  individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
  year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
  
  raw_dyads2[[i]] <- year_sri # add the subsetted dataset into a list to create the subset dataframe
}

raw_dyads2 <- rbindlist(raw_dyads2)

###########################################################################################################################
# Import the prepAssocIndex function to calculate the SRI for each dyad
source("R_Functions/prepAssocIndex.R")

############################################################
##### What is the overall CV and prop non-zero of the population and is that different from random?

# First determine the overall SRI's for each dyad
output <- as.data.table(prepAssocIndex(raw_dyads2)) # determine the SRI of each observed dyad
# Order the data so that each dataframe is consistent
output <- output[order(Indiv1, Indiv2),]
# Assign the dyad IDs the exact same way they were assigned for the potetnial dyads to ensure they match between the two datasets
output$dyads <- ifelse(as.numeric(output$Indiv1) < as.numeric(output$Indiv2), 
                       paste0(output$Indiv1, "-", output$Indiv2),
                       paste0(output$Indiv2, "-", output$Indiv1))

# Next, for any dyad that had an observed HR overlap but did not associate add them in with an SRI of 0

pot_dyads <- potential_dyads[!duplicated(potential_dyads$dyads), ]

#obs_dyad <- unique(observed$dyads)
# create a new dataframe that will add in all the missing dyads
missing <- list()
e <- 1
for(t in 1:nrow(pot_dyads)){
  
  if(nrow(output[dyads == pot_dyads[t,4]]) < 1){
    missing[[e]] <- data.table(Indiv1 = pot_dyads$ID1[t],
                               Indiv2 = pot_dyads$ID2[t],
                               Weight = 0,
                               Ya = 0,
                               Yb = 0,
                               Yab = 0,
                               SRI = 0,
                               dyads = pot_dyads[[t,4]])
    e <- e + 1
  }
}

missing <- rbindlist(missing)

# Add the missing dyads into the observed one
output2 <-rbind(output, missing)
# place the data in the same order as the observed
output2 <-  output2[order(Indiv1, Indiv2),]
# Add a new column where if there is an SRI greater than 0 it gets a 1. This is so it can be summed to determine the proportion of
# non-zero values within the dataset to test for long-term aversions
output2$NotZero <- ifelse(output2$SRI == 0, 0, 1)


observedCVs <- data.table(CV = coefficient.variation(sd = sd(output2$SRI), avg = mean(output2$SRI)),PropNZ = sum(output2$NotZero)/length(output2$NotZero) )

#==================================================================================================================================#
# Next repeat the same for all of the simulations
# This is a function that does all the same steps as above. If monthly = F it does the exact same and calculates values for the entire period.
# If monthly = T it then does it for each calendar month across the entire study period (i.e., all January's at once)
# Using the years1 column, import the monthly home range overlaps of each simulation, convert them to a dataframe as above and then stitch
# it together and into a list. Then import the simulated data co-occurrences, determine the SRI's and then assign the 0 values.
# Output of the function will be a dataframe with the SRI values of every possible association within the network per month

simulatedSRIs <- function(i, years1, monthly = F){ # if monthly = F it works out the overall values for the entire study period, else it does it for each calendar month across years
  
  # First step, import the home range overlaps for each month and convert to a dataframe.
  HRdataframe <- list() # create a list to store the values within
  for(t in 1:length(years1)){
    # If the file exists, import it into the list and then determine the dyads in the same manner as the observed data above
    if(file.exists(paste0("Data/Simulations/HR_overlap/Matrix/HR_overlaps_Simulation-",i,"-",years1[t],"_matrix.csv")) == T){
      sim_overlap <- fread(paste0("Data/Simulations/HR_overlap/Matrix/HR_overlaps_Simulation-",i,"-",years1[t],"_matrix.csv"))
      # Convert the matrix to a dataframe with each of the unique dyads in the same order as the observed data
      sim_overlap <- as.data.frame(sim_overlap) # convert to data.frame
      colnames(sim_overlap) = sim_overlap[1, ] # convert the first row back to column names
      sim_overlap = sim_overlap[-1, ] # remove the first row with the names  
      rownames(sim_overlap) <- c() # remove the current row names
      sim_overlap <- tibble::column_to_rownames(sim_overlap, var = "NA") # convert iD back into row.names
      sim_overlap <- as.matrix(sim_overlap)
      sim_overlap <- na.omit(as.data.table(as.table(sim_overlap, na="", row.names=T, col.names=T))) # converts from matrix to table
      names(sim_overlap) <- c("ID1", "ID2","HRoverlap") # Rename the columns
      # Create a dyads column so that the they are in order from lowest ID to highest to match what's in the SRI code
      sim_overlap$dyads <- ifelse(as.numeric(sim_overlap$ID1) < as.numeric(sim_overlap$ID2), 
                                  paste0(sim_overlap$ID1, "-", sim_overlap$ID2),
                                  paste0(sim_overlap$ID2, "-", sim_overlap$ID1))
      # Remove all rows that have the HRoverlap of the same individual
      sim_overlap <- sim_overlap[ID1 != ID2]
      # The code to make the df from the matrix didn't do exactly what I wanted, removing duplicated dyads
      sim_overlap <- sim_overlap[!duplicated(sim_overlap$dyads),]
      sim_overlap <- sim_overlap[HRoverlap > 0.01]
      sim_overlap$Year <- years1[t]
      sim_overlap$Month <- substr(years1[t],6,7) # add a column to identify the month
      HRdataframe[[t]] <- sim_overlap
      
    }
  }
  
  HRdataframes <- rbindlist(HRdataframe) # combine the list into a dataframe
  year <- unique(HRdataframes$Month) # extract out a new vector which just contains the observed months which had enough data
  
  # Next import the simulated co_occurrences and then determine the SRI while subsetting out the individuals who did not have enough
  # data to generate home range estimates
  sim_co_occur <- fread(paste0("Data/Simulations/Co_detections/Simulation-",i,".csv"))
  # Remove all of the receivers that have not been used within the study for generating home ranges
  receivers_out <- c("Ducie half way upper", "Jardines landing", "Pennefather AR", 
                     "Ducie 1", "Ducie 2", "Ducie 3", "DUC1", "DUC2", "DUC2 2", "Ducie mouth upstream",
                     "Ducie mouth downriver", "Namaleeda upstream", "Namaleeda downstream",
                     "PM Cullen AR", "PM2", "PM Ducie Namaleeta AR", "Ling creek 3", "Skardon AR","Catfish 1", "Alice Creek",
                     "Test receiver","Dalhunty twin creeks","Rocky hole","Skardon river","Janie Creek AR", "PM Mouth North AR",
                     "Ducie to PM AR", "PM Wenlock AR", "No name Creek" )
  
  sim_co_occur <- subset(sim_co_occur, !(Location %in% receivers_out)) # subset out the receivers that are not wanted
  
  sim_co_occur$Year <- substr(sim_co_occur$DatetimeID1,1,7)
  # For each of the observed months, go through the raw data and remove any dyads where an indivdual did not have the required data to determine a HR overlap
  # for that specific month
  sim_co_occur2 <- list()
  for(k in 1:length(years1)){
    year_sri <- sim_co_occur[Year == years1[k]] # subset out all of the co-detections for that year
    year_HR <- HRdataframes[Year == years1[k]] # extract the home range overlaps for the year of interest
    individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
    year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
    
    sim_co_occur2[[k]] <- year_sri # add the subsetted dataset into a list to create the subset dataframe
  }
  
  sim_co_occur <- rbindlist(sim_co_occur2)
  
  
  # Create a column within the dataset, which places the data into the correct years for the dynamic networks
  sim_co_occur$Month <- substr(sim_co_occur$DatetimeID1,6,7)
  
  if(monthly == F){
    # Remove all individuals who could not have a HR generated throughout the study
    individuals <- unique(c(HRdataframes$ID1, HRdataframes$ID2)) # create a list of the individuals with a HR estimate
    raw_dyads2 <- raw_dyads[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
    
    output <- as.data.table(prepAssocIndex(raw_dyads2)) # determine the SRI of each observed dyad
    # Order the data so that each dataframe is consistent
    output <- output[order(Indiv1, Indiv2),]
    # Assign the dyad IDs the exact same way they were assigned for the potetnial dyads to ensure they match between the two datasets
    output$dyads <- ifelse(as.numeric(output$Indiv1) < as.numeric(output$Indiv2), 
                           paste0(output$Indiv1, "-", output$Indiv2),
                           paste0(output$Indiv2, "-", output$Indiv1))
    
    # Next, for any dyad that had an observed HR overlap but did not associate add them in with an SRI of 0
    pot_dyads <- HRdataframes[!duplicated(HRdataframes$dyads), ]
    
    #obs_dyad <- unique(observed$dyads)
    # create a new dataframe that will add in all the missing dyads
    missing <- list()
    e <- 1
    for(t in 1:nrow(pot_dyads)){
      
      if(nrow(output[dyads == pot_dyads[t,4]]) < 1){
        missing[[e]] <- data.table(Indiv1 = pot_dyads$ID1[t],
                                   Indiv2 = pot_dyads$ID2[t],
                                   Weight = 0,
                                   Ya = 0,
                                   Yb = 0,
                                   Yab = 0,
                                   SRI = 0,
                                   dyads = pot_dyads[[t,4]])
        e <- e + 1
      }
    }
    
    missing <- rbindlist(missing)
    
    # Add the missing dyads into the observed one
    croc_sim_associations <-rbind(output, missing)
    # place the data in the same order as the observed
    croc_sim_associations <-  croc_sim_associations[order(Indiv1, Indiv2),]
    croc_sim_associations$Simulation <- paste0("Simulation-",i)
    
  } else {
    # Now determine the monthly SRI of the simulated data
    monthly_SRI <- list()
    for(q in 1:length(year)){ # only run this function for the simulated months which had enough data to generate home ranges
      year_sri <- sim_co_occur[Month == year[q]] # subset out all of the co-detections for that year
      year_HR <- HRdataframes[Month == year[q]] # extract the home range overlaps for the year of interest
      individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
      year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
      
      output <- as.data.table(prepAssocIndex(year_sri))
      # Order the data so that each dataframe is consistent
      output <- output[order(Indiv1, Indiv2),]
      # Assign the dyad IDs the exact same way they were assigned for the potetnial dyads to ensure they match between the two datasets
      output$dyads <- ifelse(as.numeric(output$Indiv1) < as.numeric(output$Indiv2), 
                             paste0(output$Indiv1, "-", output$Indiv2),
                             paste0(output$Indiv2, "-", output$Indiv1))
      output$Month <- year[q]
      # Save it into a list object
      monthly_SRI[[q]] <- output
      
    }
    monthly_SRI <- rbindlist(monthly_SRI) # combine the list into a dataframe
    
    # For each month, in the observed association dataframe add in all the potential dyads that were not observed with an
    # SRI value of 0
    monthly_SRIs <- list()
    for(u in 1:length(year)){
      # Extract the observed and potential dyads for this year
      observed <- monthly_SRI[Month == year[u]]
      potential <- HRdataframes[Month == year[u]]
      
      # Only run the below code if there are potential dyads observed (several months have 0 observed months but several potetnial ones)
      # to avoid errors introduced from two months (both March!) that had no overlaps greater than 0.1
      if(nrow(potential) > 0){
        # Create a list of the potential dyads
        #pot_dyads <- unique(potential$dyads)
        pot_dyads <- potential[!duplicated(potential$dyads), ]
        #obs_dyad <- unique(observed$dyads)
        # create a new dataframe that will add in all the missing dyads
        missing <- list()
        e <- 1
        for(t in 1:nrow(pot_dyads)){
          if(nrow(observed[dyads == pot_dyads[t,4]]) < 1){
            missing[[e]] <- data.table(Indiv1 = pot_dyads$ID1[t],
                                       Indiv2 = pot_dyads$ID2[t],
                                       Weight = 0,
                                       Ya = 0,
                                       Yb = 0,
                                       Yab = 0,
                                       SRI = 0,
                                       dyads = pot_dyads[[t,4]],
                                       Month = year[[u]])
            e <- e + 1
          }
        }
        missing <- rbindlist(missing)
        # Add the missing dyads into the observed one
        observed <-rbind(observed, missing)
        # place the data in the same order as the observed
        observed <-  observed[order(Indiv1, Indiv2),]
        observed$Simulation <- paste0("Simulation-",i)
        # Add back into the primary list
        monthly_SRIs[[u]] <- observed
      }
    }
    croc_sim_associations <- rbindlist(monthly_SRIs) # combine all of the simulations into a single dataframe
    
  }
  
  
  # Add a new column where if there is an SRI greater than 0 it gets a 1. This is so it can be summed to determine the proportion of
  # non-zero values within the dataset to test for long-term aversions
  croc_sim_associations$NotZero <- ifelse(croc_sim_associations$SRI == 0, 0, 1)
  
  return(croc_sim_associations)
}


# create the matrices
registerDoParallel(6) # set the number of cores to use
# Run the above function for each of the simulated datasets
sim_SRI <- foreach(p = 1:50, 
                   .packages = c("data.table","lubridate","plyr", "tibble", "dplyr")) %dopar%
  simulatedSRIs(p, years1 = years1, monthly = F)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems


#=============================================================================================#
# The below function than compares the observed to the simulated dataset and determines the p-value for each metric
sim_P_value <- function(observed, simulated){
  
  # Exact Feb, Mar, Apr from the analyses
  #observed$Month <- substr(observed$Year, 6,7)
  #observed <- observed[Month != "02" | Month != "03" | Month != "04"]
  #observed <- observed[Month == years2[i]]
  
  # Firstly, extract the overall CV and proportion of non_zeros from the observed data
  observedCV <- ifelse(is.na(coefficient.variation(sd = sd(observed$SRI), avg = mean(observed$SRI))) == T, 0,
                       coefficient.variation(sd = sd(observed$SRI), avg = mean(observed$SRI)))
  observed_prop_NZ <- sum(observed$NotZero)/nrow(observed)
  ObservedSD <- sd(observed$SRI)
  sim_values <- list()
  # For each of the simulations, calculate their CV and non-zero proportions
  for(p in 1:length(simulated)){
    
    sim <- simulated[[p]]
    #sim <- sim[Month == years2[i]]
    # Exact Feb, Mar, Apr from the analyses
    #sim$Month <- substr(sim$Year, 6,7)
    #sim <- sim[Month != "02" | Month != "03" | Month != "04"]
    sim_values[[p]] <- data.table(Simulation = paste0("Simulation-",p),
                                  CV = ifelse(is.na(coefficient.variation(sd = sd(sim$SRI), avg = mean(sim$SRI))) == T, 0,
                                              coefficient.variation(sd = sd(sim$SRI), avg = mean(sim$SRI))),
                                  prop_NZ = sum(sim$NotZero)/nrow(sim),
                                  SD = sd(sim$SRI))
  }
  
  sim_values <- rbindlist(sim_values)
  sim_values$CV_P <- ifelse(observedCV < sim_values$CV, 1, 0)
  sim_values$prop_NZ_P <- ifelse(observed_prop_NZ > sim_values$prop_NZ, 1, 0)
  sim_values$SD_P <- ifelse(ObservedSD < sim_values$SD, 1, 0)
  
  observed_overall <- data.table(CV = observedCV,
                                 minSCV = min(sim_values$CV, na.rm = T),
                                 maxSCV = max(sim_values$CV, na.rm = T),
                                 meanSCV = mean(sim_values$CV, na.rm = T),
                                 sdSCV = sd(sim_values$CV, na.rm = T),
                                 CV_P_value = sum(sim_values$CV_P)/nrow(sim_values),
                                 non_Zero = observed_prop_NZ,
                                 minSnon_Zero = min(sim_values$prop_NZ, na.rm = T),
                                 maxSnon_Zero = max(sim_values$prop_NZ, na.rm = T),
                                 meanSnon_Zero = mean(sim_values$prop_NZ, na.rm = T),
                                 sdSnon_Zero = sd(sim_values$prop_NZ, na.rm = T),
                                 non_Zero_P_value = sum(sim_values$prop_NZ_P)/nrow(sim_values),
                                 observedSD = ObservedSD,
                                 minSsd = min(sim_values$SD, na.rm = T),
                                 maxSsd = max(sim_values$SD, na.rm = T),
                                 SD_P_value = sum(sim_values$SD_P)/nrow(sim_values))
  
  ##############################################################################################
  # Duplicate the above, but this time for each month observed
  #Years <- unique(observed$Year)
  return(observed_overall)
  
}

# Run the above function and determine what is occurring
sim_P_value(observed = output2, simulated = sim_SRI)