# Author: Cameron J Baker
# 2022-12-12
# Email: cameron.baker@uqconnect.edu.au/cameron.baker@cdu.edu.au
# Code for the analyses contained within the manuscript "Long-term tracking reveals a dynamic crocodylian social system" published in
# the journal Animal Behaviour
# Aim: This script identifies associations between tagged crocodiles following the methods described in the "Defining social associations" method section

library(data.table)
library(parallel)

#########################################################################################################################################################
### Import and prepare the required dataset ####
croc_data <- fread("Data/Baker_et_al_data.csv")

######################################################################################################################################################
### Create the functions for identifying co-occurrences within the river system ####

coAcoustic <- function(data, n = 240, icores = 2){
  
  #####################################################################################################
  # data: The acoustic telemetry data in the VTrack format you would like to determine associations from
  # n: The number of seconds of your sampling window. Default is 4 minutes
  # icores: The number of cores you would like to use for parallel processing
  ######################################################################################################
  data2 <- as.data.table(allcrocs) # create a copy of the data to work from

  # Create empty vectors to populate below
  Receiverlist <- NULL

  
  # Create the floating window function that is used to idenitfy co-occurrences of the focal individual
  floatingWindow <- function(x, y, n = 240){
    
    ###################################################################################################
    # x: The detections of the focal individual
    # y: The detections of all other crocodiles
    # n: The number of seconds of your sampling window. Default is 4 minutes
    ###################################################################################################
    
    # create the empty vectors required to store the data outputs
    indiv1 <- NULL
    indiv2 <- NULL
    DatetimeID1 <- NULL
    DatetimeID2 <- NULL
    Location <- NULL
    t <- 1
    
    for(o in 1:nrow(x)){
      if(o == 1){
        search <- y[DATETIME >= (x$DATETIME[o] - n) & DATETIME <= (x$DATETIME[o] + n)]  # subset out the floating window of 120 seconds around the point
        
        for(e in 1:nrow(search)){
          
          if(search[, .N] >= 1){
            indiv1[[t]] <- x$TRANSMITTERID[1] #focal individual
            indiv2[[t]] <- search$TRANSMITTERID[e] # Individual they are occurring with
            DatetimeID1[[t]] <- as.character(x$DATETIME[o]) # Date time that the first individual was recorded
            DatetimeID2[[t]] <- as.character(search$DATETIME[e]) # Date time that the second indvidual was recorded
            Location[[t]] <- as.character(x$STATIONNAME[e]) # and at what receiver station
            
            t <- t + 1 # add x so that it iteratively keeps recording to the next piece and doesn't overwrite anything
          }
        }
        
        
      } else if(as.numeric(difftime(x$DATETIME[o],x$DATETIME[o-1], units = "secs")) >= n * 2){ #N multiplied by 2 to avoid overlapping search windows
        
        search <- y[DATETIME >= (x$DATETIME[o] - n) & DATETIME <= (x$DATETIME[o] + n)]  # subset out the floating window of 120 seconds around the point
        
        
        for(e in 1:nrow(search)){
          
          if(search[, .N] >= 1){
            indiv1[[t]] <- x$TRANSMITTERID[1] #focal individual
            indiv2[[t]] <- search$TRANSMITTERID[e] # Individual they are occurring with
            DatetimeID1[[t]] <- as.character(x$DATETIME[o]) # Date time that it occurred 
            DatetimeID2[[t]] <- as.character(search$DATETIME[e]) # Date time that the second indvidual was recorded
            Location[[t]] <- as.character(search$STATIONNAME[e]) # and at what receiver station
            
            t <- t + 1 # add one so that it iteratively keeps recording to the next piece and doesn't overwrite anything
          }
        }
        
        
      } else {  # Use the time plus 2 minutes of the previous level to make sure it is only searching fo co-occurences in a time window that has not been previously searched to prevent over inflation of occurences
        search <- y[DATETIME >= (x$DATETIME[o-1] + n +1) & DATETIME <= (x$DATETIME[o] + n)]  # subset out the floating window of 120 seconds around the point
        
        for(e in 1:nrow(search)){
          
          if(search[, .N] >= 1){
            indiv1[[t]] <- x$TRANSMITTERID[1] #focal individual
            indiv2[[t]] <- search$TRANSMITTERID[e] # Individual they are occurring with
            DatetimeID1[[t]] <- as.character(x$DATETIME[o]) # Date time that it occurred 
            DatetimeID2[[t]] <- as.character(search$DATETIME[e]) # Date time that the second indvidual was recorded
            Location[[t]] <- as.character(search$STATIONNAME[e]) # and at what receiver station
            
            t <- t + 1 # add one so that it iteratively keeps recording to the next piece and doesn't overwrite anything
          }
        }
      }
    }
    
    
    output <- data.frame(Indiv1 = indiv1,
                         Indiv2 = indiv2,
                         DatetimeID1 = DatetimeID1,
                         DatetimeID2 = DatetimeID2,
                         Location = Location)
    
    return(output)
    
  }
  
  # Function to parallelize the function by receiver stations
  ReceiveFloatingWindow <- function(i){
    
    ###################################################################################################
    # i: Value to select which specific receiver to determine the co-occurrences from
    ###################################################################################################
    
    p <- 1
    
    # Create an empty list to store each individuals observed co-occurrences within
    IndivList <- NULL
    
    receive <- data2[STATIONNAME == receivers[i]]
    
    individuals <- sort(unique(receive$TRANSMITTERID), decreasing = F) # creates a list of all the individuals detected and then sort in decreasing order to make sure each receiver runs through in the same order
    
    if(length(individuals) > 1){ # if there is more than one individual present, search for co-occurrences
      
      for(d in 1:length(individuals)){ # for each of the individuals present, search to see if there are any co-occurring conspecifics
        
        indiv <- receive[TRANSMITTERID == individuals[d]]
        
        remaining <- unique(receive$TRANSMITTERID)
        
        receive <- receive[TRANSMITTERID != individuals[d]] # subset out the individual to avoid duplicating co-occurrences and speed up code.
        
        if(length(remaining) <= 1){
          break
        } else{
          
          IndivList[[p]] <- floatingWindow(indiv, receive, n)
          
          p <- p + 1
        }
      }
      output <- do.call(rbind.data.frame, IndivList)
      return(output)
    }
  }
  #----------------------------------------------------------------------------------------------------------------#
  # Subset down and run the below code for each receiever station to avoid the code assigning co-occurences when individuals are at different receivers
  receivers <- unique(data2$STATIONNAME)

  # Initiate cluster for parallel processing 
  cl <- makeCluster(icores)
  
  # To ensure packages, functions and data are sent to the workers, export the names directly
  clusterExport(cl, list("ReceiveFloatingWindow","n", "floatingWindow", "data.table", "data2", "receivers"),
                envir = environment()) # makes it look in the functions environment to find the things
  
  # Run the parallel version of lapply
  Receiverlist <- parLapply(cl, 1:length(receivers), ReceiveFloatingWindow)
  
  # Combine all of the lists created for each individual per receiver station together
  output <- do.call(rbind.data.frame, Receiverlist)
  output$Weight <- 1 # assign a column with ones to allow weightings of intereactions to be determined
  
  # Code seems to be passing Falses between trues as NA currently, this removes them from the output 
  output <- na.omit(output) 
  
  #output <- aggregate(x = output[c("Weight")],
  #by = output[c("Indiv1","Indiv2","Datetime", "Location")],
  #FUN =  sum)
  # stopImplicitCluster()
  stopCluster(cl) # stop the cluster to return memory and resources to other systems
  
  return(output)
}

######################################################################################################################################################
### Identify the co-occurrences present throughout the dataset ####

interact <- coAcoustic(data = croc_data, n = 240, icores = 6)

# Save the output for future analyses
write.csv(interact, "Data/crocodile_edges_2020_4.csv", row.names = F)
