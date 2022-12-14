# Author: Cameron J Baker
# 2022-12-12
# Email: cameron.baker@uqconnect.edu.au/cameron.baker@cdu.edu.au
# Code for the analyses contained within the manuscript "Long-term tracking reveals a dynamic crocodylian social system" published in
#  the journal Animal Behaviour
# Aim: This script examines the social organisation of crocodiles following the "Defining social associations" and "Crocodile social organisation" method sections

library(data.table)
library(lubridate)
library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(MuMIn)
library(performance)
library(effects)
library(tidyr)
library(VTrack)
library(ungeviz)
library(readr)

######################################################################################################################
# Create the functon to convert from raw co-occurrences into the full associations
######################################################################################################################
coooccurence_events <- function(p, data = data, time_window = 240){
  
  #################################################################
  # p: Value to select which specific receiver to determine the co-occurrences from
  # data: The output from the coAcoustic function
  
  # First subset out the data from the receiver that is of interest
  reciever <- data[Location == stations[p]]
  
  individuals <- sort(unique(reciever$Indiv1)) # create the list of individuals to sort through
  
  events <- NULL # create the list to store the results into
  
  collate_events <- function(r){
    
    # First subset out the data of the individuals
    
    whats <- reciever[Indiv1 == individuals[r]]
    # At least one receiver had a recording exactly at midnight which just left the date. Go through and find any values for the date columns that
    # only have 10 characters, as these are at midnight and then add in the time to them
    whats$DatetimeID1 <- ifelse(nchar(whats$DatetimeID1) == 10, paste0(whats$DatetimeID1," 00:00:00"), whats$DatetimeID1)
    whats$DatetimeID2 <- ifelse(nchar(whats$DatetimeID2) == 10, paste0(whats$DatetimeID2," 00:00:00"), whats$DatetimeID2)
    start <- 0
    end <- NA
    No_Indivs <- NA
    No_detects <- NA
    IDs <- NA
    f <- 1 # a value to use to make sure that the data from above is being stored correctly
    times <- NA
    
    # Create a repeat loop to go through and group the data into co-occurence events
    # Before the loop subset out the codetection data to only the individual that you are interested in (done above in the for loop)
    while(nrow(whats) > 0) {
      #while(f <= 25) {
      
      # First step, break the data into the obvious cooccurence events
      t <- 1
      for(i in 1:nrow(whats)){
        if(i == 1){
          whats$breaker[i] <- t
        } else if (as.POSIXct(whats$DatetimeID1[i]) - time_window > whats$DatetimeID2[i - 1]){
          
          t <- t + 1
          whats$breaker[i] <- t
          
        } else{
          whats$breaker[i] <- t
        }
      }
      # Second step, subset out the first break in the list
      event <- whats[breaker == 1]
      
      # If there is only two individuals present extract the start and end date
      if(length(unique(c(event$Indiv1, event$Indiv2))) == 2){
        start[f] <- ifelse(event[1,3] > event[1,4], (event$DatetimeID1[1]), (event$DatetimeID2[1])) # cooccurences start when the second individual first arrives
        end[f] <- ifelse(event[nrow(event),3] < event[nrow(event),4], as.POSIXct(as.POSIXct(event$DatetimeID1[nrow(event)]) + time_window, origin = "1970-01-01"), as.POSIXct(as.POSIXct(event$DatetimeID2[nrow(event)]) +time_window, origin = "1970-01-01"))
        #as.POSIXct(pmin(event[nrow(event),3], event[nrow(event),4])) + time_window # The end of the event is recorded as the last time point that another codetection could have occured but didn't
        No_Indivs[f] <- 2
        IDs[f] <- paste0(event$Indiv1[1],"-",event$Indiv2[1])
        No_detects[f]<- length(unique(c(event$DatetimeID1, event$DatetimeID2)))
        # Subset out the section from the the individuals data so that the breaks can be reassigned and run again. 
        # As above, make sure to set the subset of the minimum time used to prevent the function becomming stuck indefinently
        #if((as.POSIXct(end[f], origin = "1970-01-01")-time_window) < start[f]){
        #   event1 <- whats[DatetimeID1 < (as.POSIXct(end[f], origin = "1970-01-01")-time_window) & DatetimeID2 < (as.POSIXct(end[f], origin = "1970-01-01")-time_window)]
        # } else {
        #   event1 <- whats[DatetimeID1 < start[f] & DatetimeID2 < start[f]]
        #}
        #event1 <- whats[DatetimeID1 < start[f] & DatetimeID2 < start[f]] # susbet down to all of the data before the start of the event
        #event2 <- whats[DatetimeID1 > (as.POSIXct(end[f], origin = "1970-01-01")) & DatetimeID2 > (as.POSIXct(end[f], origin = "1970-01-01"))] # susbet out any data after the event ends
        
        whats <- whats[breaker != 1] # combine the two datasets and write over the original to remove this lot from the analysis
        f <- f + 1
        
      } else { 
        # For events with more than a single conspecific present, first step susbet out the detections from each conspecifc
        # and place them within a list
        consp <- list()
        starts <- NA
        ends <- NA
        for(y in 1:length(unique(event$Indiv2))){
          indiv2 <- event[Indiv2 == unique(event$Indiv2)[y]]
          # Next while in this loop calculate the time difference between each point and the one following it. If that difference is 
          # greater than 8 minutes, that is the end of the overall cooccurence event
          for(g in 1:nrow(indiv2)){
            if(g == nrow(indiv2)){
              indiv2$Time_diff[g] <- 0
            } else {
              indiv2$Time_diff[g] <- difftime(indiv2$DatetimeID2[g+1], indiv2$DatetimeID2[g], units = "mins")
            }
          }
          consp[[y]] <- indiv2
          starts[[y]] <- indiv2$DatetimeID2[1] # extract the first detection of each individual
          ends[[y]] <- indiv2$DatetimeID2[nrow(indiv2)] # extract the last detection of each individual
        }
        
        # Add in the start and end time of the focal individual
        starts[[length(starts)+1]] <- event$DatetimeID1[1] # extract the focal individuals first detection
        ends[[length(ends)+1]] <- event$DatetimeID1[nrow(event)] # extract the last detection of the focal individual
        
        # Next step, determine the start time of the largest group present by determining the time that the last individual
        # arrived
        
        start[f] <- max(starts) # the start point of the cooccurence is when the last individual first arrives
        #----------------------------------------------------------------------------------------------------------------#
        # Next determine the end period of the cooccurence event. The end time will be defined as either the time (+4min)
        # when an individual first leaves. Alternatively, if an individal leaves and then returns within the same event,
        # the end time is when there is a data gap of greater than 8mins between any individuals recordings.
        # For this, first rbind all of the dataframes of the co-detecting individuals
        consp2 <- as.data.table(do.call(rbind.data.frame, consp))
        consp2 <- consp2[order(DatetimeID2),] # reorder it so that detections are in order based on the second individuals
        consp2 <- consp2[Time_diff > 8] # subset all points where the time difference is greater than 8 mins between an individuals consecutive detections
        #------------------------------------------------------------------------------------------------------#
        if(nrow(consp2) == 1){
          # if there is only one row of points with a greater than 8 min difference, assign it as the end
          end[f] <- as.POSIXct(as.POSIXct(consp2$DatetimeID2[1]) + time_window, origin = "1970-01-01")
          # Then check to see if it is greater than the start time, if not subset the event data to everything <= to it and
          # reassign the start time. This is to prevent the overestimation of the duration of triads and above that may occur
          # if an individual has a substanital time (i.e., 40 mins) between detections and another individual is present
          # and disappears within that time gap
          if(as.POSIXct(end[f], origin = "1970-01-01") < as.POSIXct(start[f], origin = "1970-01-01")){
            event <- event[DatetimeID2 <= as.POSIXct(end[f], origin = "1970-01-01")]
            
            # Repeat the above steps to redefine the start of the assosication even
            starts <- NA
            ends <- NA
            for(y in 1:length(unique(event$Indiv2))){
              indiv2 <- event[Indiv2 == unique(event$Indiv2)[y]]
              starts[[y]] <- indiv2$DatetimeID2[1] # extract the first detection of each individual
              ends[[y]] <- indiv2$DatetimeID2[nrow(indiv2)] # extract the last detection of each individual
            }
            # Add in the start and end time of the focal individual
            starts[[length(starts)+1]] <- event$DatetimeID1[1] # extract the focal individuals first detection
            ends[[length(ends)+1]] <- event$DatetimeID1[nrow(event)] # extract the last detection of the focal individual
            start[f] <- max(starts) # the start point of the cooccurence is when the last individual first arrives
          }
          #-------------------------------------------------------------------------------------------------#
        } else if(nrow(consp2) > 1){ 
          # Encase there are multiple instances where there is a greater than 8 minute gap, repeat assigning end[f] until
          # it is greater than the start time
          # First step, check to see if there are any points greater than the current start time. IF there is not, resubset
          # the event data and reassign the starting position as above
          consp3 <- consp2[DatetimeID2 > start[f]]
          
          if(nrow(consp3) == 0){ # if there are no points greater than the start, break the data at the minimum point and
            # reassign the start time
            end[f] <- as.POSIXct(as.POSIXct(min(consp2$DatetimeID2)) + time_window, origin = "1970-01-01")
            event <- event[DatetimeID2 <= as.POSIXct(end[f], origin = "1970-01-01")] # subset to only points less than the end time
            
            # Repeat the above steps to redefine the start of the assosication even
            starts <- NA
            ends <- NA
            for(y in 1:length(unique(event$Indiv2))){
              indiv2 <- event[Indiv2 == unique(event$Indiv2)[y]]
              starts[[y]] <- indiv2$DatetimeID2[1] # extract the first detection of each individual
              ends[[y]] <- indiv2$DatetimeID2[nrow(indiv2)] # extract the last detection of each individual
            }
            # Add in the start and end time of the focal individual
            starts[[length(starts)+1]] <- event$DatetimeID1[1] # extract the focal individuals first detection
            ends[[length(ends)+1]] <- event$DatetimeID1[nrow(event)] # extract the last detection of the focal individual
            start[f] <- max(starts) # the start point of the cooccurence is when the last individual first arrives
            #------------------------------------------------------------------------------------# 
          } else {
            
            # Search through consp3 to find the miniumum point after the start time
            end[f] <- as.POSIXct(as.POSIXct(min(consp3$DatetimeID2)) + time_window, origin = "1970-01-01")
          }
          #----------------------------------------------------------------------------#
        } else {
          # If the minimum end time is greater than the start, assign it as ther end point
          if(as.POSIXct(min(ends)) > start[f]){
            end[f] <- as.POSIXct(as.POSIXct(min(ends)) + time_window, origin = "1970-01-01")
          } else { 
            # if the minimum end point is less than the start time, subset the end points so that only times equal to or
            # greater than the start time are present and then choose the minimum value of the subset
            ends <- subset(ends, ends >= start[f]) # subset the ends to only values equal to or greater than the start
            end[f] <- as.POSIXct(as.POSIXct(min(ends)) + time_window, origin = "1970-01-01") # assign the minimum value of the subset
          }
        }
        #----------------------------------------------------------------------------------------------------------------#
        # As the detection that is 8 minutes from the last point could potentially be after the first individual left,
        # this if loop will replace it with the true min value if appropriate
        
        if(as.POSIXct(end[f], origin = "1970-01-01") > as.POSIXct(as.POSIXct(min(ends)) + time_window, origin = "1970-01-01")){
          end[f] <- as.POSIXct(as.POSIXct(min(ends)) + time_window, origin = "1970-01-01")
        }
        #----------------------------------------------------------------------------------------------------------------#
        # Now that the co-detction event with the most individuals present has been identified, subset it out from the
        # data to then determine the other metrics that it included
        # As 4 minutes are added to the last position to account for the search window used, the last detection of an
        # individual may have been before the first detection of the last animal. As such this if statement will decide 
        # which is the smallest and subset the start by that to ensure that all points are included and subsetted from the
        # dataset before continuing
        
        if((as.POSIXct(end[f], origin = "1970-01-01")-time_window) < start[f]){
          occurence <- event[DatetimeID1 >= (as.POSIXct(end[f], origin = "1970-01-01")-time_window) | DatetimeID2 >= (as.POSIXct(end[f], origin = "1970-01-01")-time_window)]
        } else {
          occurence <- event[DatetimeID1 >= start[f] | DatetimeID2 >= start[f]]
        }
        occurence <- occurence[DatetimeID1 < as.POSIXct(end[f], origin = "1970-01-01") & DatetimeID2 < as.POSIXct(end[f], origin = "1970-01-01")]
        
        # If the subset results in 0 this indicates that it is the last run present and that event = occurence
        if(nrow(occurence) == 0){
          occurence <- event
        }
        # A quick subset to see if the start time is within occurence, if it is not add one to the # of detections
        
        ids <- sort(unique(c(occurence$Indiv1, occurence$Indiv2)))
        
        No_Indivs[f] <- length(ids)
        IDs[f] <- paste(ids, collapse = '-')
        
        # Subset out the section from the the individuals data so that the breaks can be reassigned and run again. 
        # As above, make sure to set the subset of the minimum time used to prevent the function becomming stuck indefinently
        
        whats <- whats[DatetimeID2 < occurence$DatetimeID2[1] | DatetimeID2 > occurence$DatetimeID2[nrow(occurence)]]
        
        
        #if((as.POSIXct(end[f], origin = "1970-01-01")-time_window) < start[f]){
        # event1 <- whats[DatetimeID1 < (as.POSIXct(end[f], origin = "1970-01-01")-time_window) & DatetimeID2 < (as.POSIXct(end[f], origin = "1970-01-01")-time_window)]
        #} else {
        # event1 <- whats[DatetimeID1 < start[f] & DatetimeID2 < start[f]]
        #}
        #event1 <- whats[DatetimeID1 < start[f] & DatetimeID2 < start[f]] # susbet down to all of the data before the start of the event
        #event2 <- whats[DatetimeID1 > (as.POSIXct(end[f], origin = "1970-01-01")) & DatetimeID2 > (as.POSIXct(end[f], origin = "1970-01-01"))] # susbet out any data after the event ends
        #whats <- rbind(event1, event2) # combine the two datasets and write over the original to remove this lot from the analysis
        f <- f + 1
      }
    }
    
    
    
    # combine all of the data into a single dataframe
    output <-  data.table(Start = start,
                          End = as.POSIXct(end,origin = "1970-01-01"),
                          #Duration = difftime(end, start),
                          No_Indivs = No_Indivs,
                          IDs = IDs)#,
    return(output)
  }
  
  registerDoParallel(10) # set the number of cores to use
  
  events <- foreach(r = 1:(length(individuals)), 
                    .packages = c("data.table","lubridate","plyr")) %dopar%
    collate_events(r)
  
  stopImplicitCluster() # stop the cluster to return memory and resources to other systems
  
  events2 <- rbindlist(events) # combine each of the individuals events into one dataframe
  events <- events2[order(Start, -No_Indivs),] # reorder the dataframe so that it is in temporal order and by the number of individuals present
  
  events$Start <- lubridate::ymd_hms(events$Start, tz = "Australia/Brisbane")
  events$End <- lubridate::ymd_hms(events$End, tz = "Australia/Brisbane")
  
  # Create a new for loop to run through and double check that the start and end times are correct, while also deleting
  # any rows that are occurring at the same time as either the time point before or after it (i.e. a dyad being formed
  # while a triad is occurring)
  
  t <- 2 # a vector to control and stop the function
  # First step, remove any co-occurrences that have/are occurring at the same time as the previous one
  while(t <= nrow(events)-1) { # start at the second point as the first does not need checking
    
    # First step, subset out any cooccurence events that happen within a larger grouping (i.e., removing the dyad that has
    # formed within a concurrent triad)
    if(events$Start[t] >= events$Start[t-1] & events$End[t] <= events$End[t-1]){ # if the duration of this events occurs within the previous one, delete the row
      events <- events[-t]
    } 
    
    if(t <= (nrow(events)-1)){
      
      repeat{
        
        if(t == nrow(events)){
          break
        } else if (events$Start[t] <= events$Start[t+1] & events$End[t] >= events$End[t+1]){ # if the following row is within the duraton of the current one, delete it
          events <- events[-(t+1)]
        } else {
          break
        }
      }
      #if (events$Start[t] <= events$Start[t+1] & events$End[t] >= events$End[t+1]){ # if the following row is within the duraton of the current one, delete it
      #events <- events[-(t+1)]
      #}
    }
    
    t <- t + 1
  }
  # As it is possible that some of the end times for some durations is after the start of the next event because of the
  # plus 4 minutes at the end, if the end time is greater than the start of the next replace the end time with the start
  # For this to be accurate though, the direction of the shift will be decided by whether the group is gaining or losing
  # members. If it gains, then the new end time is the start of the next group. If they lose a member though, the start time
  # is when that individual has left
  t <- 2
  while(t <= nrow(events)){
    
    if(events$No_Indivs[t] - events$No_Indivs[t-1] > 0){ # if the group adds another individual, the previous endtime is this ones start
      if(events$Start[t] < events$End[t-1]){
        events$End[t-1] <- events$Start[t]
      }
    } else if(events$No_Indivs[t] - events$No_Indivs[t-1] < 0){ # if the group loses an individual, the current starttime is the previous endtime
      if(events$Start[t] < events$End[t-1]){
        events$Start[t] <- as.POSIXct(events$End[t-1], origin = "1970-01-01", tz = "Australia/Brisbane")
      }
    }
    #if(events$Start[t] < events$End[t-1]){
    # events$End[t-1] <- events$Start[t]
    #}
    t <- t + 1
  }
  # Calculate the duration of the coocurrence event
  events$Duration <- as.numeric(difftime(as.POSIXct(events$End,origin = "1970-01-01"), as.POSIXct(events$Start,origin = "1970-01-01"), units = "mins"))
  
  events <- events[Duration >= 0] # remove any erroneous occurrences that may have been made
  
  # Add in the location that these cooccurences occurred at
  events$Location <- reciever$Location[1]
  
  # Amazingly, legitimate durations of 0 minutes have been captured (twice!) in the data at Island midstream, so the dataset
  # should potentially be subset to remove all 0 durations
  #events <- events[Duration > 0] # r
  #print(paste0(stations[p]," has completed successfully"))
  # Write a copy to file encase it stops working at another receiver
  write.csv(events, paste0("Data/Network-analysis/Occurences/occur_",stations[p],".csv"), row.names = F)
  
  print(paste0("Co-occurence events for ", stations[p], " have been determined"))
  
  return(events)
}

######################################################################################################################
#### Import and prepare the data for analysis ####

# Import the raw data
croc_data <- fread("Data/Baker_et_al_data.csv")

# Import the output from the coAcoustic function
croc_codetect <- fread("Data/crocodile_edges_2020_4.csv")

# Create a column for the study year of each co-detection. 
# Study years run from August-July to reflect when crocodiles were initally tagged
croc_codetect$Year <- ifelse(month(croc_codetect$DatetimeID1) >= 8, year(croc_codetect$DatetimeID1), year(croc_codetect$DatetimeID1)-1)

# Create a unique list of the receivers that co-detections were observed
stations <- unique(allcrocs_dyads$Location)

# Set the number of cores to use
registerDoParallel(5) 
# Run the cooccurence_events function in parallel to speed up computing time
lmig <- foreach(p = 1:(length(stations)), 
                .packages = c("data.table","lubridate","plyr")) %dopar%
  coooccurence_events(p, data = allcrocs_dyads, time_window = 240)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# Combine all of the outputs from the cooccurence_events function into a single dataset
occur_events <- rbindlist(lmig)

write.csv(occur_events, "Data/Occurence_events.csv", row.names = F) # Write a copy to save having to rerun the above function again.

#######################################################################################################################
#### Determine the yearly reproductive status of crocodiles ####

# Import the list of crocodile names, along with the year they were first caught, their size and sex
croc_list <- fread("Data/Baker_et_al_init_croc_size.csv")

# Add a new column with the reproductive status of each indivdual for their first year in the study
croc_list$Repro_status <- ifelse(croc_list$Sex == "M" & croc_list$Total_Length < 3300, "NR", 
                                 ifelse(croc_list$Sex == "F" & croc_list$Total_Length < 2000, "NR", "R"))

croc_list$Year <- as.numeric(croc_list$Year)

# Create a list of all the years of the study
years <- seq(2008, 2021, 1)

# Create the function to increase the size of indivduals per year
reproductive_status <- function(i){
  
  # Subset out all of the crocodiles that have been captured until the year of interest
  crocodiles <- croc_list[Year <= years[i]]
  # For each crocodile  add in their expected growth rate multiplied by the number of years since they were first captured
  crocodiles$Total_Length <- crocodiles$Total_Length + 73.9 * (years[i] - crocodiles$Year) 
  
  crocodiles$Repro_status <- ifelse(crocodiles$Sex == "M" & crocodiles$Total_Length < 3300, "NR", 
                                    ifelse(crocodiles$Sex == "F" & crocodiles$Total_Length < 2000, "NR", "R"))
  crocodiles$Year <- years[i]
  return(crocodiles)
}

# Run the function for each of the years the study has been going for
dynamic_repro_status <- lapply(1:length(years),reproductive_status)

# Combine the output of the reproductive _status function into a single dataframe
dynamic_repro_status <- rbindlist(dynamic_repro_status)

# import the crocodile metrics file to get the transmitterid of each animal and then join that with the dynamic reproductive statuses of each individual
croc_metric <- fread("Data/Metrics/Baker_et_al_crocodile_metrics.csv") 
croc_metric <- croc_metric[,c(1,3)] # remove all of the unnecessary columns

# Join the two datasets together
dynamic_repro_status <- plyr::join(dynamic_repro_status, croc_metric, by = "Croc_Name")

# Reorder the columns and remove any unnecessary ones
dynamic_repro_status <- dynamic_repro_status[,c(6,1,3,4,5)]

# Export the data for later
write.csv(dynamic_repro_status, "Data/Dynamic_repro_status.csv", row.names = F)

#####################################################################################################################################################
### Determine the reproductive status combination for each of the different co-occurrence events ####

# This is a function that determines the reproductive status combination of each co-occurrence for each year to allow for the reproductive status
# of conspecifics to dynamically adjust as they grow
yearly_repro <- function(i){
  
  # first subset the co-occurrence data and the reproductive statuses to the year of interest
  yearly_co_occur <- co_occures[Year == unique(co_occures$Year)[i]]
  repro_year <- dynamic_repro_status[Year == unique(co_occures$Year)[i]]
  
  # extract a list of reproductive and non-reproductive individuals
  repro <- unique(repro_year[Repro_status == "R"]$TRANSMITTERID)
  non_repro <- unique(repro_year[Repro_status != "R"]$TRANSMITTERID)
  
  yearly_co_occur$Repro <- NA
  yearly_co_occur$Non_repro <- NA
  
  for(q in 1:nrow(yearly_co_occur)){
    # If there is at least one non-reproductive individual assign a yes to the non_repro column
    if(length(intersect(unlist(yearly_co_occur$Individuals[q]), non_repro)) >= 1){
      yearly_co_occur$Non_repro[q] <- "Yes"
    } else {
      yearly_co_occur$Non_repro[q] <- "No"
    }
    # If there is at least one reproductive individual assign a yes to the repro column
    if(length(intersect(unlist(yearly_co_occur$Individuals[q]), repro)) >= 1){
      yearly_co_occur$Repro[q] <- "Yes"
    } else {
      yearly_co_occur$Repro[q] <- "No"
    }
  }
  
  # Creates a column with the Maturity_status structure of the groupings. A = Adult, j = Juvenile
  yearly_co_occur$Maturity_status <- ifelse(yearly_co_occur$Repro == "Yes" & yearly_co_occur$Non_repro == "Yes", "AJ", 
                                ifelse(yearly_co_occur$Repro == "Yes" & yearly_co_occur$Non_repro == "No", "AA", "JJ"))
  
  yearly_co_occur$Repro <- yearly_co_occur$Non_repro <- NULL
  
  return(yearly_co_occur)
}

# Run the function across each of the main study years
co_occures2 <- lapply(1:11, yearly_repro)
# Combine the output into a single dataframe
co_occures2 <- rbindlist(co_occures2)

# Assign the sex combination of each association event
co_occures2$Sex <- factor(ifelse(co_occures2$Sex_Ratio == 1, "MM", ifelse(co_occures2$Sex_Ratio == 0, "FF", "MF")))

#####################################################################################################################################################
### Determine the number of detections, number of crocodiles observed and the number of associations present at receivers for each month of the study ####
associations <- function(i){
  
  # subset out the receiver of intrest from the raw data and the co-occurrences
  receiver <- crocs[STATIONNAME == receivers[i]]
  receiver_co_occur <- co_occures2[Location == receivers[i]]
  # output <- list()
  # Next determine the years that this receiver was present within the array
  years <- unique(receiver_co_occur$Year)
  yearly_values <- list()
  # For each year
  for(t in 1:length(years)){
    # Subset out the year of interest
    year_detect <- receiver[Year == years[t]]
    year_co_occur <- receiver_co_occur[Year == years[t]]
    
    # Now figure out all of the months that this receiver had detections recorded for within this month
    months <- unique(year_co_occur$Month)
    monthly_values <- list()
    # for each of the months, determine the number of co-occurrences by association type
    for(r in 1:length(months)){
      month_detect <- year_detect[Month == months[r]]
      month_co_occur <- year_co_occur[Month == months[r]] 
      
      combos <- c("AA", "AJ","JJ")
      values <- list() # create a list to store vectors in
      for(q in 1:length(combos)){
        co_occur <- month_co_occur[Maturity_status == combos[q]]
        if(nrow(co_occur) > 0){
          values[[q]] <- data.table(Maturity_status = combos[q],
                                    Year = years[t],
                                    Month = months[r],
                                    Location = receivers[i],
                                    Co_occur = nrow(co_occur),
                                    Number_crocs = length(unique(month_detect$TRANSMITTERID)),
                                    Number_detections = nrow(month_detect))
        } else {
          values[[q]] <- data.table(Maturity_status = combos[q],
                                    Year = years[t],
                                    Month = months[r],
                                    Location = receivers[i],
                                    Co_occur = 0,
                                    Number_crocs = length(unique(month_detect$TRANSMITTERID)),
                                    Number_detections = nrow(month_detect))
        }
      }
      monthly_values[[r]] <- rbindlist(values)
      
    }
    
    yearly_values[[t]] <- rbindlist(monthly_values)
  }
  
  return(rbindlist(yearly_values))
}

# Create a list of receivers that associations were detected at
receivers <- unique(co_occures2$Location)

# Run the above function for each receiver
number_co_occures <- lapply(1:length(receivers), associations)
# Combine the output into a single data frame
number_co_occures <- rbindlist(number_co_occures)

#####################################################################################################################################################
### Does time of year influence the number of detections at acoustic receivers? ####

# Create the model
num_detect_mod <- glmer(Number_detections ~ Month + (1|Location) + (1|Year), family = "poisson",  data = number_co_occures)

# Check model performance using the check-model function from performance
check_model(num_detect_mod)
# Overall model looks fine, check the models summary
summary(num_detect_mod2)
# Everything looks good and there appears to be an impact of month on the number of detections
car::Anova(num_detect_mod) # There is a significant effect present
# What months are significantly different from each other?
emmeans(num_detect_mod,  pairwise ~ Month)
# Generate and save the predictions from the model for creating Figure 2
month_pred_detect <- as.data.table(effects::predictorEffect("Month", num_detect_mod))
write.csv(month_pred_detect, "Data/Model_predictions/Month_detects.csv", row.names = F)

#####################################################################################################################################################
### Does time of year influence the number of crocodiles observed at acoustic receivers? ####

# Create the model
num_crocs_mod <- glmer(Number_crocs ~ Month + (1|Location) + (1|Year), family = "poisson",  data = number_co_occures)
# Model did not fully converge. Following the instructions from Ben Bolker (link below), restart the model and increase the number of interations
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
ss <- getME(num_crocs_mod ,c("theta","fixef"))
num_crocs_mod <- update(num_crocs_mod ,start=ss,control=glmerControl(optCtrl=list(maxfun=2e5))) # The model converged
# Check model diagnostics
check_model(num_crocs_mod)
# Overall model looks fine, check the models summary
summary(num_crocs_mod)
# Everything looks good and there appears to be an impact of time of year
car::Anova(num_crocs_mod) # There is a significant effect
# What months are significantly different from each other?
emmeans(num_crocs_mod,  pairwise ~ Month)
# Generate and save the predictions from the model for creating Figure 2
month_pred_crocs <- as.data.table(effects::predictorEffect("Month", num_crocs_mod))
write.csv(month_pred_crocs, "Data/Model_predictions/Month_crocodiles.csv", row.names = F)

#####################################################################################################################################################
### Does time of year and maturity status influence the number of associations observed at acoustic receivers? ####

# Create the model
co_occur_mod <- glmer(Co_occur ~ Number_crocs + Month +  (1|Location) + (1|Year), family = "poisson",  data = whats)
# Model failed to converge. Model trouble shooting completed follow the instructions from Ben Bolker: # https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
tt <- getME(co_occur_mod ,"theta")
ll <- getME(co_occur_mod ,"lower")
min(tt[ll==0])
# No singularity, now check the gradient calculations
derivs1 <- co_occur_mod @optinfo$derivs
sc_grad1 <- with(co_occur_mod ,solve(Hessian,gradient))
max(abs(sc_grad1))
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
# All of the calculations seem fine, try restarting the model and have it run for more iterations
ss <- getME(co_occur_mod ,c("theta","fixef"))
co_occur_mod2 <- update(co_occur_mod ,start=ss,control=glmerControl(optCtrl=list(maxfun=2e5)))
# Still no success, try with a different optimizer at the end
co_occur_mod2 <- update(co_occur_mod ,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # Success!
# Check model diagnostics
check_model(co_occur_mod2)
# Overall model looks fine, check the models summary
summary(co_occur_mod2)
# Everything looks good and there appears to be an impact of time of year
car::Anova(co_occur_mod2) # There is a significant interaction between time of year and maturity status
# Generate and save the predictions from the model for creating Figure 4
month_pred <- as.data.table(effects::predictorEffect("Month", co_occur_mod2))
write.csv(month_pred, "Data/Model_predictions/Month_co_occurs.csv", row.names = F)

#####################################################################################################################################################
### Does time of year and maturity status influence the proportion of time indivdiuals associate with conspecifics? ####

# First determine the total time that individuals were observed in the acoustic array
residencies <- RunResidenceExtraction(sInputFile = croc_data, sLocation = "STATIONNAME", iResidenceThreshold = 1, iTimeThreshold = 240)

# Create a new column with the IDs of each conspecific present in the association dataset
occur_events$Individuals <- strsplit(occur_events$IDs, "-")

# For each dataset, create a Year, Month and Unique month (i.e., 2020-08) column
residencies$Month <- month(residencies$STARTTIME)
residencies$Year <- year(residencies$STARTTIME)
residencies$Unique_months <- paste0(residencies$Year, "-", residencies$Month)

co_occures$Month <- month(co_occures$Start)
co_occures$Year <- year(co_occures$Start)
co_occures$Unique_months <- paste0(co_occures$Year, "-", co_occures$Month)

#-----------------------------------------------------------------------------------------------#
# Create the function which determines the yearly or monthly proportion of time individuals are co-occurring with conspecifics
prop_time_co_occur <- function(co_occures, residencies, yearly = T){
  
  ########################################
  # co_occurres: The association dataframe
  # residences: The output from the RunResidenceExtraction function
  # yearly: If true it calculates the proportion of time with conspecifis for the full year, if false it calculates it monthly
  
  if(yearly == T){
    
    years <- unique(residencies$Year)
    output <- list()
    for(i in 1:length(years)){
      
      #Subset out the year of interest in both datasets
      resid_year <- residencies[Year == years[i]]
      co_occures_year <- co_occures[Year == years[i]]
      repro_year <- dynamic_repro_status[Year == years[i]] # subset out the reproductive statuses for this year
      
      # Extract a list of the individuals present
      indivs <- unique(resid_year$TRANSMITTERID)
      
      indiv_list <- list()
      for(t in 1:length(indivs)){
        
        # Subset out all of the co-occurrences that the individual of interest was within
        co_occures_year$Present <- NA
        # an if else statement to identify the co-occurrences that have individuals who could not generate a home range present
        for(q in 1:nrow(co_occures_year)){
          if(length(setdiff(unlist(co_occures_year$Individuals[q]), indivs[t])) < co_occures_year$No_Indivs[q]){
            co_occures_year$Present[q] <- "Yes"
          } else {
            co_occures_year$Present[q] <- "No"
          }
          
        }
        
        co_occures_indiv <- co_occures_year[Present == "Yes"]
        resid_indiv <- resid_year[TRANSMITTERID == indivs[t]]
        # If the end reason is less than 4
        
        indiv_list[[t]] <- data.table(TRANSMITTERID = indivs[t],
                                      Year = years[i],
                                      Number_assoc = nrow(co_occures_indiv),
                                      Time_co_occur = sum(co_occures_indiv$Duration),
                                      Time_detected = sum(as.numeric(resid_indiv$Duration)), # convert from seconds to minutes
                                      Prop_time = sum(co_occures_indiv$Duration)/(sum(as.numeric(resid_indiv$Duration))),
                                      Numer_conspecifcs = length(unique(resid_year$TRANSMITTERID)))
        
      }
      output[[i]] <- rbindlist(indiv_list)
      
    }
    
  } else {
    
    years <- unique(residencies$Unique_months)
    output <- list()
    for(i in 1:length(years)){
      
      #Subset out the year of interest in both datasets
      resid_year <- residencies[Unique_months == years[i]]
      co_occures_year <- co_occures[Unique_months == years[i]]
      repro_year <- dynamic_repro_status[Year == substr(years[i],1,4)] # subset out the reproductive statuses for this year
      #quart_metrics <- social_metrics[Unique_months == years[i]]
      # Extract a list of the individuals present
      indivs <- unique(resid_year$TRANSMITTERID)
      
      indiv_list <- list()
      for(t in 1:length(indivs)){
        
        # Subset out all of the co-occurrences that the individual of interest was within
        co_occures_year$Present <- NA
        if(nrow(co_occures_year) > 0){
          # an if else statement to identify the co-occurrences that have individuals who could not generate a home range present
          for(q in 1:nrow(co_occures_year)){
            if(length(setdiff(unlist(co_occures_year$Individuals[q]), indivs[t])) < co_occures_year$No_Indivs[q]){
              co_occures_year$Present[q] <- "Yes"
            } else {
              co_occures_year$Present[q] <- "No"
            }
            
          }
        }
        
        co_occures_indiv <- co_occures_year[Present == "Yes"]
        resid_indiv <- resid_year[TRANSMITTERID == indivs[t]]
        # If the end reason is less than 4
        # create the list for each of the reproductive status combinations present
        repro_list <- list()
        if(nrow(co_occures_indiv) > 0){
          for(g in 1:length(unique(co_occures_indiv$Maturity_status))){
            co_occures_indiv_repro <- co_occures_indiv[Maturity_status == unique(co_occures_indiv$Maturity_status)[g]]
            repro_list[[g]] <- data.table(TRANSMITTERID = indivs[t],
                                          Maturity_status = unique(co_occures_indiv$Maturity_status)[g],
                                          Quarter = years[i],
                                          Number_assoc = nrow(co_occures_indiv_repro),
                                          Time_co_occur = ifelse(nrow(co_occures_indiv_repro) == 0, 0, sum(co_occures_indiv_repro$Duration)),
                                          Time_detected = sum(as.numeric(resid_indiv$Duration)), # convert from seconds to minutes
                                          Prop_time = sum(co_occures_indiv_repro$Duration)/(sum(as.numeric(resid_indiv$Duration))),
                                          Number_conspecifcs = length(unique(resid_year$TRANSMITTERID)))
          }
        } else {
          repro_list[[1]] <- data.table(TRANSMITTERID = indivs[t],
                                        Maturity_status = c(ifelse(repro_year[TRANSMITTERID == indivs[t]]$Repro_status == "R", "AA", "JJ"), "AJ"),
                                        Quarter = years[i],
                                        Number_assoc = 0,
                                        Time_co_occur = 0,
                                        Time_detected = sum(as.numeric(resid_indiv$Duration)), # convert from seconds to minutes
                                        Prop_time = 0,
                                        Number_conspecifcs = length(unique(resid_year$TRANSMITTERID)))
        }
        
        
        indiv_list[[t]] <- rbindlist(repro_list)
        
      }
      output[[i]] <- rbindlist(indiv_list)
      
    }
    
  }
  
  return(rbindlist(output))  
}


# Calculate the monthly proportion of time that individuals spent with conspecifics
monthly_prop_time <- prop_time_co_occur(co_occures, residencies, yearly = F)

# Due to some slight discrepencies with how the proportion of time spent at receivers compared to in groups is calculated, update the file so 
# if the proportion of time in a group is greater than the proportion in the array, that value is replaced
monthly_prop_time$Time_detected <- ifelse(monthly_prop_time$Time_co_occur > monthly_prop_time$Time_detected, monthly_prop_time$Time_co_occur,
                                          monthly_prop_time$Time_detected)
monthly_prop_time$Prop_time <- monthly_prop_time$Time_co_occur/monthly_prop_time$Time_detected
monthly_prop_time$Month <- factor(substr(monthly_prop_time$Quarter, 6,nchar(monthly_prop_time$Quarter)), levels = c("5","6","7","8","9","10","11","12","1","2", "3","4"))
monthly_prop_time$Month <- ifelse(monthly_prop_time$Month == "1", "Jan",
                                  ifelse(monthly_prop_time$Month == "2", "Feb",
                                         ifelse(monthly_prop_time$Month == "3", "Mar",
                                                ifelse(monthly_prop_time$Month == "4", "Apr",
                                                       ifelse(monthly_prop_time$Month == "5", "May",
                                                              ifelse(monthly_prop_time$Month == "6", "Jun",
                                                                     ifelse(monthly_prop_time$Month == "7", "Jul",
                                                                            ifelse(monthly_prop_time$Month == "8", "Aug",
                                                                                   ifelse(monthly_prop_time$Month == "9", "Sep",
                                                                                          ifelse(monthly_prop_time$Month == "10", "Oct",
                                                                                                 ifelse(monthly_prop_time$Month == "11", "Nov", "Dec")))))))))))
monthly_prop_time$Month <- factor(monthly_prop_time$Month, levels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","Dec", "Jan", "Feb", "Mar", "Apr"))
monthly_prop_time$Year <- factor(substr(monthly_prop_time$Quarter, 1, 4))
monthly_prop_time$TRANSMITTERID <- factor(monthly_prop_time$TRANSMITTERID)
# Create the model
prop_time_mod <-  glmer(cbind(round(Time_co_occur), round(Time_detected)) ~   Month*Age + (1|TRANSMITTERID) + (1|Year), family = "binomial", data = monthly_prop_time)
# All of the calculations seem fine, try restarting the model and have it run for more iterations
ss <- getME(prop_time_mod ,c("theta","fixef"))
prop_time_mod2 <- update(prop_time_mod ,start=ss,control=glmerControl(optCtrl=list(maxfun=2e5)))
# Still no success, try with a different optimizer at the end
prop_time_mod2 <- update(prop_time_mod ,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

check_model(prop_time_mod) # check model diagnostics
car::Anova(prop_time_mod) # Examine what correlations are present
summary(prop_time_mod)
# Generate the predictions from the model and save for later
month_dur <- as.data.table(effects::predictorEffect("Month", prop_time_mod))
write.csv(month_dur, "Data/Model_predictions/Proportion_of_time.csv")
write.csv(monthly_prop_time, "Data/Model_predictions/Prop_time_raw_data.csv") # also export a version of the raw data to make plotting easier