# Author: Cameron J Baker
# 2022-12-12
# Email: cameron.baker@uqconnect.edu.au/cameron.baker@cdu.edu.au
# Code for the analyses contained within the manuscript "Long-term tracking reveals a dynamic crocodylian social system" published in
# the journal Animal Behaviour
# Aim: This script examines whether the associations of crocodiles are stable through time following the "Stability of associations through time" methods section

library(plyr)
library(data.table)
library(lubridate)
library(doParallel)
library(foreach)
library(tibble)
library(asnipe)
library(ggplot2)
library(FinCal)
library(ggpubr)

######################################################################################################################
#### Import and prepare the data for analysis ####
raw_dyads <- read("Data/crocodile_edges_2020_4.csv")

# Create a column within the dataset, which places the data into the correct years for the dynamic networks
raw_dyads$Year <- substr(raw_dyads$DatetimeID1, 1, 7)
raw_dyads$Month <- month(raw_dyads$DatetimeID1)

years <- sort(unique(raw_dyads$Year))
years1 <- sort(unique(raw_dyads$Year))# a vector for importing the simulated data
#===================================================================================================================================================#
# Home ranges were constructed following the methods described in Baker et. al (2022) Crocodile social environments dictated by male philopatry. Behavioral Ecology, 33(1), 156-166. doi: https://doi.org/10.1093/beheco/arab120
# The required R code and shape files can be found at:  https://github.com/Cameron-J-Baker/Supporting-code-Baker_et-al 
# Import the observed home ranges overlap matrices produced following the code above
# Next import the home range overlaps of the population and then subset out the individuals each month who didn't have the data to generate a home
# range
files <- list.files(path = "Data/Overlap_matrices/Monthly_Wenlock/", full.names = T) # create a vector list with the names of all files

lmig <- lapply(files, fread) # Import all files using lapply

# Convert each of the files imported into matrices and move the TRANSMITTERIDs into the rownames
for(i in 1:length(lmig)){
  lmig[[i]] <- as.data.frame(lmig[[i]]) # convert to data.frame
  colnames(lmig[[i]]) = lmig[[i]][1, ] # convert the first row back to column names
  lmig[[i]] = lmig[[i]][-1, ] # remove the first row with the names  
  rownames(lmig[[i]]) <- c() # remove the current row names
  lmig[[i]] <- tibble::column_to_rownames(lmig[[i]], var = "NA") # convert iD back into row.names
  
  lmig[[i]] <- as.matrix(lmig[[i]])
  
  lmig[[i]][upper.tri(lmig[[i]], diag = T)] <- NA 
  # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
  lmig[[i]] <- na.omit(as.data.table(as.table(lmig[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
  names(lmig[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns
  #lmig[[i]]$Month <- substr(files[i], 44,45)
  lmig[[i]]$Year <- substr(files[i], 39,45)
  
}

yearlyHR <- rbindlist(lmig)

#=================================================================================================================================================#
# Import the prepAssocIndex function to calculate the SRI of each dyad
source("R_Functions/prepAssocIndex.R")
#=================================================================================================================================================#
# Subset out the individuals who did not have the required data to generate a home range and then determine the SRI of each dyad

yearly_SRI <- list()
yearSRI<- function(i){
  year_sri <- raw_dyads[Year == years[i]] # subset out all of the co-detections for that year
  year_HR <- yearlyHR[Year == years[i]] # extract the home range overlaps for the year of interest
  individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
  year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
  if(nrow(year_sri) > 0){
    output <- prepAssocIndex(year_sri)
    output <- output[,c(1,2,7)] # reduce down to just the first, second and seventh column for plotting
    output$Month <- years[i] 
    return(output) 
  }
}

what <- yearSRI(6)

yearly_SRI <- lapply(1:length(years), yearSRI)

##################################################################################################################################################
# Create a function to convert the dataframe of SRIs to a matrix

dataframe_Matrix <- function(r, data_list, totalIDs){
  
  ########################################################################
  # r: Vector to indicate which of the unique months to be selected from the list
  # data_list: The list object containing the SRI dataframes for each month
  # totalIDs: An object containing the IDs of each crocodile observed associating throughout the study
  #########################################################################
  
  # Extract the dataframe of interest from the list
  dataframe <- data_list[Year == years[r]]
  dataframe$Year <- NULL # Now remove this column so the function below can run without problem
  # If a month has no associations recorded, create and save an empty matrix of all individuals, else run the function as normal
  if(nrow(dataframe) == 0){
    
    res2 <- as.data.frame(matrix(NA, ncol = length(totalIDs), nrow = length(totalIDs)))
    # Add the column and row names to the matrix to make it easier to work with
    rownames(res2) <- totalIDs # add the ID labels to the rows
    colnames(res2) <- totalIDs # add the ID labels to the columns
    
    res2 <- rownames_to_column(res2, var="ID") # sets the rownames to a column
    res2 <- res2 %>% arrange(ID) # rearrange the rows so that they are in ascending order
    res2 <- res2[ , order(names(res2))] # place the columns into ascending order
    res2 <- column_to_rownames(res2, var = "ID") # convert ID back into row.names
    
    res2 <- data.matrix(res2, rownames.force = "ID") # convert back into a matrix
    
    return(res2)
  }
  
  # create a list of the IDs present - make sure it is in numeric order to keep all matrices consistent
  IDs <- sort(unique(c(dataframe$Indiv1, dataframe$Indiv2)))
  # create the empty matrix to fill the values into
  res <- as.data.frame(matrix(0, ncol = length(IDs), nrow = length(IDs))) # create the matrix to save the results into
  
  # Assign the SRI values to the correct positions within the matrix
  for(i in 1 : length(IDs)){
    for(t in 1:length(IDs)){
      values <- subset(dataframe, Indiv1 == IDs[i] & Indiv2 == IDs[t]) # subset out the value if the first id is the column
      values2 <- subset(dataframe, Indiv1 == IDs[t] & Indiv2 == IDs[i]) # subset out the value if the first id is the row
      
      if(nrow(values) == 1){
        res[i,t] <- 1
      } else if (nrow(values2) == 1){
        res[i,t] <- 1
      } else {
        res[i,t] <- 0
      }
    }
  }
  # Add the column and row names to the matrix to make it easier to work with
  rownames(res) <- IDs # add the ID labels to the rows
  colnames(res) <- IDs # add the ID labels to the columns
  
  res <- as.data.frame(res) # convert from a matrix to a data.frame
  res <- rownames_to_column(res, var="ID") # sets the rownames to a column
  ####=======================================================================####
  # Add in the crocodiles who are not present this period and then
  missingIDs <- subset(totalIDs, !(totalIDs %in% IDs))
  
  # For all of the crocodiles not present/HRs could not be determined for this year add a row and column of NA values
  for(h in 1 : length(missingIDs)){
    
    # Add a new row of NA values for this individual
    res[nrow(res)+1, 1] <- missingIDs[h] # can use this to create a new column with the ID of each croc not included and then everything else NA
    
    res$new <- NA # create a new column
    
    colnames(res)[colnames(res)=="new"] <- missingIDs[h] # rename the new column to the individuals name
  }
  
  # Reorder the HR overlap data.frame so that the rows are in ascending order
  
  res2 <- res %>% arrange(ID) # rearrange the rows so that they are in ascending order
  
  res2 <- res2[ , order(names(res2))] # place the columns into ascending order
  
  res2 <- column_to_rownames(res2, var = "ID") # convert ID back into row.names
  
  res2 <- data.matrix(res2, rownames.force = "ID") # convert back into a matrix
  
  return(res2)
}

# Combine all the SRI dataframes into a single one to get a list of all the individuals observed interacting with conspecifics
yearlySRI <- rbindlist(yearly_SRI)

yearlySRI$Year <- yearlySRI$Month # create a Unique month/Year column

years <- unique(yearlySRI$Year) # create a list object of all the years present

ids <- unique(c(yearlySRI$Indiv1, yearlySRI$Indiv2)) # Create a list of all the crocodiles observed associating 

### Convert the SRI dataframe into a matrix for running the LARprep function
registerDoParallel(6) # set the number of cores to use
# Run the function in parallel to speed things up
SRI_matrices <- foreach(p = 1:(length(years)), 
                        .packages = c("data.table","lubridate","plyr","tibble", "dplyr")) %dopar%
  dataframe_Matrix(p, data_list = yearlySRI, totalIDs = ids)

stopImplicitCluster()

####===========================================================================================================####
# The below function takes a list of matrices and then prepares them into the format used by the LAR function from the asnipe package

LARprep <- function(lmig, years, repro_status = c("Total","R","NR"), categories, dynamic = F, sex_split = F, opposite = F){ #years is the list of unique monthly IDs from the crocs data. Will be used to add column to each matrix to assign when they occurred
  
  #############################################################################################
  # lmig: The list of matrices
  # years: An object containing the time sample of interest in the same order as lmig
  # repro_status: Indicate what reproductive status you would like examined. "Total" gives the overall LAR. "R" subsets to just reproductive mature individuals. "NR" subsets to only non-reproductive immature individuals. Combinations of both "R" and "NR" can also be indicated
  # categories: A dataframe containing the identity and reproductive status of each animal
  # dynamic: A TRUE/FALSE indicating if you would like reprodutive status to dynamically change through time. The default is false for a constant reproductive status
  # sex_split: If TRUE the data is further subset so that only conspecifics with the same sex as the focal indivdiual are included in the analysis. Default is F
  # opposite: If TRUE a list of the opposite sex presnt to the focal indivdiual is created and the matrix is then subset to only include these conspecifics. If FALSE all of the opposite conspecifics are removed
  #############################################################################################
  
  # First for each of the matrices present, go through and remove all excess individuals and set overlaps to either 0 or 1
  # then convert each matrix to a data.frame for setting the time point for each individual.
  o <- 1
  for(i in 1 : length(lmig)){
    
    lmig[[i]][is.na(lmig[[i]])] <- 0 # replace all NA values with 0
    
    # Next, set the indviduals who are overlapping with each other as 1 to demonstrate they are apart of a group. If a
    # pair of individuals are 
    #lmig[[i]][lmig[[i]] < 0.01] <- 0 # any value less than 0.05 will be change to 0
    lmig[[i]][lmig[[i]] > 0] <- 1 # any value greater than or equal to 0.05 will be change to 0
    
    # Convert the matrix to a data.frame and then make row.names into a column named TRANSMITTERID
    lmig[[i]] <- as.data.frame(lmig[[i]]) # convert from a matrix to a data.frame
    lmig[[i]] <- lmig[[i]][rowSums(lmig[[i]][, -1])>0, ] # remove any rows that only have 0 present i.e., individuals not found that year
    
    if(nrow(lmig[[i]]) == 0){
      #lmig <- list.remove(lmig, i)
      #lmig[[i]] <- NULL # if after doing above that a month is found to have 0 overlaps, remove from the list
    } else { # else assign a column to include an individuals ID and one for what month this is
      lmig[[i]] <- rownames_to_column(lmig[[i]], var="TRANSMITTERID") # sets the rownames to a column
      lmig[[i]]$Month <- years[o]
      o <- o +1
    }
  }
  
  overlaps <- as.data.frame(Reduce(rbind, lmig)) # Merge each of the matrices in the list into a single dataframe
  overlaps$Month <- as.Date(paste(overlaps$Month, "01", sep = "-"))
  
  #----------------------------------------------------------------------------------------------------------------#
  # Create the time series file so that each individual as their own unique 0 point and time series for when they're present
  myids <- unique(overlaps$TRANSMITTERID)
  
  IDorder <- list() # create an empty list to place each individual into, essentially reordering the detections by individuals
  time_list <- list() # create an empty list to place each individuals time points into
  
  for(t in 1: length(myids)){
    
    if(length(repro_status) == 1){
      if(repro_status == "Total"){
        # IF "Total" is choosen, simply run the function as originally intended to retrieve the data for all
        # potential dyads
        IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
        # Create the time series data for the analysis
        time_between <- 0
        r <- 1
        while(r <= nrow(IDorder[[t]])){
          # create the time series for the individual from 0
          time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
          r <- r + 1
        }
        time_list[[t]] <- time_between
        IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
        IDorder[[t]]$Month <- NULL # remvoe the month column
        
      } else {
        
        if(dynamic == T){ # If it has been selected to do dynamic association rates then first subset out the individuals data and then determine what
          # years were they observed and then run through the function below
          indiv <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out the individual of interest
          year_present <- unique(substr(indiv$Month, 1, 4)) # get a list of the years an indivdual was present
          ID_year <- list() # create a list to store the yearly values into
          time_year <- list()
          for(b in 1:length(year_present)){
            # get a list of all the individuals that are the status of interest for that year
            IDs <- unique((subset(categories, Repro_status == repro_status & Year == year_present[b]))$TRANSMITTERID) 
            
            if(is.element(myids[t], IDs) == T){ # Only if the indivdual is the reproductive status of interest do the below
              ID_year[[b]] <- subset(indiv, substr(indiv$Month, 1, 4) == year_present[b]) # subset out that indivduals data for that year
              # Create the time series data for the analysis
              time_between <- 0
              r <- 1
              while(r <= nrow(ID_year[[b]])){
                # create the time series for the individual from 0
                time_between[[r]] <- as.numeric(difftime(ID_year[[b]]$Month[r],paste0(years[1],"-01")))
                r <- r + 1
              }
              time_year[[b]] <- time_between
              ID_year[[b]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
              ID_year[[b]]$Month <- NULL # remvoe the month column
              
              for(q in 1:length(names(ID_year[[b]]))){ # for all of the columns present
                if(is.element(names(ID_year[[b]][q]), IDs) == F){ # test to see if the individual for that column is in the list
                  ID_year[[b]][,q] <- 0 # if not set their values to 0
                }
              }
              
              # Finally, if sex_split = T then that means we need to then to further refine the data so that we only either have indivduals that are the 
              # same sex or opposite sex of the focal individual
              if(sex_split == T){
                # Determine the individuals sex
                indiv_sex <- subset(categories, TRANSMITTERID == myids[t])$Sex[1]
                
                # Next if opposite == T then we create a list of the opposite sex present and then subset the matrix again so that they are only present
                # Else create a list of all the same sex individuals and then remove all opposite sex ones
                if(opposite == T){
                  Sex_list <- unique(subset(categories, Repro_status == repro_status  & Year == year_present[b] & Sex != indiv_sex)$TRANSMITTERID)
                  for(q in 1:length(names(ID_year[[b]]))){ # for all of the columns present
                    if(is.element(names(ID_year[[b]][q]), Sex_list) == F){ # test to see if the individual for that column is in the list
                      ID_year[[b]][,q] <- 0 # if not set their values to 0
                    }
                  }
                } else {
                  Sex_list <- unique(subset(categories, Repro_status == repro_status  & Year == year_present[b] & Sex == indiv_sex)$TRANSMITTERID)
                  for(q in 1:length(names(ID_year[[b]]))){ # for all of the columns present
                    if(is.element(names(ID_year[[b]][q]), Sex_list) == F){ # test to see if the individual for that column is in the list
                      ID_year[[b]][,q] <- 0 # if not set their values to 0
                    }
                  }
                }
              } 
            }
          }
          if(length(ID_year) > 0){ # if there is at least 1 month present add the values into the total list
            IDorder[[t]] <- rbindlist(ID_year)
            time_list[[t]] <- Reduce(c,time_year)
          } 
        } else {
          # If we are not after the total, susbet out the group of interest
          IDs <- unique((subset(categories, Repro_status == repro_status))$TRANSMITTERID)
          if(is.element(myids[t], IDs) == T){
            #----------------------------------------------------------------------#
            # subset out the focal individuals data
            IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
            # Create the time series data for the analysis
            time_between <- 0
            r <- 1
            while(r <= nrow(IDorder[[t]])){
              # create the time series for the individual from 0
              time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
              r <- r + 1
            }
            time_list[[t]] <- time_between
            IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
            IDorder[[t]]$Month <- NULL # remvoe the month column
            
            for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
              if(is.element(names(IDorder[[t]][q]), IDs) == F){ # test to see if the individual for that column is in the list
                IDorder[[t]][,q] <- 0 # if not set their values to 0
              }
            }
          }
        }
      }
    }  else {
      
      if(dynamic == t){ # If it has been selected to do dynamic association rates then first subset out the individuals data and then determine what
        # years were they observed and then run through the function below
        indiv <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out the individual of interest
        year_present <- unique(substr(indiv$Month, 1, 4)) # get a list of the years an indivdual was present
        ID_year <- list() # create a list to store the yearly values into
        time_year <- list()
        for( b in 1:length(year_present)){
          # Subset out the two lists of interest
          IDs1 <- unique((subset(categories, Repro_status == repro_status[1] & Year == year_present[b]))$TRANSMITTERID)
          IDs2 <- unique((subset(categories, Repro_status == repro_status[2] & Year == year_present[b]))$TRANSMITTERID)
          #----------------------------------------------------------------------#
          if(is.element(myids[t], IDs1) == T | is.element(myids[t], IDs2) == T){
            # subset out the focal individuals data
            ID_year[[b]] <- subset(indiv, substr(indiv$Month, 1, 4) == year_present[b]) # subset out that indivduals data for that year
            # Create the time series data for the analysis
            time_between <- 0
            r <- 1
            while(r <= nrow(ID_year[[b]])){
              # create the time series for the individual from 0
              time_between[[r]] <- as.numeric(difftime(ID_year[[b]]$Month[r],paste0(years[1],"-01")))
              r <- r + 1
            }
            time_year[[b]] <- time_between
            ID_year[[b]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
            ID_year[[b]]$Month <- NULL # remvoe the month column
            #----------------------------------------------------------------#
            # Determine which list the individual of interest belongs to
            # This is important as we want this individuals matrix to reflect how it spatially overlaps with only
            # individuals in the other strategy to it
            if(is.element(myids[b], IDs1) == F){
              # If the focal individual is not in the first list, add it's ID to the first list to create the list
              # of IDs for which to keep the overlaps (1) within the matrix
              IDs1[length(IDs1)+1] <- myids[b] # add the ID of interest to the previous column
              for(q in 1:length(names(ID_year[[b]]))){ # for all of the columns present
                if(is.element(names(ID_year[[b]][q]), IDs1) == F){ # test to see if the individual for that column is in the ID list to keep
                  ID_year[[b]][,q] <- 0 # if not set their values to 0
                }
              }
            } else {
              # If that individual is in that list, add it's ID to the other list and then filter the data
              IDs2[length(IDs2)+1] <- myids[b] # add the ID of interest to the previous column
              for(q in 1:length(names(ID_year[[b]]))){ # for all of the columns present
                if(is.element(names(ID_year[[b]][q]), IDs2) == F){ # test to see if the individual for that column is in the ID list to keep
                  ID_year[[b]][,q] <- 0 # if not set their values to 0
                }
              }
            }
          }
        }
        if(length(ID_year) > 0){ # if there is at least 1 month present add the values into the total list
          IDorder[[t]] <- rbindlist(ID_year)
          time_list[[t]] <- Reduce(c,time_year)
        } 
        
      } else {
        # Subset out the two lists of interest
        IDs1 <- unique((subset(categories, Repro_status == repro_status[1]))$TRANSMITTERID)
        IDs2 <- unique((subset(categories, Repro_status == repro_status[2]))$TRANSMITTERID)
        #----------------------------------------------------------------------#
        if(is.element(myids[t], IDs1) == T | is.element(myids[t], IDs2) == T){
          # subset out the focal individuals data
          IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
          # Create the time series data for the analysis
          time_between <- 0
          r <- 1
          while(r <= nrow(IDorder[[t]])){
            # create the time series for the individual from 0
            time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
            r <- r + 1
          }
          time_list[[t]] <- time_between
          IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
          IDorder[[t]]$Month <- NULL # remvoe the month column
          #----------------------------------------------------------------#
          # Determine which list the individual of interest belongs to
          # This is important as we want this individuals matrix to reflect how it spatially overlaps with only
          # individuals in the other strategy to it
          if(is.element(myids[t], IDs1) == F){
            # If the focal individual is not in the first list, add it's ID to the first list to create the list
            # of IDs for which to keep the overlaps (1) within the matrix
            IDs1[length(IDs1)+1] <- myids[t] # add the ID of interest to the previous column
            for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
              if(is.element(names(IDorder[[t]][q]), IDs1) == F){ # test to see if the individual for that column is in the ID list to keep
                IDorder[[t]][,q] <- 0 # if not set their values to 0
              }
            }
          } else {
            # If that individual is in that list, add it's ID to the other list and then filter the data
            IDs2[length(IDs2)+1] <- myids[t] # add the ID of interest to the previous column
            for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
              if(is.element(names(IDorder[[t]][q]), IDs2) == F){ # test to see if the individual for that column is in the ID list to keep
                IDorder[[t]][,q] <- 0 # if not set their values to 0
              }
            }
          }
        }
      }
    }
    
    
    #time_list[[t]] <- time_between
    
    #time_list[[t]] <- seq(0, by = 365, length.out = nrow(IDorder[[t]])) # create the time sequence that they are present
  }
  
  
  IDorder <- data.frame(Reduce(rbind, IDorder)) # convert the list to a data.table
  IDlist <- unique(IDorder$TRANSMITTERID) # create a list of the ID's in the order they're in the data frame
  #IDorder$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
  #IDorder$Month <- NULL
  
  IDorder <- as.matrix(IDorder) # then convert ID order into a matrix to be used in the LAR function
  
  # Combine the list of vectors for the time data
  time_list <- Reduce(c,time_list)
  
  # combine the times vector and the overlap matrix into a list to be output from the function
  output <- list(IDorder, time_list, IDlist)
  
  return(output) 
  # Add a new string of vectors to the times list to indicate this time point
  #timepoint <- nrow(lmig[[i]])
  #time_list[[i]] <- rep(time_vec[i], timepoint)
}


# Import the dataframe with the dynamic reproductive status of individuals through time created in 2_Crocodile_social_organisation
dynamic_repro_status <- fread("Data/Metrics/Dynamic_repro_status.csv")

###################################################################################################################################
#### Calculate the LAR for each possible maturity status combination of interest

# Mature same sex dyads (i.e., MM, FF)
R_associations <- LARprep(lmig = SRI_matrices, years = years, repro_status = c("R"), categories = dynamic_repro_status, dynamic = T, sex_split = T, opposite = F)
R_lagged_rate <- LAR(R_associations[[1]],R_associations[[2]],30, identities = ids) %>%
  as.data.frame()
# Mature opposite sex dyads (MF)
R_associations_MF <- LARprep(lmig = SRI_matrices, years = years, repro_status = c("R"), categories = dynamic_repro_status, dynamic = T, sex_split = T, opposite = T)
R_lagged_rate_MF <- LAR(R_associations_MF[[1]],R_associations_MF[[2]],30, identities = ids) %>%
  as.data.frame()
# Immature dyads
NR_associations <- LARprep(lmig = SRI_matrices, years = years, repro_status = c("NR"), categories = dynamic_repro_status, dynamic = T)
NR_lagged_rate <- LAR(NR_associations[[1]],NR_associations[[2]],30, identities = ids) %>%
  as.data.frame()
# For some reason there are NaN's present in this dataframe. As they are surrounded by 0 it appears to be that nothing is occurring, or missing data was present. 
# Replace the NaN's with 0 to visualise things better
NR_lagged_rate$V2 <- ifelse(is.na(NR_lagged_rate$V2) == T, 0, NR_lagged_rate$V2)
# Mature-Immature dyads
RNR_associations <- LARprep(lmig = SRI_matrices, years = years, repro_status = c("R", "NR"), categories = dynamic_repro_status, dynamic = T)
RNR_lagged_rate <- LAR(RNR_associations[[1]],RNR_associations[[2]],30, identities = ids) %>% 
  as.data.frame()


######################################################################################################################################
#### Duplicate the above steps for each of the random simulations

# This fuction duplicates the above steps for calculating the dynamic LAR for each random simulation. It automatically reads in the simulated datasets and home range overlaps
simulatedSRIs <- function(i, years1, class = "Total", cat = croc_metric, dyn = F, sex = F, opp = F){
  
  ################################################################################################################
  # i: Vector indicating which set of random simulation you would like selected
  # years1: An object containing all of the potential unique months for the duration of the study in the same order as the observed data
  # class: The reproductive status of interest. Values the same as used in the LARprep function "Total", "R", "NR"
  # dyn: A TRUE/FALSE indicating if you would like reprodutive status to dynamically change through time. The default is false for a constant reproductive status
  # sex: If TRUE the data is further subset so that only conspecifics with the same sex as the focal indivdiual are included in the analysis. Default is F
  # opp: If TRUE a list of the opposite sex presnt to the focal indivdiual is created and the matrix is then subset to only include these conspecifics. If FALSE all of the opposite conspecifics are removed
  ################################################################################################################
  
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
      sim_overlap[upper.tri(sim_overlap, diag = T)] <- NA 
      sim_overlap <- na.omit(as.data.table(as.table(sim_overlap, na="", row.names=T, col.names=T))) # converts from matrix to table
      names(sim_overlap) <- c("ID1", "ID2","HRoverlap") # Rename the columns
      sim_overlap$Month <- years1[t]
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
  # Create a column within the dataset, which places the data into the correct years for the dynamic networks
  sim_co_occur$Month <- substr(sim_co_occur$DatetimeID1,1,7)
  
  # Now determine the monthly SRI of the simulated data
  monthly_SRI <- list()
  monthSRI <- function(q){ # only run this function for the simulated months which had enough data to generate home ranges
    year_sri <- sim_co_occur[Month == year[q]] # subset out all of the co-detections for that year
    year_HR <- HRdataframes[Month == year[q]] # extract the home range overlaps for the year of interest
    individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
    year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
    
    output <- as.data.table(prepAssocIndex(year_sri))
    output <- output[,c(1,2,7)] # reduce down to just the first, second and seventh column for plotting
    output$Month <- year[q] 
    # Save it into a list object
    return(output)
    
  }
  monthly_SRI <- lapply(1:length(year), monthSRI)
  monthly_SRI <- rbindlist(monthly_SRI) # combine the list into a dataframe
  monthly_SRI$Year <- monthly_SRI$Month
  monthly_SRI$Month <- substr(monthly_SRI$Month, 6,7)
  #monthly_SRI <- monthly_SRI[Month != "02" & Month != "03" & Month != "04"]
  monthly_SRI$Month <- NULL # remove the month column
  years <- unique(monthly_SRI$Year) # create a list object of all the years present
  
  ids <- unique(c(monthly_SRI$Indiv1, monthly_SRI$Indiv2))
  
  SRI_matrices <- list()
  for(f in 1:length(years)){
    SRI_matrices[[f]] <- dataframe_Matrix(f, data_list = monthly_SRI, totalIDs = ids)
  }
  
  full_associations <- LARprep(lmig = SRI_matrices, years = years, repro_status = class, categories = cat, dynamic = dyn, sex_split = sex, opposite = opp)
  groupsdata <- full_associations[[1]] # extract the groups
  time_points <- full_associations[[2]] # extract the time series of the data
  lagged_rate <- LAR(groupsdata,time_points,30, identities = ids)
  
  lagged_rate <- as.data.frame(lagged_rate)
  lagged_rate$Simulation <- paste0("Simulation-",i)
  
  return(lagged_rate)
}


###################################################################################################################################################
#### Calculate the simulated LAR values to act as the Null model and then plot the results with the observed data #####

##-------------------------------------------------------------------------------------##
# Mature same sex dyads (MM, FF)
registerDoParallel(6) # set the number of cores to use

sim_LAR_R <- foreach(p = 1:50, 
                     .packages = c("data.table","lubridate","plyr", "tibble", "dplyr","asnipe")) %dopar%
  simulatedSRIs(p, years1 = years1, class = c("R"), cat = dynamic_repro_status, dyn = T, sex = T, opp = F)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

sim_LAR_R <- rbindlist(sim_LAR_R) # combine all of the simulations into a single dataframe
sim_LAR_R <- sim_LAR_R[,.(MeanLAR = mean(V2, na.rm = T), seLAR = sd(V2, na.rm = T)), by = V1] # determine the average and SD of the simulated LARs to form the null
sim_LAR_R <- as.data.frame(sim_LAR_R)
sim_LAR_R$lower <- ifelse(sim_LAR_R$MeanLAR - sim_LAR_R$seLAR < 0, 0, sim_LAR_R$MeanLAR - sim_LAR_R$seLAR)
sim_LAR_R$upper <- sim_LAR_R$MeanLAR + sim_LAR_R$seLAR

Repro_repro <- ggplot()+
  geom_line(data = sim_LAR_R, aes(x = exp(sim_LAR_R[,1])/12, y = sim_LAR_R[,2]), colour = "red")+
  geom_ribbon(data = sim_LAR_R, aes(x = exp(sim_LAR_R[,1])/12, ymin = sim_LAR_R[,4], ymax = sim_LAR_R[,5]), fill = "red", alpha = 0.3)+
  geom_line(data = R_lagged_rate, aes(x = exp(R_lagged_rate[,1])/12, y = R_lagged_rate[,2]))+
  ylim(0, 1) + ylab("Lagged assocation rate") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  annotate("text",x = 4.3, y = 1, 
           label = "Mature - Mature",
           size = 4)+
  labs(tag = "(a)")
Repro_repro

##-------------------------------------------------------------------------------------##
# Mature opposite sex dyads (MF)
registerDoParallel(6) # set the number of cores to use

sim_LAR_R_MF <- foreach(p = 1:50, 
                        .packages = c("data.table","lubridate","plyr", "tibble", "dplyr","asnipe")) %dopar%
  simulatedSRIs(p, years1 = years1, class = c("R"), cat = dynamic_repro_status, dyn = T, sex = T, opp = T)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

sim_LAR_R_MF <- rbindlist(sim_LAR_R_MF)
sim_LAR_R_MF <- sim_LAR_R_MF[,.(MeanLAR = mean(V2, na.rm = T), seLAR = sd(V2, na.rm = T)), by = V1]
sim_LAR_R_MF <- as.data.frame(sim_LAR_R_MF)
sim_LAR_R_MF$lower <- ifelse(sim_LAR_R_MF$MeanLAR - sim_LAR_R_MF$seLAR < 0, 0, sim_LAR_R_MF$MeanLAR - sim_LAR_R_MF$seLAR)
sim_LAR_R_MF$upper <- sim_LAR_R_MF$MeanLAR + sim_LAR_R_MF$seLAR

Repro_repro_MF <- ggplot()+
  geom_line(data = sim_LAR_R_MF, aes(x = exp(sim_LAR_R_MF[,1])/12, y = sim_LAR_R_MF[,2]), colour = "red")+
  geom_ribbon(data = sim_LAR_R_MF, aes(x = exp(sim_LAR_R_MF[,1])/12, ymin = sim_LAR_R_MF[,4], ymax = sim_LAR_R_MF[,5]), fill = "red", alpha = 0.3)+
  geom_line(data = R_lagged_rate_MF, aes(x = exp(R_lagged_rate_MF[,1])/12, y = R_lagged_rate_MF[,2]))+
  ylim(0, 1) + ylab("Lagged assocation rate") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  annotate("text",x = 4.3, y = 1, 
           label = "Mature - Mature MF",
           size = 4)+
  labs(tag = "(b)")
Repro_repro_MF

##-------------------------------------------------------------------------------------##
# Mature-Immature dyads
registerDoParallel(6) # set the number of cores to use

sim_LAR_NR <- foreach(p = 1:50, 
                      .packages = c("data.table","lubridate","plyr", "tibble", "dplyr","asnipe")) %dopar%
  simulatedSRIs(p, years1 = years1, class = c("R","NR"), cat = dynamic_repro_status, dyn = T)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

sim_LAR_NR <- rbindlist(sim_LAR_NR)
sim_LAR_NR <- sim_LAR_NR[,.(MeanLAR = mean(V2, na.rm = T), seLAR = sd(V2, na.rm = T)), by = V1]
sim_LAR_NR <- as.data.frame(sim_LAR_NR)
sim_LAR_NR$lower <- ifelse(sim_LAR_NR$MeanLAR - sim_LAR_NR$seLAR < 0, 0, sim_LAR_NR$MeanLAR - sim_LAR_NR$seLAR)
sim_LAR_NR$upper <- sim_LAR_NR$MeanLAR + sim_LAR_NR$seLAR

Repro_nonrepro <- ggplot()+
  geom_line(data = sim_LAR_NR, aes(x = exp(sim_LAR_NR[,1])/12, y = sim_LAR_NR[,2]), colour = "red")+
  geom_ribbon(data = sim_LAR_NR, aes(x = exp(sim_LAR_NR[,1])/12, ymin = sim_LAR_NR[,4], ymax = sim_LAR_NR[,5]), fill = "red", alpha = 0.3)+
  geom_line(data = RNR_lagged_rate, aes(x = exp(RNR_lagged_rate[,1])/12, y = RNR_lagged_rate[,2]))+
  ylim(0, 1) + ylab("Lagged assocation rate") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  annotate("text",x = 4.3, y = 1, 
           label = "Mature - Immature",
           size = 4)+
  labs(tag = "(c)")
Repro_nonrepro

##-------------------------------------------------------------------------------------##
#Immature dyads
registerDoParallel(6) # set the number of cores to use

sim_LAR_NRNR <- foreach(p = 1:50, 
                        .packages = c("data.table","lubridate","plyr", "tibble", "dplyr","asnipe")) %dopar%
  simulatedSRIs(p, years1 = years1, class = c("NR"), cat = dynamic_repro_status, dyn = T)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

sim_LAR_NRNR <- rbindlist(sim_LAR_NRNR)
sim_LAR_NRNR <- sim_LAR_NRNR[,.(MeanLAR = mean(V2, na.rm = T), seLAR = sd(V2, na.rm = T)), by = V1]
sim_LAR_NRNR <- as.data.frame(sim_LAR_NRNR)
sim_LAR_NRNR$lower <- ifelse(sim_LAR_NRNR$MeanLAR - sim_LAR_NRNR$seLAR < 0, 0, sim_LAR_NRNR$MeanLAR - sim_LAR_NRNR$seLAR)
sim_LAR_NRNR$upper <- sim_LAR_NRNR$MeanLAR + sim_LAR_NRNR$seLAR

nonepro_nonrepro <- ggplot()+
  geom_line(data = sim_LAR_NRNR, aes(x = exp(sim_LAR_NRNR[,1])/12, y = sim_LAR_NRNR[,2]), colour = "red")+
  geom_ribbon(data = sim_LAR_NRNR, aes(x = exp(sim_LAR_NRNR[,1])/12, ymin = sim_LAR_NRNR[,4], ymax = sim_LAR_NRNR[,5]), fill = "red", alpha = 0.3)+
  geom_line(data = NR_lagged_rate, aes(x = exp(NR_lagged_rate[,1])/12, y = NR_lagged_rate[,2]))+
  ylim(0, 1) + ylab("Lagged assocation rate") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  annotate("text",x = 4.3, y = 1, 
           label = "Immature - Immature",
           size = 4)+
  labs(tag = "(d)")
nonepro_nonrepro

