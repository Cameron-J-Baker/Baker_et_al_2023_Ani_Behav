# Author: Cameron J Baker
# 2022-12-12
# Email: cameron.baker@uqconnect.edu.au/cameron.baker@cdu.edu.au
# Code for the analyses contained within the manuscript "Long-term tracking reveals a dynamic crocodylian social system" published in
#  the journal Animal Behaviour
# Aim: This script examines the social structure of crocodiles following the "Crocodile social structure" method section

library(data.table)
library(igraph)
library(tidygraph)
library(ggraph)
library(plyr)
library(lubridate)
library(tnet)
library(nlme)
library(rcompanion)
library(brms)
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(mcmcplots)
library(lme4)
library(performance)
library(doParallel)
library(foreach)
library(ggpubr)

######################################################################################################################
#### Import and prepare the data for analysis ####

# Import the output from the coAcoustic function
croc_codetect <- fread("Data/crocodile_edges_2020_4.csv")

# Create a column which denotes each unique study month (i.e., 2020-08) to allow monthly networks to be created
croc_codetect$Unique_month <- substr(croc_codetect$DatetimeID1, 1,7)

# Create a vector with all of the unique months in order
Unique_months <- sort(unique(croc_codetect$Unique_month))

# Import the crocodile metrics dataframe to act as the nodes for the social networks
croc_metric <- fread("Data/Metrics/Baker_et_al_crocodile_metrics.csv") 

##########################################################################################################################
# Home ranges were constructed following the methods described in Baker et. al (2022) Crocodile social environments dictated by male philopatry. Behavioral Ecology, 33(1), 156-166. doi: https://doi.org/10.1093/beheco/arab120
# The required R code and shape files can be found at:  https://github.com/Cameron-J-Baker/Supporting-code-Baker_et-al 
# Import the observed home ranges overlap matrices produced following the code above

files <- list.files(path = "Overlap_matrices/Monthly_Wenlock/", full.names = T) # create a vector list with the names of all files

lmig <- lapply(files, fread) # Import all files using lapply

# This for loop first converts the imported datafiles into a matrix of home range overlap and then creates a dataframe with the ID of individuals and 
# their observed spatial overlap
for(i in 1:length(lmig)){
  lmig[[i]] <- as.data.frame(lmig[[i]]) # convert to data.frame
  colnames(lmig[[i]]) = lmig[[i]][1, ] # convert the first row back to column names
  lmig[[i]] = lmig[[i]][-1, ] # remove the first row with the names  
  rownames(lmig[[i]]) <- c() # remove the current row names
  lmig[[i]] <- tibble::column_to_rownames(lmig[[i]], var = "NA") # convert iD back into row.names
  
  lmig[[i]] <- as.matrix(lmig[[i]])
  # Convert the matrix to a dataframe for analysis
  lmig[[i]][upper.tri(lmig[[i]], diag = T)] <- NA # To avoid duplicating overlaps, assign the top triangle of the matrix NA
  # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
  lmig[[i]] <- na.omit(as.data.table(as.table(lmig[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
  names(lmig[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns
  
}

yearlyHR <- plyr::compact(lmig)

######################################################################################################################
#### Determine the monthly social networks ####

#=======================================================================================================================#
# import the functions that will be required
#=======================================================================================================================#
source("R_Functions/prepAssocIndex.R")

# Next determine the social network for each of the unique months in the study with associations present
yearly_SRI <- list() # Create an empty list to store the month SRI values and networks
no_indivs <- list() # Create an empty list to store the total observed crocodiles in each network
for(i in 1:length(years)){
  year_sri <- raw_dyads[Year == years[i]] # subset out all of the co-detections for that year
  year_HR <- yearlyHR[[i]] # extract the home range overlaps for the year of interest
  individuals <- unique(c(year_HR$ID1, year_HR$ID2)) # create a list of the individuals with a HR estimate
  year_sri <- year_sri[Indiv1 %in% individuals & Indiv2 %in% individuals] # subset out the individuals without enough data
  
  output <- prepAssocIndex(year_sri)
  
  output <- output[,c(1,2,7)] # reduce down to just the first, second and seventh column for plotting
  colnames(output)[colnames(output)=="Indiv1"] <- "from"
  colnames(output)[colnames(output)=="Indiv2"] <- "to"
  colnames(output)[colnames(output)=="SRI"] <- "weight"
  
  # Run the prepAssocIndex function and then save it into a list
  yearly_SRI[[i]] <- output
  no_indivs[[i]] <- data.table(Year = years[i],
                               No_indivs_network = length(unique(c(output$from, output$to))),
                               No_indivs_total = length(individuals))
}
no_indivs <- rbindlist(no_indivs)


##################################################################################################################################################
#### Determine the network level measures of the social structure of crocodile populations throughout each month #####
##################################################################################################################################################
# First step, create a function which will extract for each month calculate the network transitivity, number of clusters and modularity

network_metrics_Extract <- function(i){
  print(i)
  year_HR <- yearlyHR[[i]]
  
  ids <- unique(c(year_HR$ID1, year_HR$ID2))
  
  year_connect <- yearly_SRI[[i]]
  
  #year_connect <- year_connect[weight > 0]
  
  nodes <- croc_metric[id %in% ids]
  interactions <- graph_from_data_frame(d = year_connect, vertices=nodes, directed=F)
  
  community <- cluster_leading_eigen(interactions, options=list(maxiter=1000000))
  
  output <- data.table(Year = years[i],
                       Transitivity = transitivity(interactions, type = "global", isolates = "zero"),
                       Number_clusters = length(community),
                       Modularity = modularity(community))
  return(output)
}

# Run the above function through each of the observed months SRI data
network_metrics <- lapply(1: length(years), network_metrics_Extract)
# Combine the outputs of the above function into a single dataframe
network_metrics <- rbindlist(network_metrics)
# Combine with the number of individuals output from above for each month
network_metrics <- plyr::join(network_metrics, no_indivs, by = "Year")

network_metrics <- na.omit(network_metrics) # Exclude any potential NAs that have formed in the dataframe
# Determine the mean and SD transtivity and modularity of the population through time
mean(network_metrics$Transitivity) 
sd(network_metrics$Transitivity)
mean(network_metrics$Modularity)
sd(network_metrics$Modularity)
# Create separate month and year columns from the Unique month column for analyses below
network_metrics$Month <- as.factor(substr(network_metrics$Year, 6,nchar(network_metrics$Year)))
network_metrics$Year <- as.factor(substr(network_metrics$Year, 1,4))

############################################################################################################################################
#### Does time of year influence the transitivity of the population? ####

# Create the model
mod1_Trans <- lmer((Transitivity) ~ Month + (1|Year), data = network_metrics)

check_model(mod1_Trans) # Check the diagnostic plots of the model
car::Anova(mod1_Trans) # There is a significant effect of time of year on the transitvity of the population
summary(mod1_Trans)
# Generate and save the predictions from the model for creating Figure 6
month_trans <- as.data.table(effects::predictorEffect("Month", mod1_Trans))
month_trans$lower <- ifelse(month_trans$lower < 0, 0, month_trans$lower) # As transitivity can not be below 0, set the lower bounds to 0

# Export a copy of the model predictions for plotting later on
write.csv(month_trans, "Data/Model_predictions/Transitivity.csv")

############################################################################################################################################
#### Does time of year influence the modularity of the population? ####

# Create the model
mod1_modul <- lmer(Modularity ~ Month + (1|Year), data = network_metrics)

check_model(mod1_modul) # Check the diagnostic plots of the model
car::Anova(mod1_modul)# There is a significant effect of time of year on the transitvity of the population
summary(mod1_modul)
# Generate and save the predictions from the model for creating Figure 6
month_modul <- as.data.table(effects::predictorEffect("Month", mod1_modul))
month_modul$lower <- ifelse(month_modul$lower < 0, 0, month_modul$lower) # As modularity can not be below 0, set the lower bounds to 0

# Export a copy of the model predictions for plotting later on
write.csv(month_modul, "Data/Model_predictions/Modularity.csv")