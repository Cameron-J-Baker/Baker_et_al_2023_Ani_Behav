# This function calculates the simple ratio index (SRI) from the output of the coAcoustic function
prepAssocIndex <- function(allco){
  
  ##################################################
  # allco: The output from the coAcoustic function
  ##################################################
  
  # Generate the empty vectors to place the data into
  Ya <- NULL
  Yb <- NULL
  Yab <- NULL
  icount <- 1
  
  # Aggregate down to only unique associations between individuals with Weight being the number of interactions between a pair
  edge <- allco[,.(Weight = sum(Weight)), by = c("Indiv1","Indiv2")]
  edge <- as.data.frame(edge)
  # Calculate Ya, Yb and Yab
  while(icount <= nrow(edge)){
    
    # For individual A
    test <- subset(edge, Indiv1 == edge[icount,1] | Indiv2 == edge[icount,1])
    
    # If an individual has only been seen with one other 
    if(nrow(test) == 1){
      
      Ya[[icount]] <- 0 
      
    } else {
      
      test <- subset(test, Indiv2 != edge[icount,2]) # remove individual B to calculate the number of times A was interacting with others
      
      Ya[[icount]] <- sum(test$Weight) # sum together every other time individual A was observed interacting
      
    }
    #---------------------------------------------------------------------------------------------------------------------#   
    # For individual B
    #---------------------------------------------------------------------------------------------------------------------#
    test2 <- subset(edge, Indiv1 == edge[icount,2] | Indiv2 == edge[icount,2])
    
    # If an individual has only been seen with one other 
    if(nrow(test2) == 1){
      
      Yb[[icount]] <- 0 
      
    } else {
      
      test2 <- subset(test2, Indiv2 != edge[icount,1]) # remove individual A to calculate the number of times B was interacting with others
      
      Yb[[icount]] <- sum(test2$Weight)# sum together every other time individual A was observed interacting
      
    }
    #---------------------------------------------------------------------------------------------------------------------#
    # Calculate Yab
    #---------------------------------------------------------------------------------------------------------------------#
    # step 1 subset out all of the detections of indiviual A from the timed dataset
    indivA <- subset(allco, Indiv1 == edge[icount,1] | Indiv2 == edge[icount,1])
    
    # from this remove any interactions recorded between A and B
    indivA <- subset(indivA,  Indiv1 != edge[icount,2])
    indivA <- subset(indivA,  Indiv2 != edge[icount,2])
    
    # Do this same thing for individual B
    indivB <- subset(allco, Indiv1 == edge[icount,2] | Indiv2 == edge[icount,2])
    
    # from this remove any interactions recorded between A and B
    indivB <- subset(indivB, Indiv1 != edge[icount,1] )
    indivB <- subset(indivB, Indiv2 != edge[icount,1] )
    
    # Search for any matches where individual A and B were both observed but not interacting
    result <- indivA$Datetime %in% indivB$Datetime
    
    Yab[[icount]] <- sum(result)
    #---------------------------------------------------------------------------------------------------------------------#
    
    icount <- icount + 1
    
  }
  
  # Add Ya, Yb and Yab to edge 
  edge$Ya <- as.numeric(Ya)
  edge$Yb <- as.numeric(Yb)
  edge$Yab <- as.numeric(Yab)
  
  # Calculate the simple association ratio index for each pair
  edge$SRI <- edge$Weight/(edge$Ya + edge$Yb + edge$Yab + edge$Weight)
  
  return(edge)
}
