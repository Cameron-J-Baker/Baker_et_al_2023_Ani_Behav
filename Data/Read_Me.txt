Baker_et_al_data.csv - This file contains the raw acoustic detections of the 166 estuarine crocodiles (Crocodylus porosus) that have been tagged and detected in 
the Wenlock River from August 2010 - August 2020. As some crocodiles were fitted with multiple acoustic transmitters throughout the duration of this study, this 
data has been processed so that each individual has been assigned a single identification code. The data has been processed to remove all tag detections from 
outside of the period they were implanted within an inidividual and has been reformated for analysis using the 'ReadInputData' function for the VTrack R package. 
The unprocessedraw data supporting this study can be sourced from the Acoustic Animal Tracking Database (https://animaltracking.aodn.org.au) of the Integrated 
Marine Observing System (IMOS, www.imos.org.au) - IMOS is a national collaborative research infrastructure supported by Australian Government. The database is 
a centralised acoustic telemetry data repository maintained by the IMOS Animal Tracking Facility and the Australian Ocean Data Network (AODN, https://portal.aodn.org.au/). 

Baker_et_al_acoustic_hydrophone_location.csv - The location (latitude and longitiude) of each of the acoustic hydrophone stations that have been used as part of the 
Wenlock and Ducie river acoustic telemetry array. Included with the dataset is also the distance (km) that each hydrophone station is from the mouth of Port Musgrave in 
both eulidean (sDISTANCE) and distance along the river system via least cost paths (rDISTANCE). 

Baker_et_al_hydrophone_distance_matrix.csv - A matrix of the distances along the river system via least cost paths between each of the acoustic hydrophone stations throughout
the Wenlock and Ducie River systems in Cape York Queensland.

Baker_et_al_crocodile_metrics.csv - A dataframe detailing the unique identitiy and body size metrics of each of the 166 estuarine crocodiles that have been tagged since 
2008.

Baker_et_al_init_croc_size.csv - A dataframe indicating the year and initial body size of each estuarine crocodiles when they where initally captured and tagged as part 
of this long-term study. 