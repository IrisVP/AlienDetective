
##################################################################################
# This script is an alternative script for calculating sea distances and fly distances
##################################################################################

############################################################################################
# LOAD PACKAGES
############################################################################################
library("gdistance")
library("dplyr")
require("geosphere")
require("rgbif")
library("ggplot2")
library("tidyr")
library("worrms")
library("sf")
library("sp")
library("raster")
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library("FRK")

# instructions to download the rnaturalearthhires
# Install the devtools package if it's not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the rnaturalearthhires package from GitHub
devtools::install_github("ropensci/rnaturalearthhires")

# Load the necessary library
library(rnaturalearthhires)

############################################################################################
# LOAD DATA
############################################################################################
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Read a species list
df <- read.csv("Output_preparation/Species_Location_defence.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")


############################################################################################
# RESHAPE DF WITH LOCATIONS AND SPECIES FROM WIDE TO LONG FORMAT
############################################################################################

long <- pivot_longer(df, !Specieslist)
# from tidyr package, reshape df from wide to long format
long <- long[long$value > 0, ]

#####################################################################################
# FIRST CHECK IF FILE WITH OCCURRENCE DATA IN OCCURRENCEDATA DIRECTORY EXISTS
# IF NOT: MAKE ONE AND GET DATA FROM GBIF
# limit in occ_data() function is changeable for personal preference
####################################################################################

# Function to ensure all required columns are present (used in coming for loop)
ensure_columns <- function(df, required_columns) {
  missing_columns <- setdiff(required_columns, colnames(df))
  for (col in missing_columns) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  return(df)
}
fetch_data_in_batches <- function(species_name, basisOfRecord, batch_size = 10000) {
  #start <- 0
  combined_data <- data.frame()
    res <- occ_data(scientificName = species_name, 
                    hasCoordinate = TRUE, 
                    limit = batch_size,
                    basisOfRecord = basisOfRecord,
                    continent = "europe")
  combined_data <- rbind(combined_data, res$data)
  Sys.sleep(2)  # Adding a small delay to be polite to the server
  
  return(combined_data)
}

required_columns <- c('decimalLongitude', 'decimalLatitude', 'year', 'month', 'country')

for (species_name in unique(long$Specieslist)){
  print(paste0("species_name: ", species_name))
  # if file exists:        put content into variable res
  # if file doesn't exist: make one, get data, and put into variable res
  occurrence_coord <- paste0("OccurrenceData/", species_name, ".csv")
  if (file.exists(occurrence_coord) == TRUE) {
    res <- read.csv(occurrence_coord, header = TRUE)
  } else {
    data_list <- list(
      fetch_data_in_batches(species_name, "OBSERVATION"),
      fetch_data_in_batches(species_name, "MACHINE_OBSERVATION"),
      fetch_data_in_batches(species_name, "HUMAN_OBSERVATION"),
      fetch_data_in_batches(species_name, "MATERIAL_SAMPLE"),
      fetch_data_in_batches(species_name, "LIVING_SPECIMEN"),
      fetch_data_in_batches(species_name, "OCCURRENCE"))

    # Initialize an empty list to store processed data frames
    processed_data <- list()
  
    # Loop over each data frame in the list, ensure columns, and select required columns
    for (i in seq_along(data_list)) {
      temp_df <- data_list[[i]]
    
      if (!is.null(temp_df) && nrow(temp_df) > 0) { # check if temp_df is not NULL and not empty
        # Ensure the required columns are present
        temp_df <- ensure_columns(temp_df, required_columns)
        temp_df <- temp_df[, required_columns]
        processed_data[[i]] <- temp_df
        
      } else {
        print(paste0("Data frame ", i, " is empty."))
      }
    }
  }
    
  if (length(processed_data) > 0) {
    res_total <- do.call(rbind, processed_data)
    print(paste0("Total NA values in res_total: ", sum(is.na(res_total))))
  } else {
    print(paste0("No data to combine for species: ", species_name))
  }
  #rename the column names
  colnames(res_total) <- c('Longitude', 'Latitude', 'year', 'month', 'country')
  # Remove occurrences where longitude or latitude is NA
  res_total <- res_total[!is.na(res_total$Latitude) & !is.na(res_total$Longitude),]
  # check if there's no information for a species
  if (nrow(res_total) == 0) {
    error_message <- paste0("No information found on GBIF for ", species_name)
    
    # check if directory with error messages exists, if it doesn't: make one
    if (!dir.exists("Output_calculations/errors")){
      dir.create("Output_calculations/errors", recursive = TRUE)
      error_file_name <- paste0("Output_calculations/errors/error_", species_name, ".csv")
      # error files written to Output_calculations/errors/
      write.csv(error_message, file = error_file_name)
      return(FALSE)
    
      # write error file to the directory
    } else {
      error_file_name <- paste0("Output_calculations/errors/error_", species_name, ".csv")
      # error files written to Output_calculations/errors/
      write.csv(error_message, file = error_file_name)
      return(FALSE)
    }
  }
  write.csv(res_total, file = occurrence_coord)
}


############################################################################################
# REVISION: CALCULATE DISTANCES
############################################################################################

# Iterate over species_name and location_name
Calculation_seadistance <- function(species_name, species_location){
  # Initialize an empty list to store error messages
  error_messages <- list()
  
  # Define a helper function to add error messages
  add_error_message <- function(message) {
    error_messages <<- c(error_messages, message)
  }
  # Print species name and location
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  
  # Try to read occurrence data
  OccurrenceData <- tryCatch({
    occurrence_coord <- paste0("OccurrenceData/", species_name, ".csv")
    read.csv(occurrence_coord, header = TRUE)
  }, error = function(e) {
    add_error_message(paste("Error reading occurrence data for", species_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(OccurrenceData)) return(list(result = NA, error_messages = error_messages))
  
  # Try to get coordinates for ARMS location
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(species_location, Coordinates$Observatory.ID)
  if (length(location_row_index) == 0) {
    add_error_message(paste("Location not found in Coordinates for", species_location))
    return(list(result = NA, error_messages = error_messages))
  }
  # save longitude and latitude for the row that you selected with grep
  longitude <- Coordinates[location_row_index, "Longitude"]
  latitude <- Coordinates[location_row_index, "Latitude"]
  # make a dataframe out of the longitude and latitude called samplelocation
  samplelocation <- data.frame(Latitude = latitude, Longitude = longitude)
  
  # check if OccurrenceData has coordinates
  if (nrow(OccurrenceData) < 1) {
    add_error_message("OccurrenceData has no coordinates")
    return(list(result = NA, error_messages = error_messages))
  }
  
  ##########################################################################
  # DISTANCES CALCULATION
  ##########################################################################

  # Load world data and prepare the raster
  world <- ne_countries(scale = "medium", returnclass = "sf")  # Load medium or large scale natural earth countries as an sf (simple features) object
  r <- raster(extent(-180, 180, -90, 90), res = 0.1)           # Create a raster object with a global extent and resolution of 0.1 degrees
  r <- rasterize(world, r, field = 1, fun = max, na.rm = TRUE) # Rasterize the 'world' sf object, assigning a value of 1 to cells with country presence
  costs <- reclassify(r, cbind(1, Inf))                        # Reclassify the raster: convert all values of 1 to Inf (infinity)
  costs[is.na(costs)] <- 1    # Replace NA values in the 'costs' raster with 1
  
  
  # Initialize lists to store distances
  sea_distances <- c()
  flying_distances <- c()
  
  # for loop to iterate over OccurrenceData
  for (row in 1:nrow(OccurrenceData)) {
    print(paste0("Calculating latitude: ", OccurrenceData[row, 3], " and longitude: ", OccurrenceData[row, 2]))
    print(paste0("for ", species_location, " latitude, longitude: ", samplelocation[,1], " ", samplelocation[,2]))
    
    #################
    ## SEA DISTANCE##
    #################
    
    transition_matrix <- "transitMatrix.rds"
    if (!file.exists(transition_matrix)) {
      # Create a transition object for adjacent cells
      transitMatrix <- transition(costs, transitionFunction = function(x) 1/mean(x), directions = 16)
      # Set infinite costs to NA to prevent travel through these cells
      transitMatrix <- geoCorrection(transitMatrix, scl = TRUE)
      # Save/Load transition matrix
      saveRDS(transitMatrix, file = "transitMatrix.rds")
      
    } else {
      transitMatrix <- readRDS(file = "transitMatrix.rds")
    }
    
    
    # Define points using correct projection
    point1 <- SpatialPoints(cbind(samplelocation$Longitude, samplelocation$Latitude), proj4string = CRS(proj4string(r)))
    point2 <- SpatialPoints(cbind(OccurrenceData[row, 2], OccurrenceData[row, 3]), proj4string = CRS(proj4string(r)))
    
    # Check if the OccurrenceData point is on land, if so, skip this iteration
    if (!is.na(raster::extract(r, point2))) {
      add_error_message(paste("Point on land detected for ", species_name, " at ", OccurrenceData[row, 2], "", OccurrenceData[row, 3]))
      cat("Point on land detected for species at ", OccurrenceData[row, 2], OccurrenceData[row, 3], "\n")
      sea_distances <- append(sea_distances, Inf)
      next
    }

    sea_distance <- tryCatch({
      # Coerce points to SpatialPointsDataFrame for compatibility with gdistance
      point1_df <- SpatialPointsDataFrame(coords = point1, data = data.frame(id = 1), proj4string = CRS(proj4string(r)))
      point2_df <- SpatialPointsDataFrame(coords = point2, data = data.frame(id = 2), proj4string = CRS(proj4string(r)))
      # Compute cost distance
      cost_distance <- costDistance(transitMatrix, point1_df, point2_df)
      # Calculate the shortest path
      shortest_path <- shortestPath(transitMatrix, point1_df, point2_df, output = "SpatialLines")
      
      # Plotting the shortest path and the world map
      #plot(r, main = "Shortest Water Path")
      #plot(world, add = TRUE, col = "grey")
      #plot(shortest_path, add = TRUE, col = "blue", lwd = 2)
      #points(point1_df, col = "red", pch = 20)
      #points(point2_df, col = "green", pch = 20)
  
      # Assuming 'shortest_path' is your SpatialLines object from the shortestPath function
      # First, ensure the CRS is set on the original SpatialLines object
      crs_info <- proj4string(shortest_path)  # or use crs(shortest_path) if using `sp`
  
      # If it's not set, set it here, assuming the original data was in WGS 84 (EPSG:4326)
      if (is.na(crs_info)) {
        proj4string(shortest_path) <- CRS("+init=epsg:4326")
      }
  
      # Convert SpatialLines to sf object
      shortest_path_sf <- st_as_sf(shortest_path)
      
      # Confirm CRS is set for sf object, if not, set it:
      if (is.na(st_crs(shortest_path_sf))) {
        st_crs(shortest_path_sf) <- 4326  # EPSG code for WGS 84
      }
  
      # Transform to a suitable projected CRS for distance calculation (e.g., UTM zone 33N)
      shortest_path_utm <- st_transform(shortest_path_sf, 32633)  # UTM zone 33N

      # Calculate the length in meters
      path_length <- st_length(shortest_path_utm)
      path_length <- as.numeric(path_length)
      
      # Print the length
      print(paste0("distance through sea in m: ", path_length))
      sea_distances <- append(sea_distances, path_length)
      
    }, error = function(e) {
      add_error_message(paste("An error occurred during sea distance calculation for", species_name, "in", species_location, ":", e$message))
      return(NA)
    
    }) # trycatch() closed
  } # iteration over OccurrenceData stopped
  
  # for loop to iterate again over OccurrenceData for fly distances
  for (row in 1:nrow(OccurrenceData)) {
    
    # for this calculation, longitude comes first and then latitude!!!
    # Define points using correct projection
    point1 <- SpatialPoints(cbind(samplelocation$Longitude, samplelocation$Latitude), proj4string = CRS(proj4string(r)))
    point2 <- SpatialPoints(cbind(OccurrenceData[row, 2], OccurrenceData[row, 3]), proj4string = CRS(proj4string(r)))
    
    # Check if the OccurrenceData point is on land
    if (!is.na(raster::extract(r, point2))) {
      add_error_message(paste("Point on land detected for ", species_name, " at ", OccurrenceData[row, 2], "", OccurrenceData[row, 3]))
      flying_distances <- append(flying_distances, Inf)
      next
    }
    
    #####################
    ## FLYING DISTANCE ##
    #####################
    
    # Calculate the straight-line distance (accounting for the Earth's curvature)
    straight_line_distance <- distHaversine(coordinates(point1), coordinates(point2))
    print(paste("Straight line distance:", straight_line_distance, "meters"))
    flying_distances <- append(flying_distances, straight_line_distance)
    
  } # iteration over OccurrenceData stopped
  
  # Create data frame if lengths match
  create_data_frame <- function(distances, year, month, country) {
    if (length(distances) == length(year) && length(year) == length(month) && length(month) == length(country)) {
      return(data.frame(
        x = distances,
        year = year,
        month = month,
        country = country
      ))
    } else {
      stop("Lengths of vectors do not match. Please ensure all vectors have the same length.")
    }
  }
  ### CHECKS IF NECESSARY ###
  #cat("Length of flying_distances: ", length(flying_distances), "\n")
  #cat("content of flying_distances: ", flying_distances, "\n")
  #cat("Length of sea_distances: ", length(sea_distances), "\n")
  #cat("content of sea_distances: ", sea_distances, "\n")
  #cat("Length of OccurrenceData$year: ", length(OccurrenceData$year), "\n")
  #cat("Length of OccurrenceData$month: ", length(OccurrenceData$month), "\n")
  #cat("Length of OccurrenceData$country: ", length(OccurrenceData$country), "\n")
  
  # Create sea data frame
  sea_data <- create_data_frame(sea_distances, OccurrenceData$year, OccurrenceData$month, OccurrenceData$country)
  
  # Create fly data frame
  fly_data <- create_data_frame(flying_distances, OccurrenceData$year, OccurrenceData$month, OccurrenceData$country)
  
  # Define file paths
  sea_distance_file <- paste0("Output_calculations/sea_distances/", species_name, "_distancesTo_", species_location, ".csv")
  fly_distance_file <- paste0("Output_calculations/fly_distances/", species_name, "_distancesTo_", species_location, ".csv")
  
  # Create directories if they do not exist
  if (!dir.exists("Output_calculations/sea_distances")) {
    dir.create("Output_calculations/sea_distances", recursive = TRUE)
  }
  
  if (!dir.exists("Output_calculations/fly_distances")) {
    dir.create("Output_calculations/fly_distances", recursive = TRUE)
  }
  
  # Write data frames to CSV files
  write.csv(sea_data, file = sea_distance_file, row.names = FALSE)
  write.csv(fly_data, file = fly_distance_file, row.names = FALSE)
  # Return the result and the list of error messages
  list(result = list(sea_distances = sea_distances, flying_distances = flying_distances), error_messages = error_messages)
  
}

results <- lapply(seq_len(nrow(long)), function(i) Calculation_seadistance(long$Specieslist[i], long$name[i]))

