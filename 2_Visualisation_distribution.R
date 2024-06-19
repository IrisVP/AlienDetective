####################################################
### LOAD LIBRARIES

library("tidyr")
library("ggplot2")
library("dplyr")
library("poliscidata")
####################################################

####################################################
### PREPARE LOCATION FILE
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# first load in species_location file
species_location <- read.csv("Output_Preparation/Species_Location.csv")
# Reshape data from wide to long format
long <- pivot_longer(species_location, !Specieslist)
# Filter rows where 'value' is greater than 0
long <- long[long$value > 0, ]

####################################################

# iterate over species_location file again using function and map()
# and read csv file per species

Distribution_seadistance <- function(species_name, species_location){
  # read csv file per species
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  distance_file <- paste0("Output_calculations/sea_distances/", species_name, "_distancesTo_", species_location, ".csv")
  
  distances <- read.csv(distance_file, header = TRUE)
  # clean dataframe from rows with Inf in them
  distances <- distances[is.finite(distances$x), ]
  # calculates km
  distances$x <- distances$x/1000
  # select distances
  distances <- subset(distances, x < 40000)
  
  ###############################################################################
  # make histograms of distances per species, with filtering on distance limit 40000
  ###############################################################################
  dist_plot <- ggplot(distances, aes(x = x)) +
    geom_histogram(binwidth = 50, fill = "blue", color = "black", boundary = 0) +
    labs(title = "Histogram of Distances", x = "Distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement with hjust
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
    
  #print(dist_plot)
  # create a directory to save the plots in
  
  if(!dir.exists("Output_calculations/sea_distribution_plots")) {
    dir.create("Output_calculations/sea_distribution_plots")
    ggsave(filename = paste0("Output_calculations/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("Output_calculations/sea_distribution_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = dist_plot, width = 50, height = 30, units = "cm", dpi = 900)
  }
}
################################# END FUNCTION

# error lists initiations to store error messages in
error <- c()

# executing iteration for the long file!
plot <- Map(Distribution_seadistance, long$Specieslist, long$name)

# write error files
print(error)

file_error <- paste0("Output_calculations/errors_graph_sea_distances.csv",)
write.table(error, file = file_error, append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)



#############################################################
# Make combined histogram of sea distances and fly distances
#############################################################

Distribution_combDistance <- function(species_name, species_location){
  # make variables with filenames in which distances are saved
  sea_distance_file <- paste0("Output_calculations/sea_distances/", species_name, "_distancesTo_", species_location, ".csv")
  fly_distance_file <- paste0("Output_calculations/fly_distances/", species_name, "_distancesTo_", species_location, ".csv")
  
  # read the files
  sea_distances <- read.csv(sea_distance_file, header = TRUE)
  fly_distances <- read.csv(fly_distance_file, header = TRUE)
  # transform Inf to NA
  sea_distances <- sea_distances[is.finite(sea_distances$x),]
  fly_distances <- fly_distances[is.finite(fly_distances$x),]
  
  # convert each dataframe from m to km and select the distances smaller than 40000km
  sea_distances$x <- sea_distances$x/1000
  sea_distances <- subset(sea_distances, x < 40000)
  sea_distances$type <- "sea_distance" # type will be a column in combined dataframe
  fly_distances$x <- fly_distances$x/1000
  fly_distances <- subset(fly_distances, x < 40000)
  fly_distances$type <- "fly_distance" # type will be a column in combined dataframe
  
  # make a combined dataframe to use for the graph
  combined_distances <- rbind(sea_distances, fly_distances)
  
  p <- combined_distances %>%
    ggplot( aes(x=x, fill=type)) +
    geom_histogram(binwidth = 50, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species_name, " from ", species_location)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 16),  # Set font size for axis numbers
      axis.title = element_text(size = 20), # Set font size for title
      legend.title = element_text(size = 18, face="bold"), # Settings for legend title
      legend.text = element_text(size = 16)) +  # settings for legend text
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +  # settings for x axis
    scale_y_continuous(expand = c(0, 0)) +
    # used expand to make sure the axes are on the lines of the axes and not above them floating
    coord_cartesian(xlim = c(0, 3000)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    scale_fill_discrete(
      name = "Distance type",
      breaks = c("sea_distance", "fly_distance"),
      labels = c("Sea distance", "Fly distance")) +
    labs(x = "Distance in km", y = "Frequency")
  
  #print(p)
  
  # check if directory exists, save the plot
  if(!dir.exists(paste0("Output_calculations/combined_distribution_plots"))) {
    dir.create(paste0("Output_calculations/combined_distribution_plots"))
    ggsave(filename = paste0("Output_calculations/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("Output_calculations/combined_distribution_plots/plotcomb_", species_name, "_from_", species_location,".png"), 
           plot = p, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Distribution_combDistance, long$Specieslist, long$name)

#############################################################
# Make histograms of locations
#############################################################

Location_histograms <- function(species_name, species_location){
  
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  # make variable with filename
  distance_file <- paste0("Output_calculations/sea_distances/", species_name, "_distancesTo_", species_location, ".csv")
  # read csv file per species
  distance_file <- read.csv(distance_file, header = TRUE)
  # convert to meters
  distance_file$x <- distance_file$x/1000
  # clean dataframe from rows with Inf in them
  distance_file <- distance_file[is.finite(distance_file$x), ]
  # select only distances below 40000km
  distances <- subset(distance_file, x < 40000)
  
  ### PLOT
  
  country_plot <- ggplot(distances, aes(x = x, fill = country)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species_name," in ", species_location),
                        x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  #print(country_plot)
  # create a directory to save the plots in
  if(!dir.exists("Output_calculations/sea_distribution_country_plots")) {
    dir.create("Output_calculations/sea_distribution_country_plots")
    ggsave(filename = paste0("theoretical_data/plot_", species_name, "_from_", species_location, ".png"), 
           plot = country_plot, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("Output_calculations/sea_distribution_country_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = country_plot, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Location_histograms, long$Specieslist, long$name)



#############################################################
# Make histograms of year categories
#############################################################

Year_histograms <- function(species_name, species_location){
  
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))
  # make variable with filename
  distance_file <- paste0("Output_calculations/sea_distances/", species_name,
                          "_distancesTo_", species_location, ".csv")
  
  # read the csv files
  distance_file <- read.csv(distance_file, header = TRUE)
  # convert to meters
  distance_file$x <- distance_file$x/1000
  # clean dataframe from rows with Inf in them
  distance_file <- distance_file[is.finite(distance_file$x), ]
  distances <- subset(distance_file, x < 40000)
  # make a function which assigns years to specific year categories
  assign_category <- function(year) {
    if (is.na(year)) {
      return(NA)   # return NA when year is not present
    }
    for (category in year_categories) {
      range <- as.numeric(unlist(strsplit(category, "-"))) # save years as numeric without "-"
      if (year >= range[1] & year < range[2]) {  # if the year falls into this category
        return(category) # return this category
      }
    }
    return(NA) # If year doesn't fall into any category, return NA
  }
  # make a variable with preferenced year categories
  year_categories <- c("1965-1985","1985-1990",
                       "1990-1995","1995-2000","2000-2005","2005-2010","2010-2015",
                       "2015-2020","2020-2025")
  # Apply function to create new column with year categories
  distances$year_category <- sapply(distances$year, assign_category)
  
  #################################
  ### PLOT
  
  year_plot <- ggplot(distances, aes(x = x, fill = year_category)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species_name," in ", species_location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "YlOrRd", na.value = "black") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 16),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # set legend title size
          legend.text = element_text(size = 12),    # set legend text size
          legend.key.size = unit(1.5, "lines")) +   # set legend key size
    scale_x_continuous(breaks = seq(0, 3000, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 3000)) # Use coord_cartesian for setting limits
  
  #print(year_plot)
  # create a directory to save the plots in
  if(!dir.exists("Output_calculations/sea_distribution_year_plots")) {
    dir.create("Output_calculations/sea_distribution_year_plots")
    ggsave(filename = paste0("Output_calculations/sea_distribution_year_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = year_plot, width = 60, height = 30, units = "cm", dpi = 900)
  } else {
    ggsave(filename = paste0("Output_calculations/sea_distribution_year_plots/plot_", species_name, "_from_", species_location, ".png"), 
           plot = year_plot, width = 60, height = 30, units = "cm", dpi = 900)
  }
}

# executing iteration for the long file!
plot <- Map(Year_histograms, long$Specieslist, long$name)

