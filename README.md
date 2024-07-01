# AlienDetective
Species of maritime fauna all over the world are known to travel great distances in the oceans and seas. The goal of this workflow is to be able to detect them by calculating sea distances (going around land) from a sample location to the occurrence data from this species (from GBIF). This workflow focuses on occurrence data in Europe.
![logo3](https://github.com/IrisVP/AlienDetective/assets/151626670/21dd7508-bd81-448a-a096-db07bace2515)
### Input

In the 'Input' folder, there are examples of input files necessary for this workflow. 
The files are in .csv format:
- Coordinates.csv: Contains the names of the samplelocations in the column 'Observatory.ID', along with their latitudes and longitudes. <br />
- Metadata.csv: Contains metadata per sample, including sample region, country, latitudes, longitudes, sample dates, etc. Note: Metadata should be prepared for personal use in this workflow. The 'Coordinates.csv' file is an example of this preparation. <br />
- Species_Location.csv: Contains a 'Specieslist' column with all the species names. The other columns are sample location names, with the data in these columns representing the number of samples taken from a species for a specific location.<br />

Note: The Coordinates.csv and Species_Location.csv files should be prepared exactly like the examples for use in 1_Calculation_sea_distances.R.

### Short description of R scripts

For detailed documentation: see the 'Documentation' directory for .Rmd files

#### 1_Calculation_sea_distances.R
/!\ Run this on a server. This is computationally expensive. <br />

1. Download and install libraries, then load data.
2. Retrieve occurrence data from GBIF, focusing on Europe. The settings for fetching data can be personalized here.
3. Calculate distances: both flying and sea distances are calculated and saved in files inside the 'Output_calculations' directory.

##### Output
The output are two directories in the 'Output_calculations' directory named 'sea_distances' and 'fly_distances'. Separate csv files are made for each organism and their corresponding samplelocation name.

#### 2_Visualisation_distribution.R

Histograms are created showing the frequencies of species occurrences per sample location and species. These graphs are saved in the 'Output_calculations' directory, with a separate directory for each type of graph:

1. The first graph shows the sea distances per species and sample location.
2. The second graph is a combined graph of flying distances and sea distances.
3. The third graph is similar to the first graph but colored by country.
4. The fourth graph is also similar to the first graph but colored by year category.

Patterns can be observed in these histograms. The x-axis represents distances in km, and the y-axis represents frequency. An alien species can be detected when there is a spike in frequency at a large distance. Using the country plot, the origin region of the species can be identified. Using the year plot, the year in which the species migrated to that location can be determined.

#### 1_Preparation.R

This script is located in the directory called 'Preparation_BOLDigger'. This script can be used when files from the BOLDigger workflow are used. Be sure to check if the headers of your files are the same as the example BOLDigger file in 'Input'. 

### /!\ Warnings /!\

The distance to a location inside a bay or canal next to the coast line will most likely not be calculated. The world map used in the calculation does not have a good enough resolution to differentiate bays or canals from land areas. If you want calculation to these locations, place the coordinates outwards in the sea or ocean. Fresh water sources are also not included in these calculations.
