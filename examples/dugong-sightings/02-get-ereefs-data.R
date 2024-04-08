# In this example we demonstrate extracting eReefs data for simulated 
# dugong sightings. This is a random subset and randomisation (time and
# space) of observations for demonstration purposes only. The original
# data is not made available in this code repository, only the example
# simulated data.

# Just need basic map with minimal dependencies to check the output
if (!require(maps)) {
  install.packages("maps")
  library(maps)
} else {
  message("The 'maps' package is already installed.")
}

source('../../netcdf_points_downloader.R')

# We want the finest time and spatial resolution possible and so we are
# going to use the GBR1 Hydro 2.0 model available from NCI. This gives
# us hourly data.
opendap_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'

# File with the location and times to extract data from eReefs.
input_csv_path <- 'sim-dugong-sightings-nGBR2023.csv'

# Choose a path to save the final combined data. This path is also used
# to save the cache RDATA file.
output_csv_path <- paste0('temp/sim-dugong-nGBR2023_gbr1_hydro_v2.csv')

# We want the tidal information
variables = c('eta')

verbosity = 2 # Want to see all the progress messages. Set to 0 for no messages.

# Download the points. This will be slow for lots of points, but progress
# is saved and the script will resume from where it left off.
# Both eReefs and the Dugong data is in Aust/Brisbane UTC+10 timezone.
# This is the default so we don't need to specify this. 
data_with_ereefs <- netcdf_points_downloader(
    variables, input_csv_path, opendap_url, output_csv_path, verbosity)
    


# --------------- Plot the results -----------------
# Start capturing plot output to a PNG file
png("sim-dugong-sightings-map.png", width=1000, height=1500, res=200)

# Define colors from blue (for low values) to red (for high values)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Scale eta values to the range of 1 to 100 to match the colors vector
scaled_eta <- round(((data_with_ereefs$eta + 2) / 4) * 99) + 1

# Handle NA values in eta by assigning a default color (e.g., grey)
colors_with_na <- ifelse(is.na(data_with_ereefs$eta), "grey", colors[scaled_eta])

# Basic map of the area of interest, e.g., Great Barrier Reef or a broader area if necessary
map("world", region = "Australia", 
  xlim = c(min(data_with_ereefs$long), max(data_with_ereefs$long)), 
  ylim = c(min(data_with_ereefs$lat), max(data_with_ereefs$lat)))

# Add points to the map with colors based on eta value
points(data_with_ereefs$long, data_with_ereefs$lat, col = colors_with_na, pch = 20)


# Define legend colors based on selected positions in the color ramp
legend_colors <- c(colors[1], colors[25], colors[50], colors[75], colors[100])

# Define legend labels
legend_labels <- c("<-1.5", "-1.5 to -0.5", "-0.5 to 0.5", "0.5 to 1.5", ">1.5")

# Add a legend to the plot
legend("topright", legend = legend_labels, fill = legend_colors, title = "Eta Value")

# Adding a title
title("Tide level for simulated Dugong Sightings")

# Add a border around the plot
box()

# Close the PNG device, saving the file
dev.off()
