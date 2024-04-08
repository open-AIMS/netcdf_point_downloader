# In this example we compare the tide levels calculated in eReefs with the
# tide gauge in Townsville. This tests whether there is good alignment with
# time zones. 
TIDE_GAUGE <- 'QLD-Gov_TSV-near-real-time-gauge/Townsville_tides.txt'
TEMP_PATH <- 'temp'

# This is the samples that we will extract data from eReefs from.
# We create this CSV file using this script at the start.
# We only create the CSV file once if it doesn't already exists. 
# 
input_csv_path <- paste0(TEMP_PATH,'/tides_30min.csv')

if (!file.exists(input_csv_path)) {
  # Tide gauge data is in local time, as verified against Willy Weather
  # The input file uses double spaces for separators for most columns. When
  # we load the data we end up with columns filled with NA.
  tidesRaw <- read.csv(file = TIDE_GAUGE, sep=' ', header=FALSE, stringsAsFactors = FALSE, skip=18)

  # Just keep the columns of interest.
  tides <- data.frame(dateTime = as.POSIXct(tidesRaw[,1],format = '%d%m%Y%H%M'), tideHeight=tidesRaw[,3])

  # Function to check if a time is at a 30-minute interval
  is_half_hour <- function(time) {
    as.integer(format(time, "%M")) %% 30 == 0
  }

  # Filter rows that fall exactly on 30-minute intervals
  tides_30min <- tides[sapply(tides$dateTime, is_half_hour), ]

  # Filter out any rows with NA dateTime values before proceeding
  tides_30min <- tides_30min[!is.na(tides_30min$dateTime), ]

  # Plotting
  plot(tides_30min$dateTime, tides_30min$tideHeight, type = 'l', xlab = "Time", ylab = "Tide Height (m)", 
    main = "Tide Height Over Time at 30-minute Intervals")

  # Format the dateTime column to ensure consistent datetime format
  # Without this times corresponding to mid night would have the 
  # time removed, saving just the date. This would make the date time
  # format inconsistent.
  tides_30min$dateTime <- format(tides_30min$dateTime, "%Y-%m-%d %H:%M:%S")

  # Location of the tide gauge. This is covered by the GBR1 model
  tides_30min$lat <- -19.25
  tides_30min$lon <- 146.833333

  # Ensure the temporary directory exists
  if (!dir.exists(TEMP_PATH)) {
    dir.create(TEMP_PATH, recursive = TRUE)
  }

  # Each time we update this file get_ereefs_by_loc_time_depth will notice that
  # the modification date is different and it will redownload the data. We
  # therefore ensure that we only calculate this CSV file once.
  write.csv(tides_30min, input_csv_path, row.names = FALSE)
} else {
  print("Already created tide data, skipping generation")
}


source('../../netcdf_points_downloader.R')

# We want hourly hydro data and so we need the raw model data
# on NCI
opendap_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'

# Saved output after download, cache is also saved based on same path.
output_csv_path <- paste0(TEMP_PATH,'/tides_30min_ereefs.csv')

# Tides (model is Mean Sea Level)
variables = c('eta')

# Print all messages
verbosity = 2

# Here we read in the townsville tide data from the generated CSV file.
# Alternatively we could have skipped this and just passed in the tides_30min 
# data frame. We are doing it here just to test this feature.
tides_with_ereefs <- netcdf_points_downloader(
    variables, input_csv_path, opendap_url, output_csv_path, verbosity, 
      timezone_input_data_src = 'Etc/GMT-10',
      timezone_download_data_src = 'Etc/GMT-10')


# --------  Plot the results -----------------
# Start capturing plot output to a PNG file
png("tsv-tide-gauge-vs-ereefs.png", width=2000, height=1400, res=200)

# Convert dateTime from string to POSIXct if not already (assumed necessary conversion)
tides_with_ereefs$dateTime <- as.POSIXct(tides_with_ereefs$dateTime, format = "%Y-%m-%d %H:%M:%S", tz = 'Etc/GMT-10')

# Adjust to LAT as eReefs is Mean Sea Level. I used the closest listed (Magnetic Island) at:
# https://www.qld.gov.au/environment/coasts-waterways/beach/storm/storm-sites/townsville
tides_with_ereefs$eta_adj = tides_with_ereefs$eta + 1.92

# Plot tideHeight
plot(tides_with_ereefs$dateTime, tides_with_ereefs$tideHeight, type = 'l', col = 'blue', 
     xlab = 'DateTime', ylab = 'Tide (m LAT)', 
     main = 'Tide Height (Station 100447) and Eta (GBR1 Hydro v2.0) for Townsville')

# Add eta to the existing plot
lines(tides_with_ereefs$dateTime, tides_with_ereefs$eta_adj, col = 'red')

# Add a legend
legend("topright", legend = c("Tide Height", "Eta (eReefs tide)"), col = c("blue", "red"), lty = 1)

# Close the PNG device, saving the file
dev.off()