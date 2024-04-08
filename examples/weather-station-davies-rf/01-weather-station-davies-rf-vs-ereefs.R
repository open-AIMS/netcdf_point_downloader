# This script compares the weather station data with eReefs GBR1 Hydro 2.0,
# the BARRA2 models and the IMOS GHRSST remote sensing satellite sea surfact temperature.
# The weather station data was downloaded separately.
# This tool is relatively slow at downloading data and so it not suitable for
# downloading long and detailed time series data. For this the eReefs Data Extraction
# tool is probably more suitable https://extraction.ereefs.aims.gov.au/. The
# ereefs_points_downloader will work however for as much as you are willing to wait,
# with a practical upper limit of about 50,000 data points.

TEMP_PATH <- 'temp'
WS_WIND_PATH <- paste0(TEMP_PATH,'/Davies Reef Weather Station Wind Speed (Scalar avg 10 min) from 2020-01-01 to 2020-03-01.csv')
WS_WATER_TEMP_PATH <- paste0(TEMP_PATH,'/Davies Reef Weather Station Water Temperature @4m from 2020-01-01 to 2020-03-01.csv')
WS_AIR_TEMP_PATH <- paste0(TEMP_PATH,'/Davies Reef Weather Station Air Temperature from 2020-01-01 to 2020-03-01.csv')

# Path to the zip file
zip_file_path <- "567689ec-e084-4a0f-98d2-2045ade51d85.zip"

# Unzip the file
unzip(zipfile = zip_file_path, exdir = TEMP_PATH)

load_data <- function(path) {
  # Load the data.
  data <- read.csv(path)

  # Convert the 'time' column to Date format
  data$time <- as.POSIXct(data$time, format="%Y-%m-%d %H:%M:%S", tz="UTC")

  # Subset data to include only observations on the hour
  data <- data[format(data$time, "%M") == "00", ]
  return (data)
}

ws_wind_speed_data <- load_data(WS_WIND_PATH)
ws_water_temp_data <- load_data(WS_WATER_TEMP_PATH)
ws_air_temp_data <- load_data(WS_AIR_TEMP_PATH)



# -----------------------------------------
# Let's fetch the eReefs data from GBR1
source('../../netcdf_points_downloader.R')

# Ensure the temporary directory exists
if (!dir.exists(TEMP_PATH)) {
  dir.create(TEMP_PATH, recursive = TRUE)
}

# We can find the names of the variables from here:
# https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml.html
ereefs_nci_gbr1_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'

verbosity <- 2  # Print all the messages


# ----------- eReefs GBR1.0 NCI Data Service -------------
# Process the wind and the temperature separately. This is more convenient because
# the wind and temp data are in separate rows, and the get_ereefs_by_loc_time_depth
# adds the new attributes to each row as new columns. It doesn't make sense to
# add the ereefs temp data to the wind rows.
# There is also no guarantee that the sample times for the temperatue and wind
# match up with no gaps.
# Here we rely on the auto column matching to pick up the appropriate
# latitude, longitude, date time columns and date time format. If they hadn't
# worked then we could have specified them manually.
temp_csv_path <- paste0(TEMP_PATH,'/davies_ereefs_temp.csv')
ws_water_temp_data$depth <- -4   # The data needs a depth value for the temperature variable
cat("---------- eReefs temperature ------------\n")
water_temp_with_ereefs <- netcdf_points_downloader(
    var_names=c('temp'), 
    ws_water_temp_data, ereefs_nci_gbr1_url, temp_csv_path, verbosity, stop_on_error = TRUE)

cat("---------- eReefs wind ------------\n")
wind_csv_path <- paste0(TEMP_PATH,'/davies_ereefs_wind.csv')
wind_with_ereefs <- netcdf_points_downloader(
    var_names=c('wspeed_u', 'wspeed_v'), 
    ws_wind_speed_data, ereefs_nci_gbr1_url, wind_csv_path, verbosity, stop_on_error = TRUE)
    
# Calculate eReefs wind speed magnitude and convert from m/s to 
# km/hr to match the weather station data
wind_with_ereefs$ereefs_wind <- sqrt(wind_with_ereefs$wspeed_u^2 + wind_with_ereefs$wspeed_v^2)*3600/1000

# ---------- eReefs AIMS Daily mean, time series aggregate data service ----------
# The goal is to check whether the results seem sensible compared with the original
# hourly data.
ereefs_aims_series_url <- 'https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc'
water_temp_ereefs_mean_series <- netcdf_points_downloader(
    var_names=c('temp'), 
    ws_water_temp_data, ereefs_aims_series_url, 
    output_csv_path = paste0(TEMP_PATH,'/davies_ereefs-mean_series_temp.csv'), 
    verbosity, stop_on_error = TRUE)
    
# ---------- eReefs AIMS Daily mean, individual files -------------
# This corresponds to the daily aggregate regridded data files, not combined into
# a single time series. 
# We are doing this to test whether this configuration of individual files works with this library
ereefs_aims_indiv_url <- 'https://thredds.ereefs.aims.gov.au/thredds/dodsC/ereefs/gbr1_2.0/daily-daily/EREEFS_AIMS-CSIRO_gbr1_2.0_hydro_daily-daily-%(year)-%(month)-%(day).nc'
water_temp_ereefs_mean_indiv <- netcdf_points_downloader(
    var_names=c('temp'), 
    ws_water_temp_data, ereefs_aims_indiv_url, 
    output_csv_path = paste0(TEMP_PATH,'/davies_ereefs-mean_indiv_temp.csv'),
    verbosity, stop_on_error = TRUE)

# Accessing the data through individual files, or the time series service should not change the data
if (identical(water_temp_ereefs_mean_series$temp, water_temp_ereefs_mean_indiv$temp)) {
  stop("The water_temp_ereefs_mean_series was not the same as water_temp_ereefs_mean_indiv")
} else {
  print("** water_temp_ereefs_mean_series and  water_temp_ereefs_mean_indiv match")
}


# ----------- BARRA2 data wind --------------
# Here we are going for the surface winds
# https://dapds00.nci.org.au/thredds/catalog/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/sfcWind/v20231001/catalog.html
# http://www.bom.gov.au/research/publications/researchreports/BRR-067.pdf
# This services doesn't have an aggregate data service on THREDDS, but individual
# end points for each file. We therefore need to connect to the dates using a 
# template URL.

cat("---------- BARRA2 wind ------------\n")
barra2_url_template <- 'https://dapds00.nci.org.au/thredds/dodsC/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/sfcWind/v20231001/sfcWind_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_%(year)%(month)-%(year)%(month).nc'
temp_csv_path <- paste0(TEMP_PATH,'/davies_barra2_wind.csv')
barra2_wind <- netcdf_points_downloader(
    var_names=c('sfcWind'), 
    ws_water_temp_data, barra2_url_template, temp_csv_path, verbosity, stop_on_error = TRUE,
    timezone_download_data_src = 'UTC')

# Convert to km/hr
barra2_wind$wind_kmhr <- barra2_wind$sfcWind *3600/1000

# ----------- BARRA2 data air temp --------------
# Here we are going for the surface temperature
# https://dapds00.nci.org.au/thredds/catalog/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/ts/v20231001/catalog.html
# Turns out that the air temperature, in its variable metadata in the NetCDF file has the following:
# 202002:
# add_offset: 268.0
# scale_factor: 2.44140625E-4
# 202003:
# add_offset: 251.5
# scale_factor: 2.44140625E-4
cat("---------- BARRA air temp ------------\n")
barra2_url_template <- 'https://dapds00.nci.org.au/thredds/dodsC/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/ts/v20231001/ts_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_%(year)%(month)-%(year)%(month).nc'
temp_csv_path <- paste0(TEMP_PATH,'/davies_barra2_air-temp.csv')
barra2_air <- netcdf_points_downloader(
    var_names=c('ts'), 
    ws_air_temp_data, barra2_url_template, temp_csv_path, verbosity, stop_on_error = TRUE,
    timezone_download_data_src = 'UTC')

barra2_air$ts_celcius <- barra2_air$ts - 273.1


# ---------- IMOS GHRSST L3SM-1d ----------
# https://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/SST/ghrsst/L3SM-1d/day/2020/catalog.html
ghrsst_url_template <- 'https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3SM-1d/day/%(year)/%(year)%(month)%(day)032000-ABOM-L3S_GHRSST-SSTskin-MultiSensor-1d_day.nc'
temp_csv_path <- paste0(TEMP_PATH,'/davies_ghrsst_sst.csv')
ghrsst_sst <- netcdf_points_downloader(
    var_names=c('sea_surface_temperature'), 
    ws_water_temp_data, ghrsst_url_template, temp_csv_path, verbosity, stop_on_error = TRUE,
    timezone_download_data_src = 'UTC')
    



ghrsst_sst$sst_celcius <- ghrsst_sst$sea_surface_temperature - 273.1

# -------- Plot Water Temperature ---------
# Start capturing plot output to a PNG file
png("water-temp-ws-davies-ereefs-ghrsst-temp.png", width=2000, height=1200, res=200)

# Determine the combined range of both temperature datasets
combined_temp_range <- range(c(water_temp_with_ereefs$raw_value, water_temp_with_ereefs$temp, ghrsst_sst$sst_celcius), na.rm = TRUE)


# Initial Plot Setup with the first data series
# This creates the primary plot including the main x and y axes
plot(water_temp_with_ereefs$time, water_temp_with_ereefs$raw_value, type="l", col="blue", 
     xlab="Time", ylab="Temperature (°C)", ylim=combined_temp_range, main="Water Temperature Comparison")

# Overlay the second data series
lines(water_temp_with_ereefs$time, water_temp_with_ereefs$temp, type="l", col="purple")

lines(ghrsst_sst$time, ghrsst_sst$sst_celcius, type="l", col="orange")

lines(water_temp_ereefs_mean_series$time, water_temp_ereefs_mean_series$temp, type="l", col="green")

# Adding a right-side y-axis for the second data series (if scales differ)
#axis(side=4, at=pretty(range(c(water_temp_with_ereefs$raw_value, water_temp_with_ereefs$temp))))
#mtext("Temperature (°C)", side=4, line=3)

# Adjusting legend to accurately reflect plotted data
legend("topright", legend=c("Water Temp (Station)", "Water Temp (eReefs GBR1 hr)", "SST (GHRSST L3SM-1d)", "Water Temp (eReefs GBR1 daily)"), 
       col=c("blue", "purple", "orange", "green"), lty=1, cex=0.8)
       
# Close the PNG device, saving the file
dev.off()


# -------- Plot Air Temperature ---------
# Start capturing plot output to a PNG file
png("air-temp-ws-davies-barra2.png", width=2000, height=1200, res=200)

# Determine the combined range of both temperature datasets
combined_temp_range <- range(c(barra2_air$raw_value, barra2_air$ts_celcius), na.rm = TRUE)

# Initial Plot Setup with the first data series
# This creates the primary plot including the main x and y axes
plot(barra2_air$time, barra2_air$raw_value, type="l", col="blue", 
     xlab="Time", ylab="Temperature (°C)", ylim=combined_temp_range, main="Air Temperature Comparison")

# Overlay the second data series
lines(barra2_air$time, barra2_air$ts_celcius, type="l", col="purple")


# Adjusting legend to accurately reflect plotted data
legend("topright", legend=c("Air Temperature (Station)", "Air Temperature (BARRA2)"), 
       col=c("blue", "purple"), lty=1, cex=0.8)
       
# Close the PNG device, saving the file
dev.off()


# -------- Plot Wind ---------
png("weather-st-davies-vs-ereefs-vs-barra2-wind.png", width=2000, height=1200, res=200)


# Plotting
#par(mar=c(5, 4, 4, 4) + 0.3)  # Adjust margins to accommodate two y-axes

# Plot Wind Speed (using left y-axis)
plot(wind_with_ereefs$time, wind_with_ereefs$raw_value, type="l", col="red", 
     xlab="Time", ylab="Wind Speed (Scalar avg 10 min)", 
     main="Wind Speed and Water Temperature over Time at Davies Reef")
lines(wind_with_ereefs$time, wind_with_ereefs$ereefs_wind, type="l", col="orange")
lines(barra2_wind$time, barra2_wind$wind_kmhr, type="l", col="purple")

# Adjusting legend to reflect all plotted data
legend("topright", legend=c("Wind Speed (Station)", "Wind Speed (eReefs)", "Wind Speed (BARRA2)"),
       col=c("red", "orange", "purple"), lty=1, cex=0.8)
       
# Close the PNG device, saving the file
dev.off()


# -------- Plot Cross correlation ---------
# Plot the cross correlation between the weather station data and the eReefs
# data to determine the best alignment. If there is a time zone bug then
# the data should be about 10 hours offset.
png("ws-ereefs-cross-correlation.png", width=2000, height=1200, res=200)
ccf(water_temp_with_ereefs$raw_value, water_temp_with_ereefs$temp, lag.max = 32, main = "Cross correlation between weather station and eReefs")
dev.off()