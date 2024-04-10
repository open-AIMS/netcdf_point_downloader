# Copyright 2024 Australian Institute of Marine Science - Eric Lawrey
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Version history
# 1.0   2024-04-08  Initial release.

library(ncdf4)
library(tools)

# This script uses a data structure to keep track of the downloading of the data.
# Parts of this structure are persisted to disk so that the scripts can be resumed
# and the script can pick up from where it left off.
# This script is also tolerant of various naming of input data attributes. This is
# so that it can more easily be used with existing CSV files and cope with changes
# to date formats that Excel tends to do with CSV.
# The data_cache also contains a record of all the matching of column names and
# date formatting.
# data_cache has the following structure:
# data_cache$df - data frame: corresponds to the original input locations, times and depths
#       and includes additional attributes for tracking the download progress. Persisted
#       to disk.
# data_cache$df_metadata - list: variables to keep track of. Persisted to disk. 
#     $data_cache_path      Where to store the disk cache
#     $download_data_src_url_template:         Path of the data source.
#     $var_names:           Names of the variables being downloaded
#     $input_latitude_name:     Name of the latitude column of input data
#     $input_longitude_name:    Name of the longitude column of input data
#     $input_datetime_name:     Name of the datetime column of input data
#     $input_depth_name:        Name of the depth column of input data
#     $input_datetime_format:     Format string for parsing the input datetime column
#     $download_latitude_name:     Name of the latitude column of download source
#     $download_longitude_name:    Name of the longitude column of input data
#     $download_datetime_name:     Name of the datetime column of input data
#     $download_depth_name:        Name of the depth column of input data
#     $variable_attributes: List of variables being downloaded and whether they need time or depth
#         $eta            (example)
#           $needs_time
#           $needs_depth
# data_cache$ds:            Open netcdf connection to the data source.
# data_cache$enable_disk_cache: Whether to persist the cache data to disk.

# Because we allow connecting to time series that are made up from multiple
# data end points described with a template that describes the set of URLs,
# this can present a bit of a catch 22 when combined with time zone differences.
# To align the dates of the input data with the data source URLs we need
# to know the time zone of the data sources as this will determine the matching
# date boundaries. i.e. input '2023-01-01 06:00 +10' would be in the file
# corresponding to 2022-12-31 if it is in UTC time. The problem is that we
# can't know the time zone of the data source without connecting to it.
# 
# To break this we select a date in the middle of the input time series,
# so we avoid a potential time zone issue right at the start or the end
# going off the end of the series. We then connect to that time slice
# to get key metadata. We assume that the time series is consistent, i.e.
# the structure and time zones are shared across all files.
#


# Global variables used as defaults for the names search for in the input data source
# and the download data source. These can be overriden to handle cases outside
# this default list.

netcdf_points_downloader_globals <- list ( 
  input_lat_names = c('lat', 'latitude', 'LAT', 'LATITUDE', 'Y', 'y', 'i'),
  input_lon_names = c('long', 'longitude', 'lon','LONG', 'LONGITUDE', 'LON', 'X', 'x', 'j'),
  input_datetime_names = c('datetime', 'dateTime', 'date', 'Date', 'time', 'Time'),
  input_depth_names = c('depth', 'zc', 'Depth', 'depth_m', 'k'),
  # Note: as.POSIXct only understands time zones if they are specified in the 
  # follow format: +0800, -1000. Other formats are simply ignored. So for example
  # '2021-01-01 10:00 +10' the '+10' will be ignored. It is for this reason 
  # we get the user to manually specify the time zone.
  input_datetime_formats = c('%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M', 
                        '%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M',
                        '%Y/%m/%d %H:%M:%S', '%Y/%m/%d %H:%M',
                        '%d-%m-%Y %H:%M:%S', '%d-%m-%Y %H:%M', 
                        '%d/%m/%Y %H:%M:%S', '%d/%m/%Y %H:%M'),
  # Allow for the spacing between grid cells to be bigger than a regular
  # grid. This is set just big enough that for GBR1 the points between the
  # stretched cells in the north are not marked as out of bounds. The trade
  # off is that mid way down the GBR1 grid on the outer edge the checking 
  # lets through points that are up to 3 pixels past the outer bounds.
  # These points will return NA when they are downloaded from the data source.
  curvilinear_grid_squish = 4.5
  
)

# Use the same default for the download source
# These are the names that the code will look for in the download data source to
# line up variables. 
netcdf_points_downloader_globals$download_lat_names <- netcdf_points_downloader_globals$input_lat_names
netcdf_points_downloader_globals$download_lon_names <- netcdf_points_downloader_globals$input_lon_names
netcdf_points_downloader_globals$download_datetime_names <- netcdf_points_downloader_globals$input_datetime_names
netcdf_points_downloader_globals$download_depth_names <- netcdf_points_downloader_globals$input_depth_names
netcdf_points_downloader_globals$download_datetime_names <- netcdf_points_downloader_globals$input_datetime_names

# This is a manual patch for eReefs where the mapping of the dimensions to the variables
# needed is not obvious from the files. Here we associate i and j with specific
# latitude and longitudes which is strictly not quite right as for a curvilinear grid
# the i and j and indices into both the latitude and longitude grids. In this case we
# are flagging that if we see an 'i' as a dimension then we need to get the 'latitude'
# variable for the lat_array and 'j' we need to get 'longitude' variable for the 
# long_array. This works because i and j also go together and the ordering is always
# consistent. eReefs is the only example I have seen of curvilinear grids and so
# this approach may not work in a more general setting. 
netcdf_points_downloader_globals$download_dim_map <- list (
  i = 'latitude',
  j = 'longitude',
  k = 'zc'
)

# Maximum time tolerance for matching input date times to times in the NetCDF
# download source if the source uses a template to form a series. In this case
# the script cannot work out the spacing between time slices.
netcdf_points_downloader_globals$single_end_point_time_tolerance_sec <- 24*3600

# The cache stores the progressive download that allows the
# extraction to resume after an interruption.
.save_data_cache <- function(data_cache) {
  if (!data_cache$enable_disk_cache) return()

  data_cache_path <- data_cache$df_metadata$data_cache_path
  
  dir.create(dirname(data_cache_path), recursive = TRUE, showWarnings = FALSE)
  metadata_path <- sub("\\.csv$", ".rds", data_cache_path)
  
  # Create a deep copy so that we can safely remove data_cache$ds variable 
  # prior to the save without affecting the data_cache variable
  data_cache_copy <- unserialize(serialize(data_cache, connection = NULL))
  # Remove the ds component from the copy
  data_cache_copy$ds <- NULL
  saveRDS(data_cache_copy, file = metadata_path)
}


#' This function is helps with adjusting to data sources that are 
#' spread across multiple end points, rather than a single aggregate
#' service. 
#' This function adjusts the base_template to the specified date,
#' substituting in the year, month and day into the template
#' to create the URL of the data source for the specified date.
#' For example for using GBR1 eReefs without the aggregate data services
#' it has individual daily files, with end point data URLs like:
#' https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_2024-01-17.nc
#' We can change this to a template such as:
#' https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_%(year)-%(month)-%(day).nc
#' @param time_date POSIXct date object. Assumes the time zone has been set appropriately.
#' @param timezone_download_data_src string: time zone of the download source. 'ETC/GMT-10' for example.
#' @param input_data_src_template string - URL of the data source.
generate_data_source_url <- function(time_date, input_data_src_template) {
  # Format date to extract year, month, and day components
  year <- format(time_date, "%Y")
  month <- format(time_date, "%m")
  day <- format(time_date, "%d")
  # Define patterns for placeholders
  year_pattern = "%\\(year\\)"
  month_pattern = "%\\(month\\)"
  day_pattern = "%\\(day\\)"
  
  # Substitute year, month, and day if their placeholders are found
  if (grepl(year_pattern, input_data_src_template)) {
    input_data_src_template <- gsub(year_pattern, year, input_data_src_template)
  }
  if (grepl(month_pattern, input_data_src_template)) {
    input_data_src_template <- gsub(month_pattern, month, input_data_src_template)
  }
  if (grepl(day_pattern, input_data_src_template)) {
    input_data_src_template <- gsub(day_pattern, day, input_data_src_template)
  }
  
  return(input_data_src_template)
}


# Have a common open dataset function so the messaging can be consistent.
# This ensures that the connection to the data source is setup. This indicates
# if the data source changed.
# @return list(data_cache = data_cache, download_data_src_changed = download_data_src_changed)
.open_ds <- function(data_cache, row_index, verbosity) {
  
  download_data_src_url <- data_cache$df$download_data_src_url[row_index]
  if (is.na(download_data_src_url)) {
    stop("BUG: .open_ds called before input_data_src_urls setup in data_cache$df")
  }
  download_data_src_changed <- FALSE  # Initialize the flag to FALSE
  
  # Check for an existing data source connection and whether it should be reused
  if (!is.null(data_cache$ds)) {
    if (data_cache$ds$filename == download_data_src_url) {
      # The current connection matches the requested data source URL
      # Return the existing connection along with the unchanged flag
      return(list(data_cache = data_cache, download_data_src_changed = download_data_src_changed))
    } else {
      # The current connection does not match the requested data source URL
      # Close the existing connection and mark the data source as changed
      nc_close(data_cache$ds)
      download_data_src_changed <- TRUE
    }
  }
  
  # Connect to the new data source
  if (verbosity > 0) cat(sprintf("Connecting to OpenDAP: %s\n", download_data_src_url))
  data_cache$ds <- nc_open(download_data_src_url)
  
  # Return the updated data_cache and the data source change flag
  return(list(data_cache = data_cache, download_data_src_changed = download_data_src_changed))
}
  
#' If the cache files already exist then this sets up the data cache
#' from these files. If not then we set up the cache from scratch.
#' The data cache starts as a data frame with all the sampling points to
#' be downloaded. The variables to be downloaded are added as empty columns.
#' The cache also contains the calculated index into the data, i.e.
#' what row (lat_idx), col (lon_idx), depth (depth_idx) and time (time_idx) 
#' in the grid that is the nearest location to the sample data.
#' @param input_df - data frame: all the samples to download from the data dource.
#' @param data_cache_path - String: location where the cache should be saved. If enable_disk_cache
#'          is FALSE then this string is ignored.
#' @param var_names - vector of strings: names of the variables to download. For eReefs
#'          some examples are: 'temp', 'botz', 'salt'.
#' @param download_data_src_url_template - string: URL of the data source OpenDAP end point. Can also be a file path.
#'          This can be a template that specifies a range of end points based on the
#'          date.
#' @param output_csv_path - string: Path to final download. Use this to check if the whole 
#'            download is complete and we don't need to reconnect to the data source.
#' @param timezone_input_data_src - string: Time zone string of input sample data
#' @param timezone_download_data_src - string: Time zone string of data source we are
#'            downloading from.
#' @param verbosity - numeric: 0 - 2, amount of messaged printed, 0 is none 2 is all.
#' @param enable_disk_cache - boolean: If TRUE then the cache will be saved and read from disk.
#'          If FALSE then there is no effective cache as each run will forget download
#'          progress. This might still be useful if you are download a small number of 
#'          points and don't want to clutter the disk with cache files.
#' @param input_overrides - list: See netcdf_points_downloader
#' @param download_overrides - list: See netcdf_points_downloader
#' @return List: data_cache
.setup_data_cache <- function(input_df, data_cache_path, var_names, download_data_src_url_template, 
            output_csv_path, 
            timezone_input_data_src, timezone_download_data_src, 
            verbosity = 1, enable_disk_cache = TRUE) {
  
  # This loads any existing data_cache. It doesn't renew the open connection 
  # to the OpenDAP service.
  read_data_cache <- function(data_cache_path) {
    # If we can't find any of the cache files then they must not
    # have been set up yet.
    if (!file.exists(data_cache_path)) {
      return(NULL)
    }
    data_cache <- readRDS(file = data_cache_path)
    return(data_cache)
  }
  
  
  if (enable_disk_cache && (!is.null(data_cache_path))) {
    data_cache <- read_data_cache(data_cache_path)
  } else {
    # If no disk cache then start fresh each time
    data_cache <- NULL
  }
  
  # We have a cache, but it might be stale
  if (!is.null(data_cache)) {
    # Test if the input points data is the same as our cache. If it is not
    # then invalidate the cache and setup from scratch.
    # Also check that the variables to be extracted match and that the 
    # download source matches.
    
    # Get the common columns that exist in both input_df and data_cache$df
    common_columns <- intersect(names(input_df), names(data_cache$df))

    # Subset data_cache$df to only include columns that are common
    sub_data_cache_df <- data_cache$df[, common_columns]

    # Check if all values in the subset of data_cache$df match those in input_df
    data_change <- !identical(sub_data_cache_df, input_df)

    variables_change <- !identical(sort(data_cache$df_metadata$var_names), sort(var_names))
    data_url_change <- data_cache$df_metadata$download_data_src_url_template != download_data_src_url_template
    
    # If there is a change then the cache is stale. Print some helpful messages
    # indicating why the cache is invalid.

    if (data_change || variables_change || data_url_change) {
      data_cache <- NULL
      if (verbosity > 0) {
        message <- "Cache is invalid due to changes:"
        if (data_change) {  
          # Assuming sub_data_cache_df and input_df are already defined
          # Step 1: Identify different columns
          diff_cols <- sapply(names(input_df), function(col_name) {
            !identical(sub_data_cache_df[[col_name]], input_df[[col_name]])
          })

          changed_cols <- names(diff_cols[diff_cols])
          cat(sprintf("Cache is invalid due to data changes in columns: %s.\n", 
            paste(changed_cols, collapse = ", ")))

        }
        if (variables_change) {
          cat("Cache is invalid due to change in variables. Changed variables are:\n")
          # Identify which variables have changed
          changed_vars <- setdiff(union(data_cache$df_metadata$var_names, var_names), 
                                  intersect(data_cache$df_metadata$var_names, var_names))
          cat(paste(changed_vars, collapse = ", "), "\n")
        }
        if (data_url_change) {
          cat("Cache is invalid due to change in data URL.\n")
          cat("Previous URL:", data_cache$df_metadata$download_data_src_url_template, "\nNew URL:", download_data_src_url_template, "\n")
        }
      }
    }
    
  }
  
  # Invalid cache or first time, so setup the data_cache from scratch
  if (is.null(data_cache)) {
    if (verbosity > 0 && enable_disk_cache) cat("Setting up the data cache: ", data_cache_path, "\n")
    
    # Start the cache from the initial data frame
    df <- input_df
    # Initialise columns for recording progress
    for (col in var_names) {
      df[[col]] <- NA  # Assuming variables are numerical
    }
    for (col in c('download_data_src_url', 'time_idx', 'lat_idx', 'lon_idx', 'depth_idx')) {
      df[[col]] <- NA  
    }
    for (col in c('index_complete', 'vars_download')) {
      df[[col]] <- FALSE  
    }
    # Remember all the important metadata that corresponds to the this download
    # this information is used to verify if things have changed after a restart.
    # This information is saved in a R data file.
    df_metadata <- list(
      var_names = var_names, 
      data_cache_path = data_cache_path, 
      download_data_src_url_template = download_data_src_url_template)
    
    
    # Create a starting data cache
    # Record enable_disk_cache in the cache so that it knows whether to save to disk or not.
    # We don't want this variable save to disk so we don't have it in the df_metadata
    data_cache <- list(df = df, df_metadata = df_metadata, enable_disk_cache = enable_disk_cache)
    
    # Work out the alignment between column names in the input data source and internal names
    # This adds additional values to df_metadata
    data_cache = .identify_column_names_and_format(data_cache, 
        timezone_input_data_src, timezone_download_data_src, verbosity)
    #.save_data_cache(data_cache)
    
    return(data_cache)
  } else {
    # Cache is valid
    
    # Setup variables in the cache that are not saved and restored from disk
    
    # This should logically be always true as we should only get here if the 
    # disk cache is enabled
    data_cache$enable_disk_cache <- enable_disk_cache
    
    # Connect to the data source. Start with the first row of data.
    #result <- .open_ds(data_cache, 1,  verbosity)
    #data_cache <- result$data_cache
    return(data_cache)
  }
}




#' Utility function for identifying the datetime format of a set samples.
#' This is used so that the script can automatically adjust to changes in
#' input date time formats. This often can happen if a CSV is loaded and
#' manipulated in Excel.
#' This function will stop if no matching format is found.
#' 
#' @param datetime_samples String vector: sample date time values in the data 
#'        to check against.
#' @param datetime_formats: String vector: datetime formats to try and match
#'        against the datetime_samples.
#' @return the matching format for the samples.
#' @examples
#' # This will return '%Y/%m/%d'
#' identify_datetime_format(c('2024/12/08','2022/06/01'), c('%Y/%m/%d', '%Y-%m-%d'))
# Returns:
#  - The matching datetime format string. 
identify_datetime_format <- function(datetime_samples, datetime_formats) {
  
  # Take a subset of the dataset to check against 
  valid_format_found <- FALSE
  validated_format <- ""
  
  for (fmt in datetime_formats) {
    # Attempt to parse the sample dates with the current format
    parsed_dates <- as.POSIXct(datetime_samples, format = fmt)
    # If none of the parsed dates are NA, the format is likely correct
    if (all(!is.na(parsed_dates))) {
      validated_format <- fmt
      valid_format_found <- TRUE
      break
    }
  }
  
  if (valid_format_found) {
    return(validated_format)
  } else {
    stop(paste("No match found for datetime format. Must be one of ", 
               paste(datetime_formats, collapse=", ")))
    return(NA)
  }
}



#' This function aligns the internal variable and attributes with the names used
#' in the input CSV data. This will stop if no match is found. This function will
#' by default automatically search through the attributes in the input data frame
#' to find attributes that correspond to the latitude, longitude, datetime, and 
#' depth. The names of these will be saved in the cache metadata. 
#' This contains parameters to override the automatic detection of names. These
#' are NULL by default that indicates that the automatic name lookup should be 
#' performed.
#' @param data_cache list: with df, df_metadata and ds 
.identify_column_names_and_format <- function(data_cache, 
            timezone_input_data_src, timezone_download_data_src, 
            verbosity) {
  
  variables_to_download <- data_cache$df_metadata$var_names

  df <- data_cache$df
  
  # --------- Match up column names with input data ------------
  # Find the column name used in the input from the possible list of
  # expected column names.
  find_column <- function(column_names, df) {
    for (name in column_names) {
      if (name %in% colnames(df)) {
        return(name)
      }
    }
    return(NULL)
  }
  

  data_cache$df_metadata$input_latitude_name <- find_column(netcdf_points_downloader_globals$input_lat_names, df)
  data_cache$df_metadata$input_longitude_name <- find_column(netcdf_points_downloader_globals$input_lon_names, df)
  data_cache$df_metadata$input_datetime_name <- find_column(netcdf_points_downloader_globals$input_datetime_names, df)
  data_cache$df_metadata$input_depth_name <- find_column(netcdf_points_downloader_globals$input_depth_names, df)
  
  # --------- Determine date time format and calculate download URLs ---------
  # Identify datetime format of the input file by iteratively trying to convert
  # the input date with a specific datetime format. If it fails then it is not
  # a match.

  if (!is.null(data_cache$df_metadata$input_datetime_name)) {
    datetime_samples = df[[data_cache$df_metadata$input_datetime_name]]

    formats_to_check <- netcdf_points_downloader_globals$input_datetime_formats
    data_cache$df_metadata$input_datetime_format <- identify_datetime_format(datetime_samples, formats_to_check)
 
    # Convert the input string dates into date objects with the correct
    # time zones set.
    # We need to convert the input date times to match the timezone of the
    # download data source so that the path names will work.
    # Save the converted form so we don't have to keep doing conversions.
    # POSIXct saves the datetime as seconds from 1 Jan 1970 UTC and has an
    # attribute tzone that remembers the time zone that was used in the 
    # conversion. It will use this for converting from the UTC time when 
    # displaying the date.

    data_cache$df$input_datetime_posixct <- as.POSIXct(datetime_samples, 
        format = data_cache$df_metadata$input_datetime_format, tz = timezone_input_data_src)
    data_cache$df$download_datetime_posixct <- data_cache$df$input_datetime_posixct
    # Make the download timezone dates display in its timezone
    attr(data_cache$df$download_datetime_posixct, "tzone") <- timezone_download_data_src
 
 
    # Calculate the URL paths for each of the points to be downloaded. For
    # datasets that have an aggregate data source this will all be the same
    # path, but for data that is split across multiple files or URLs then
    # this will vary.
    for (row_index in seq_len(nrow(df))) {
        data_cache$df$download_data_src_url[row_index] <- generate_data_source_url(
              data_cache$df$download_datetime_posixct[row_index], 
              data_cache$df_metadata$download_data_src_url_template)
    }
    if (verbosity > 1) {
      cat("Precalculated URL end points. For row 1: ", data_cache$df$download_data_src_url[1], "\n")
    }
  } else {
    # If there is no time dimension then the download_data_src_url for all requests will be 
    # the same as the url template.
    data_cache$df$download_data_src_url <- data_cache$df_metadata$download_data_src_url_template
  }
  
  # ------ Check the variables available for download and the dimensions of variables -----
  #ds <- open_ds(initial_input_data_src_url, verbosity)
  #data_cache$ds <- ds
  # Ensure we have a connection to the data source. Since we just need to get
  # metadata for the time series we can get this from the first row of data
  result <- .open_ds(data_cache, 1,  verbosity)
  data_cache <- result$data_cache
  
  # Verify variables to download exist in the download dataset
  vars_available <- names(data_cache$ds$var)
  missing_vars <- setdiff(variables_to_download, vars_available)
  if (length(missing_vars) > 0) {
    stop(paste("The following variables are not available in dataset:", 
      paste(missing_vars, collapse = ", "), 
      " available vars: ",paste(vars_available, collapse = ", ")))
  }
  

  
  # @param variable name of the variable to download
  # @param match_dim_names_over_vars - matching dimension names over the variables
  #       processed to date. They should all be the same if the data source is
  #       consistent.
  # @param dim_name - name of the dimension (download data source names)
  stop_if_not_unique <- function(variable, match_dim_names_over_vars, dim_name) {
    
    if (length(match_dim_names_over_vars) == 0) {
      stop("Dimension name in download source not recognised: variable: ", variable, 
      ", dim_name: ", dim_name, 
      ". See netcdf_points_downloader_globals for recognised dimension names.")
    }
    # If we find multiple different names in the download data source for
    # logically the same dimension across multiple variables, i.e. something 
    # like: 
    # temp using dimensions: lat, long, depth, time
    # salt using dimensions: latitude, long, depth, time
    # Here the name of the latitude dimension is not consistent. This would 
    # also mean there are multiple latitude dimensions in the data source.
    # This script is not made to handle this case so stop.
    if (length(unique(match_dim_names_over_vars)) != 1) {
      stop("Inconsistent variable dimensions. This script doesn't support this. var: ", 
        variable, " dim set: ", paste(unique(dim_names), collapse=", "))
    }
  }
    
  # Make sure that if we have variables that have time or depth dimensions
  # that we have an input datetime and/or depth column. If not print an error.
  # Also record if a variable needs time or depth dimensions so that we
  # know how to correctly extract the data later.
  
  # Check if time or depth dimensions are required for any of the variables to be downloaded
  time_vars <- c()
  depth_vars <- c()
  
  dim_names_found <- c()    # Record for error reporting
  latitude_names <- c()     # Collect all the latitude names
  longitude_names <- c()
  datetime_names <- c()
  depth_names <- c()
    
  for (var in variables_to_download) {
    
    # Initialize an empty list for this variable if it doesn't exist yet
    if (!var %in% names(data_cache$df_metadata$variable_attributes)) {
      data_cache$df_metadata$variable_attributes[[var]] <- list(needs_time = FALSE, needs_depth = FALSE, dim_index_names_vector = NA)
    }
    
    # Collect all the dimension names for the variable
    dim_names <- sapply(data_cache$ds$var[[var]]$dim, function(d) d$name)
    
    
    # Create a list contain the mapping between the dimension names in the download
    # source and standard internal dimension names.
    # In this case we need to map from the name specified in the download service
    # to its meaning. We also want to retain the ordering of the dimensions as
    # this is needed for the start_indices for downloading the data.
    # Here we assume that dimension names are consistent over the dataset, i.e.
    # if one variable uses dims 'lat' and 'lon' then another doesn't use 'latitude' 
    # and 'longitude'. While this is technically possible the rest of the application
    # will not cope with this as we use this assumption to save on repeated 
    # calculations in the download. Instead we go for detecting this edge case
    # and stopping if it occurs.
    # 
    # Here we add an exception for eReefs because the dimension names don't match
    # the variable names needed to get the information about that variable.
    # For example: depth has a dimension called k, but a matching variable called zc
    # For curvilinear data the dimensions are i and j that align with variables
    # latitude and longitude. I couldn't think of a way of automatically determining
    # these associations and so they are hard coded.
    dim_index_names_vector <- c() # translated index names for start_indices 

    for (i in seq_along(dim_names)) {
      name <- NULL
      dim_name <- dim_names[i]
      
      # If we find an eReefs dimensions for a variable we assume we have eReefs
      # data and we swap in the name of the variable that is needed for that
      # dimension. For example if we see a dimension 'k' then we know the
      # variable needed is 'zc'. Swap out the name so it will splice into
      # the next section and be filed as though the dimension names matched
      # the variable names.
      ereefs_dim <- netcdf_points_downloader_globals$download_dim_map[[dim_names[i]]]
      if (!is.null(ereefs_dim)) {
        dim_name <- ereefs_dim
      }
      
      # Using latitude as an example:
      # For each dimension associated with a variable work out if it is associated
      # with the known 'latitude' names, if so assume it is a latitude dimension
      # and record the name of the index variable ('lat_idx') we need to use
      # to send download request (as part of start_indices). Also record the name of the
      # variable in the download source that is needed for accessing latitude
      # information. For example a download data source with a 'lat' dimension will
      # match in the $download_lat_names, we therefore know that this is a latitude
      # dimension and that its name is 'lat'. We save it in the list of latitude names
      # found in the variables (latitude_names). In theory the latitude name should
      # be consistent across all variables, unless we have a super complex NetCDF
      # source that we can't handle. By recording all the latitude names we can detect
      # if this is the case. After each addition we check if the newly added name
      # doesn't match the previously seen names. 
      
      if (dim_name %in% netcdf_points_downloader_globals$download_lat_names) {
        name <- 'lat_idx'
        latitude_names <- c(latitude_names, dim_name)
        stop_if_not_unique(var, latitude_names, dim_name)
      } else if (dim_name %in% netcdf_points_downloader_globals$download_lon_names) {
        name <- 'lon_idx'
        longitude_names <- c(longitude_names, dim_name)
        stop_if_not_unique(var, longitude_names, dim_name)
      } else if (dim_name %in% netcdf_points_downloader_globals$download_datetime_names) {
        name <- 'time_idx'
        datetime_names <- c(datetime_names, dim_name)
        stop_if_not_unique(var, datetime_names, dim_name)
        
        data_cache$df_metadata$variable_attributes[[var]]$needs_time <- TRUE
        time_vars <- c(time_vars, var)
      } else if (dim_name %in% netcdf_points_downloader_globals$download_depth_names) {
        data_cache$df_metadata$variable_attributes[[var]]$needs_depth <- TRUE
        depth_vars <- c(depth_vars, var) 
        name <- 'depth_idx'
        depth_names <- c(depth_names, dim_name)
        stop_if_not_unique(var, depth_names, dim_name)
      }
      
      # NULL means we didn't find a matching dimension
      if (!is.null(name)) {
        dim_index_names_vector <- c(dim_index_names_vector, name)
        dim_names_found <- c(dim_names_found, dim_name)
      }
    }
    
    # Save for later so we can calculate the start_indices when fetching data
    data_cache$df_metadata$variable_attributes[[var]]$dim_index_names_vector <- dim_index_names_vector
    
    # dim_index_names_vector contains the matching internal index names, such as c('lat_idx','lon_idx','time_idx')
    # If this is not the same length as the dimensions for this variable, say ('lat','lon','datetime', 'weird-depth-name'), 
    # then we missed a match.
    if (length(dim_names) > length(dim_index_names_vector)) {
      stop("Input dimension name that is not recognised: ", setdiff(dim_names, dim_names_found), ". See netcdf_points_downloader_globals for recognised dimension names.")
    }
  }
  
  # Save the name of the variable needed to be downloaded to get information about
  # the dimension.
  # Will be NULL if unused dimension. We use this to indicate that the dimension is not needed.
  # We can pick the first in vector because we will only
  # get here if all values in each of the latitude_names are the same because of stop_if_not_unique()
  data_cache$df_metadata$download_latitude_name <- latitude_names[1]  
  data_cache$df_metadata$download_longitude_name <- longitude_names[1]
  data_cache$df_metadata$download_datetime_name <- datetime_names[1]
  data_cache$df_metadata$download_depth_name <- depth_names[1]

  
  
  # Make sure that if time / depth is needed that it is specified in the input data
  errors <- character()
  if (length(time_vars) > 0 && is.null(data_cache$df_metadata$input_datetime_name)) {
    time_vars_joined <- paste(time_vars, collapse = ', ')
    errors <- c(errors, paste("Variables (", time_vars_joined, ") need a time dimension, but no time column found in input CSV.", sep = ""))
  }
  if (length(depth_vars) > 0 && is.null(data_cache$df_metadata$input_depth_name)) {
    depth_vars_joined <- paste(depth_vars, collapse = ', ')
    errors <- c(errors, paste("Variables (", depth_vars_joined, ") need a depth dimension, but no depth column found in input CSV.", sep = ""))
  }
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"))
  }
  
  # Make sure that the input CSV has columns for latitute and longitude
  errors <- character()
  if (is.null(data_cache$df_metadata$input_latitude_name)) {
    errors <- c(errors, "latitude column")
  }
  
  if (is.null(data_cache$df_metadata$input_longitude_name)) {
    errors <- c(errors, "longitude column")
  }
  
  if (length(errors) > 0) {
      stop("Match for ", paste(errors, collapse = ", "), " in input CSV not found.")
  }
  
  return(data_cache)
}

# Function to estimate distance threshold for curvilinear grids
# This is used to flag points that are out of bounds. Since we match
# the specified locations with the grid array by looking at the closest
# point in the grid it is ambiguous whether the point it inside the grid.
# Since we don't want to hard code anything in, i.e. baking in the
# distance for GBR1 and GBR4, we estimate the grid spacing. Here we assume
# that the grid is not overly distorted. Here we are calculate the 
# grid spacing assuming that the grid is not stretched, i.e. it is a regular
# grid. We then set the threshold to be 3 times this limit. This effectively
# allows grid warping of cells. If the threshold is too small then valid
# points will be falsely triggered as outside of the grid. If the threshold
# is made too large then points near the edge of the grid will be falsely
# indicated as valid. Luckily the eReefs grids has empty cells around the edge.
.estimate_distance_threshold <- function(lat_array, long_array) {
  # Calculate lat and long range
  lat_range <- max(lat_array, na.rm = TRUE) - min(lat_array, na.rm = TRUE)
  long_range <- max(long_array, na.rm = TRUE) - min(long_array, na.rm = TRUE)
  
  # Determine the number of grid points (rows for lat, cols for long)
  num_rows <- dim(lat_array)[1]
  num_cols <- dim(lat_array)[2] # Assuming lat_array and long_array have the same dimensions
  
  # Estimate average grid spacing
  avg_lat_spacing <- lat_range / num_rows
  avg_long_spacing <- long_range / num_cols
  
  # Use the smaller of the two spacings as a conservative estimate of grid spacing
  grid_spacing <- min(avg_lat_spacing, avg_long_spacing)
  
      
                          
  distance_threshold <- grid_spacing * netcdf_points_downloader_globals$curvilinear_grid_squish
  return(distance_threshold)
}

.is_regular_grid <- function(lat_array, long_array) {
  return (length(dim(lat_array)) == 1 && length(dim(long_array)) == 1)
}

.is_curvilinear_grid <- function(lat_array, long_array) {
  return (is.matrix(lat_array) && is.matrix(long_array))
}


  # Function to check if any variable needs the time dimension
  .any_variables_need_time <- function(df_metadata) {
    # Iterate through each variable in variable_attributes
    for (var_name in names(df_metadata$variable_attributes)) {
      # Check if this variable needs the time dimension
      if (df_metadata$variable_attributes[[var_name]]$needs_time) {
        # As soon as we find one, we can return TRUE
        return(TRUE)
      }
    }
    return(FALSE)
  }

  # Function to check if any variable needs the depth dimension
  .any_variables_need_depth <- function(df_metadata) {
    # Iterate through each variable in variable_attributes
    for (var_name in names(df_metadata$variable_attributes)) {
      # Check if this variable needs the time dimension
      if (df_metadata$variable_attributes[[var_name]]$needs_depth) {
        # As soon as we find one, we can return TRUE
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  # Define a function to find the matching dimension and variable names in the data source
  # This is needed to match up to the latitude and longitudes. We look through both
  # the dimensions and the variables. For regular eReefs grid 'latitude' is a dimension
  # where as for a curvilinear eReefs the dimension is 'i' and the 'latitude' is a
  # variable. This search is also not case sensitive.
  .find_dim_name_deprecated <- function(ds, possible_names) {
    dim_names <- names(ds$dim)
    var_names <- names(ds$var)
    
    # Combine the dimension names and variable names into one vector
    combined_names <- c(dim_names, var_names)
    
    # Normalize the names for comparison
    combined_names_lower <- tolower(combined_names)
    
    # Check for duplicates after normalisation
    if(length(unique(combined_names_lower)) != length(combined_names)) {
      stop("Error: Duplicate dimension or variable names found after combining.")
    }

    # Search for the target dimension or variable name
    name_match <- NULL
    for (name in possible_names) {
      if (tolower(name) %in% combined_names_lower) {
        name_match <- combined_names[combined_names_lower == tolower(name)]
        break
      }
    }
    
    # If not found, throw an error
    #if (is.null(name_match) || length(name_match) == 0) {
    #  stop(paste("Error: Target dimension or variable names", 
    #    paste(possible_names, collapse = ", "), "not found in the dataset.\n"))
    #}
    
    return(name_match)
  }
  

# Gets the base date from the time 'units' attribute. This is needed to calculate
# the times from the floating point time values.
.get_base_date <- function (ds, time_name, timezone){
    # Get the time metadata. The data is typically provided as a floating point number
    # that is a certain number of dates since a base date. We look for the base date in
    # the metadata 'units' of the download data source.
    # Examples: 
    # https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_2024-01-17.nc.html
    # 'days since 1990-01-01 00:00:00 +10' (eReefs)
    # https://dapds00.nci.org.au/thredds/dodsC/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/sfcWind/v20231001/sfcWind_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_202308-202308.nc.html
    # 'days since 1949-12-01' (BARRA2)
    #
    # IMOS Uses seconds as units
    # https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3SM-1d/day/2020/20200101032000-ABOM-L3S_GHRSST-SSTskin-MultiSensor-1d_day.nc.html
    # 'seconds since 1981-01-01 00:00:00'
    time_units <- ncatt_get(ds, time_name, "units")
    
    if (!time_units$hasatt) {
      stop("Input data source has no time units metadata")
    }

    split_units <- strsplit(time_units$value, " since ", fixed = TRUE)[[1]]
    time_unit <- split_units[1]  # "days" or "seconds"
    base_date <- split_units[2]  # "1990-01-01 00:00:00 +10"
    
    if (time_unit == "days") {
      seconds_per_unit <- 24*3600
    } else if (time_unit == "seconds") {
      seconds_per_unit <- 1
    } else {
      stop("This script only handles time units in 'seconds' and 'days' but in this data source it is: ", time_unit)
    }
    
    
    # Convert the base date to POSIXct
    base_date_posix <- as.POSIXct(base_date, tz = timezone)
    return(list(base_date_posix = base_date_posix, seconds_per_unit = seconds_per_unit))
}


#' This function goes through all the sample locations to be downloaded and 
#' calculates their indicies in the eReefs grid. This is needed for the 
#' data request. We precalculate this because if there are points that are out
#' bounds we can flag these earlier in the process. Additionally the indicies
#' can be reused when downloading from multiple variables. As indicies are
#' calculated we save them to the disk cache.
#' @param data_cache list: this function will add columns to data_cache$df
#'        for each of the dimensions relevant for the download.
#' @param timezone_download_data_src: Timezone of the date times specified in the input data source 
#'        (points to download the data). For example:
#'        'Etc/GMT-10' refers to Queensland time. While this seems counter 
#'        intuitive Etc/GMT-X refers to a time zone which is X hours ahead of UTC.
#'        This is only used if it is not specified in the data times themselves.
#' @param verbosity integer: amount of messages to print, 0 - None, 2 - all.
#' @return data_cache
.precalculate_indices <- function(data_cache, timezone_download_data_src, verbosity) {
  df <- data_cache$df
  df_metadata <- data_cache$df_metadata
  
  
  # Check if all indices have already been calculated
  if(all(df$index_complete, na.rm = TRUE)) {
    if (verbosity > 1) cat("All indices have already been calculated.\n")
    return(data_cache)
  }
  
  
  # If we have previously done some processing (a restart) indicate where we are 
  # resuming from. 
  total_rows <- nrow(df)
  already_processed <- sum(df$index_complete, na.rm = TRUE)
  if ((already_processed > 0) && (verbosity > 1)) {
    cat(sprintf("Resuming index calculation from row %d of %d.\n", already_processed, total_rows))
  } 

  start_time <- Sys.time()
  

  
  ds <- data_cache$ds
  if (is.null(ds)) {
    stop("data_cache$ds should already been connected")
  }
  
  # ---- Match the download data source dimension names to internal names ----
  # If a dimension is unused then its name will be NULL

  any_time_needed <- .any_variables_need_time(df_metadata)
  any_depth_needed <- .any_variables_need_depth(df_metadata)

  if (any_depth_needed) {
    #data_cache$df_metadata$download_depth_name <- .find_dim_name(ds, netcdf_points_downloader_globals$download_depth_names)
    depth_array <- ncvar_get(ds, data_cache$df_metadata$download_depth_name) 
  }

  # The goal is to get the arrays for the latitude and longitude dimensions.
  # For AIMS eReefs regular grid arrays, latitude and longitude are dimension
  # values. For CSIRO curvilinear grid the latitude and longitude are
  # not dimensions, but variables. The dimensions are i and j which correspond
  # to indices into the latitude and longitude variables.
  # For other data sources the names of the latitude and longitude are sometimes
  # different. For example IMAS 
  # https://thredds.imas.utas.edu.au/thredds/dodsC/IMAS/SeamapAus_Bathymetry_National_50m/bathy_01_flt_clp.nc
  # It is a regular grid, with dimensions LATITUDE and LONGITUDE. Note all caps.
  
  # If the data source is a regular grid then the latitude and longitude is
  # specified in the dimensions. If it is curvilinear then the latitude and
  # longitude are in the variables.
  # We also know that the names are not always consistent across data sources.

  # Possible names for latitude and longitude dimensions/variables
  # This is case insenstive so "latitude" will match "LATITUDE"

  lat_name <- data_cache$df_metadata$download_latitude_name
  lon_name <- data_cache$df_metadata$download_longitude_name
  
  if (verbosity > 1) cat("Data source latitude name: ", lat_name, " longitude name: ", lon_name, "\n")
  
  # Now use the actual names to download the grid details
  # lat_array Numeric array: 1D or 2D array of latitudes representing 
  #      the grid. For curvilinear grids this will be a 2D array, corresponding
  #      to the locations of each of the slightly different latitudes of all the 
  #      points in the grid. For a regular grid this corresponds to the linear
  #      array of points corresponding to the latitudes of the rows of the grid.
  # long_array Numeric array: 1D or 2D array of longitudes describing the
  #      grid.
  lat_array <- ncvar_get(ds, lat_name)
  long_array <- ncvar_get(ds, lon_name)

  if (!is.array(lat_array) | !is.array(long_array)) {
    stop("lat_array and long_array must both be either array (1D) or matrices (2D)")
  }
  
  # Sanity check incase there is a grid type we haven't considered and set up
  # variables for the loop.
  if (.is_regular_grid(lat_array,long_array)) {
    regular_grid <- TRUE
    distance_threshold <- NA
    min_lat_array <- min(lat_array, na.rm = TRUE)
    max_lat_array <- max(lat_array, na.rm = TRUE)
    min_long_array <- min(long_array, na.rm = TRUE)
    max_long_array <- max(long_array, na.rm = TRUE)
    if (verbosity > 1) cat("Regular grid bounds:",
      " Latitude: ", min_lat_array, " - ", max_lat_array,
      " Longitude: ", min_long_array, " - ", max_long_array, "\n")
  } else if (.is_curvilinear_grid(lat_array,long_array)) {
    regular_grid <- FALSE
    # distance_threshold is only relevant for curvilinear grid.
    # If the closest grid cell to the point is further than the distance
    # then we consider the point out of bounds.
    distance_threshold <- .estimate_distance_threshold(lat_array, long_array)
  } else {
    stop("Unknown grid type. Doesn't match tests for regular grid and curvilinear grids\n")
  }
  
  if (verbosity > 1) cat("Calculating indices for input locations\n")
  
  first_time <- TRUE
  
  
  # ----------- Processing per row -------------------
  # Find the matching indices for each of the input locations. Saves these
  # in lat_idx, lon_idx, time_idx and depth_idx columns.
  for (row_index in seq_len(nrow(df))) {
    
    # Calculate indices for row if not already calculated
    if (!df$index_complete[row_index]) {
      
      # Look up the latitude and longitude of the row that we are
      # downloading. We use the column names in the data file that
      # we determined earlier and saved in df_metadata.
      input_lat <- df[[df_metadata$input_latitude_name]][row_index]
      input_lon <- df[[df_metadata$input_longitude_name]][row_index]

      if (is.na(input_lat)) {
        stop(sprintf("Input row %d, latitude is NA. All rows must have a valid latitude.", row_index))
      }
      if (is.na(input_lon)) {
        stop(sprintf("Input row %d, longitude is NA. All rows must have a valid longitude.", row_index))
      }
      # Find the latitude and longitude indices. This finds the nearest point in the 
      # lat and long of the cell with the smallest distance squared. This approach 
      # is inefficient as it needs to calculate the distance to every location in the 
      # lat and long arrays, but it is simple and robust. The problem is with points
      # that are outside the grid. This method will still find the closest in grid
      # cell, which is not what we want. For regular grid arrays we check that the
      # point is inside the bounds of the grid. For a curvilinear grid we check that
      # the distance between the point and the closest grid is less than the
      # distance threshold, which we estimate previously outside the loop.      

      lat_idx <- NA
      long_idx <- NA
      lat_val <- NA
      long_val <- NA

      if (regular_grid) {
        # Regular grid
        lat_idx <- which.min(abs(lat_array - input_lat))
        long_idx <- which.min(abs(long_array - input_lon))
        # Check that the point is within the bounds of the array.
        
        if (input_lat >= min_lat_array && input_lat <= max_lat_array && 
          input_lon >= min_long_array && input_lon <= max_long_array) {
          lat_val <- lat_array[lat_idx]
          long_val <- long_array[long_idx]
        } else {
          msg <- sprintf(
            "Requested location (input_lat: %f, input_lon: %f) is outside the grid bounds (input_lat: %f to %f, input_lon: %f to %f).", 
            input_lat, input_lon, 
            min_lat_array, 
            max_lat_array, 
            min_long_array, 
            max_long_array)
          warning(msg)
          # Indicate that the download should be skipped
          df$vars_download[row_index] <- TRUE
        }
      } else {
        # Curvilinear grid
        difference_array <- (lat_array - input_lat)^2 + (long_array - input_lon)^2
        min_distance <- sqrt(min(difference_array, na.rm = TRUE))  # Ignore NaN values
        
        if (verbosity > 1) cat(sprintf("min_distance: %f, distance_threshold: %f\n", min_distance, distance_threshold))
        if (min_distance <= distance_threshold) {
          index_flat <- which.min(difference_array)
          index <- arrayInd(index_flat, dim(difference_array))
          lat_idx <- index[1]
          long_idx <- index[2]
          lat_val <- lat_array[lat_idx, long_idx]
          long_val <- long_array[lat_idx, long_idx]
        } else {
          warning(sprintf(
            "Location is outside the acceptable distance (input_lat: %f, input_lon: %f, dist: %f). Nearest distance: %f.\n", 
            input_lat, input_lon, distance_threshold, min_distance))
          # Indicate that the download should be skipped
          df$vars_download[row_index] <- TRUE
        }
      } 

      # Save the indices 
      df$lat_idx[row_index] <- lat_idx
      df$lon_idx[row_index] <- long_idx
      df$lat_val[row_index] <- lat_val
      df$lon_val[row_index] <- long_val
      

      # Create place holders, so later printout is consistent
      time_idx <-  NA
      time_val <- NA


      # Find time index
      if (any_time_needed) {

        # If the time series is made up from multiple end points, such as monthly files
        # with daily data, and the monthly files are at different URL end points in a
        # series.
        # As we iterate through the input data to download we will need to switch 
        # data sources. This will retain the existing data source connection if
        # the end point doesn't change.

        result <- .open_ds(data_cache, row_index, verbosity)
        data_cache <- result$data_cache
        
        # Only perform this calculation each time the data source changes, such
        # as switching from file to file. This is so that we don't perform unnecessary
        # ncvar_get requests for each time step. For most datasets the base_date is constant
        # throughout the dataset, but for some (BARRA2) the base date is different for some files.
        if (result$download_data_src_changed || first_time) {
          first_time <- FALSE
          time_name <- data_cache$df_metadata$download_datetime_name
          # Get the time values for the current data source end point. 
          time_vals <- ncvar_get(ds, time_name)
          
          base_date_list <- .get_base_date(ds, time_name, timezone_download_data_src)
          
          # Convert the time values to seconds (from days) and add to the base POSIXct date
          # Note: as.POSIXct numeric origin is "1970-01-01", but here we use direct arithmetic with base_date_posix
          # time_array now holds POSIXct objects representing each time value in ds$dim$time$vals
          time_array <- base_date_list$base_date_posix + (time_vals * base_date_list$seconds_per_unit)  # 86400 seconds in a day
          
          # Do calculations that only need updating each time time_array is updated
          
          # Work out the acceptable time range for observations.
          # There is a complication here. If we are extracting from a set of end points
          # that are time aggregations then often the time range in each file is only
          # a point in time. For example AIMS eReefs GBR1 Daily (individual files) the time range for
          # an individual end point is: 2020-01-03 00:00:00 - 2020-01-03 00:00:00.
          # This is also the case with GHRSST which correspond to daily data in a daily
          # file, with only a single time slice. For the day time temperatures the time corresponds
          # to 03:20 UTC.
          # We therefore use the following process:
          # 1. If there is only one time slice in the time_array then allow a default
          #    time tolerance of 1 day.
          # 2. If the time_array has multiple times, we work out the average gap between 
          #    time slices and allow through observations that are one average time gap
          #    before the start and after the end of the time_array.
          # 3. Anything else flag as out of bounds.
          max_time <- max(time_array)
          min_time <- min(time_array)
          if (length(time_array) == 1) {
            # If we only have a single time slice we can't work out the spacing in 
            # time samples and thus what an acceptable amount of tolerance to the 
            # time matching before generating a warning, so use a default value.
            time_tolerance <- netcdf_points_downloader_globals$single_end_point_time_tolerance_sec
          } else {
            # Time series in file or end point. Estimate spacing between time slices.
            time_tolerance <- (max_time - min_time)/(length(time_array) - 1)
          }
          time_range_min <- min_time - time_tolerance
          time_range_max <- max_time + time_tolerance
        }
    
#        if (verbosity > 1) {
#          cat(sprintf("Input obs time: %s, tzone: %s; adj to tzone %s: %s; data src tzone: %s\n", 
#            df$input_datetime_posixct, attr(observation_time, "tzone"), 
#            timezone, datetime_character, 
#            attr(min(time_array), "tzone")))
#        }
        observation_time <- df$download_datetime_posixct[row_index]
        
        
        
        
        if (observation_time < time_range_min || observation_time > time_range_max) {
          warning(sprintf("Time (%s) for row %d is out of range (%s - %s). Download was skipped.", 
            format(min(observation_time), df_metadata$input_datetime_format), 
            row_index, 
            format(min(time_array), df_metadata$input_datetime_format), 
            format(max(time_array), df_metadata$input_datetime_format)))
          # Indicate that the download should be skipped
          df$vars_download[row_index] <- TRUE
          df$time_idx[row_index] <- NA
          df$time_val[row_index] <- NA
        } else {
          # Find the closest time index
          time_idx <- which.min(abs(time_array - observation_time))
          df$time_idx[row_index] <- time_idx
    
          # Format the POSIXct time value to a string with the specified format
          time_val <- format(time_array[time_idx], df_metadata$input_datetime_format)
          df$time_val[row_index] <- time_val
        }
      }
      
      depth_idx <- NA
      depth_val <- NA
      # Find depth index
      if (any_depth_needed) {
        observation_depth <- df[[df_metadata$input_depth_name]][row_index]
        # Check if the observation depth is within the range of the depth_array
        if (observation_depth < min(depth_array) || observation_depth > max(depth_array)) {
          warning(sprintf("Depth (%f) for row %d is out of range (%f - %f). Download was skipped.", 
            observation_depth,
            row_index,
            min(depth_array),
            max(depth_array)))
          df$vars_download[row_index] <- TRUE
          df$depth_idx[row_index] <- NA
          df$depth_val[row_index] <- NA
        } else {
          observation_depth <- df[[df_metadata$input_depth_name]][row_index]
          # Find the closest depth index.
          depth_idx <- which.min(abs(depth_array - observation_depth))
          df$depth_idx[row_index] <- depth_idx
          depth_val <- depth_array[depth_idx]
          df$depth_val[row_index] <- depth_val
        }
      }

      if (verbosity > 1) {
        cat(sprintf("row %d: lat, long, time, depth: indices: %d, %d, %d, %d, values: %f, %f, %s, %f\n",
          row_index, lat_idx, long_idx, time_idx, depth_idx, lat_val, long_val, time_val, depth_val))
      }
      
      # Record that we completed the index calculations for this row.
      df$index_complete[row_index] <- TRUE
      
      # Periodically save the results.
      processed_count <- row_index
      if (processed_count %% 2000 == 0 || processed_count == nrow(df)) {
        .save_data_cache(data_cache)
        if (verbosity > 0) cat(sprintf("Calculated indicies for %d of %d rows. Saving...\n", processed_count, total_rows))
      }
    }
  }

  elapsed_time <- difftime(Sys.time(), start_time, units = "secs")
  if (verbosity > 1) cat(sprintf("Elapsed time: %.2f sec. Processing look up table index complete.\n", elapsed_time))
  
  # While this should be unnecessary because df should be a reference in data_cache, this isn't
  # the case here for some reason. For this reason we assign to ensure the modifications are saved.
  data_cache$df <- df
  return(data_cache)
}

# This maps the input indices to the order needed for making the data request. It 
# maps it to the order specified in the dimensions of the variable. We previously
# mapped for each variable the types and order into a dim_index_names_vector.
# For example for dim_index_names_vector: c("longitude", "latitude", "time")
# the dim_index_names_vector would be c("lon_idx", "lat_idx", "time_idx").
# This substitutes in the indices provided. For this example with an input:
# lat_idx = 12, lon_idx = 3, depth_idx = 4, time_idx = 298
# the output would be:
# c(3, 12, 298)
.get_start_indices <- function(data_cache, variable, lat_idx, lon_idx, depth_idx, time_idx, verbosity) {
  # Retrieve the dimension index names vector for the specified variable
  dim_index_names_vector <- data_cache$df_metadata$variable_attributes[[variable]]$dim_index_names_vector
  
  # Check if dim_index_names_vector is properly setup
  if (is.null(dim_index_names_vector)) {
    stop("BUG: dim_index_names_vector is null")
  }
  
  # Initialize an empty vector to hold the start indices
  start_indices <- numeric(length(dim_index_names_vector))
  
  # Iterate over dim_index_names_vector to map each dimension name to its corresponding index
  for (i in seq_along(dim_index_names_vector)) {
    dim_name <- dim_index_names_vector[i]
    if (dim_name == "lat_idx") {
      start_indices[i] <- lat_idx
    } else if (dim_name == "lon_idx") {
      start_indices[i] <- lon_idx
    } else if (dim_name == "depth_idx") {
      start_indices[i] <- depth_idx
    } else if (dim_name == "time_idx") {
      start_indices[i] <- time_idx
    } else {
      stop(paste("BUG: Unrecognized dimension name. This should have been caught earlier. '", 
        dim_name, "' for variable '", variable, "'.", sep = ""))
    }
  }
  
  # Return the start indices in the order specified by dim_index_names_vector
  return(start_indices)
}








#' This function is intended to be run after the indicies have been calculated and
#' saved in the df. This then starts download the eReefs variables, one row and
#' one variable at a time.
#' This assumes that points that are out of bounds have been dealt with earlier.
#' @param data_cache list: 
#' @param verbosity integer: Amount of print messaging 0 - 2, 0 - no messages, 2 - most verbose. 
#' @param download_rate_limit_delay - numeric: Number of seconds pause between requests
#'          sent to the server. Setting this to 0 will make the download as fast as 
#'          possible, at the expense of potentially heavily loading the server that the
#'          data is being loaded from. Having a small detail ensures that we don't 
#'          hog the server, or overload it.
#' @param stop_on_error boolean: If TRUE then an error during the download will stop
#'          the application. If FALSE then a message will be reported and the variable
#'          saved as NA. Default is FALSE. Setting this to TRUE is useful for 
#'          running tests.
.download_data_with_precalculated_indices <- function(data_cache, verbosity,
    download_rate_limit_delay = 0.1, stop_on_error = FALSE) {

  df <- data_cache$df
  df_metadata <- data_cache$df_metadata

  variable_needs_time <- function(df_metadata, variable) {
    return (df_metadata$variable_attributes[[variable]]$needs_time) 
  }

  variable_needs_depth <- function(df_metadata, variable) {
    return (df_metadata$variable_attributes[[variable]]$needs_depth) 
  }

  # Number of points to download between reporting on progress
  report_rate <- 25
  
  # Track the start of the download so we can estimate how long it will take to finish
  start_time <- Sys.time()
  
  # Variables to download
  var_names <- df_metadata$var_names
  
  # Number of location / time / depth points to download.
  total_rows <- nrow(data_cache$df)
  
  # Determine the starting point for downloads
  already_downloaded <- sum(data_cache$df$vars_download, na.rm = TRUE)
  
  # Report resuming or starting anew
  if (already_downloaded > 0) {
    if (verbosity > 0) cat("Resuming download at row ", already_downloaded, " of ", total_rows, "\n")
  } else {
    if (verbosity > 1) cat("Starting data download using precalculated indices...\n")
  }
  
  
  # As we progressively save the download, we also save an indication that
  # we have downloaded the variables for a given row. When we resume 
  # need to only download from the rows that are no marked as complete.
  
  # Filter indices of rows that haven't been downloaded yet (vars_download is FALSE)
  indices_to_download <- which(!data_cache$df$vars_download)

  # From names(ds$dim)
  # For AIMS regular grid NetCDF dimensions are: k, latitutde, longitude, time
  # For CSIRO curvilinear grid NetCDF dimensions are: i, j, k, time
  
  # Loop through the row indices of those that need downloading
  for (index in indices_to_download) {
  
    # Make sure we are connected to the right end point for the row. If this
    # row uses the same end point as the previous row then not much will 
    # happen. If it is different then it will switch the data source end point.
    result <- .open_ds(data_cache, index, verbosity)
    data_cache <- result$data_cache
    ds <- data_cache$ds
    
    # Access the row by its original index in df
    row <- data_cache$df[index, ]
  
    lat_idx <- as.integer(row$lat_idx)
    lon_idx <- as.integer(row$lon_idx)
    
    # If there are any variables that require time or depth then 
    # get the precalculated index. Have NULL to indicate that
    # the dimension is no required. This allows us to have a consistent
    # number of parameters to get_start_indices
    time_idx <- NULL
    if (!is.null(df_metadata$input_datetime_name)) {
      time_idx <- as.integer(row$time_idx)
    }
    

    depth_idx <- NULL
    if (!is.null(df_metadata$input_depth_name)) {
      depth_idx <- as.integer(row$depth_idx)
    }
    

    # Download each variable of each row. Save in df in data_cache
    for (variable in var_names) {
      start_indices <- .get_start_indices(data_cache, variable, lat_idx, lon_idx, depth_idx, time_idx, verbosity)
      
      
      tryCatch({
        
        # When there is a network issue the nc_var_get doesn't fail with a useful
        # error message. The message is "C function R_nc4_get_vara_double returned error"
        # which can also occur if the data is out of bounds.
        # It does however print out:
        # Error:curl error: Couldn't connect to server
        # curl error details: 
        # Warning:oc_open: Could not read url
        # Error in Rsx_nc4_get_vara_double: NetCDF: I/O failure
        # Possibly to Stderror.
        
        # Redirect stderr to a text connection. This is to try and capture the error.

        
        # Now use start_indices in ncvar_get
        var_value <- ncvar_get(ds, variable, start = start_indices, count = rep(1, length(start_indices)))

        # In some cases the values in the NetCDF are compacted into a lower fit format
        # using a scale factor and offset. This information is saved in the variable metadata
        # An example of this is the ts variable in BARRA2: 
        # https://dapds00.nci.org.au/thredds/dodsC/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/1hr/ts/v20231001/ts_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_1hr_202002-202002.nc.html
        # In this case the offset changes in every file.
        # The NetCDF library seems to automatically apply this scaling. 
        # I have not done enough testing on multiple files to determine if 
        # scaling is automatically handled robustly. So I am leaving some
        # inital code that was developed to apply the scaling for later reference:
        #print("Getting scale factor")
        #scale_factor <- ncatt_get(ds, variable, "scale_factor")
        #print(scale_factor)
        #if (scale_factor$hasatt) {
        #  print("Getting add_offset")
        #  add_offset <- ncatt_get(ds, variable, "add_offset")
        #  print(add_offset)
        #  if (add_offset$hasatt) {
        #    print(var_value)
        #    adj_value <- var_value * scale_factor$value + add_offset$value
        #    if (verbosity > 1) {
        #      cat("Applying scale factor to ", variable, " raw value: ", var_value, " adjusted: ", adj_value, "\n")
        #    }
        #    var_value <- adj_value
        #  }
        #}
        
        # Update the dataframe with the downloaded value
        data_cache$df[[variable]][index] <- var_value
        
        # Indicate that the download for this row is done.
        data_cache$df$vars_download[index] <- TRUE
      }, error = function(e) {
        
        print(e$message)
        msg <- sprintf("Error downloading %s for row %d, msg: %s", variable, index, e$message)
        if (stop_on_error) {
          stop(msg)
        } else {
          message(msg)
        }
        data_cache$df[[variable]][index] <- NA  # Set to NA in case of error
      })
    }
    
    
    
    Sys.sleep(download_rate_limit_delay) # Pause for download_rate_limit_delay seconds
    
    # Calculate progress and estimate completion time every 'report_rate' rows
    if ((index - already_downloaded) %% report_rate == 0) {
      current_progress <- index
      
      # Estimate completion time
      elapsed_time <- difftime(Sys.time(), start_time, units = "secs")
      start_time <- Sys.time() # Reset start time for the next period
      rows_left <- total_rows - current_progress
      estimated_time_remaining <- as.numeric(elapsed_time) * rows_left / report_rate
      estimated_completion_time <- Sys.time() + as.difftime(estimated_time_remaining, units = "secs")
      
      if (verbosity > 0) {
        if (data_cache$enable_disk_cache) {
          msg <- "Saving cache"
        } else {
          msg <- "Downloading"
        }
        cat(sprintf("%s: %d of %d rows complete. Estimated completion time: %s\n", msg,
                  current_progress, total_rows, format(estimated_completion_time, "%Y-%m-%d %H:%M:%S")))
      }
      # Save progress
      .save_data_cache(data_cache)
    }
  }
  .save_data_cache(data_cache)
  return(data_cache)
}







#' Convert to the final data frame, stripping out the extra fields created to track the download.
.convert_cache_to_final_csv <- function(data_cache, output_csv_path = NULL, verbosity = 1) {

  # Columns to remove
  columns_to_remove <- c('lat_idx', 'lon_idx', 'time_idx', 'depth_idx', 'vars_download', 
    'index_complete', 'attempt_status', 'download_data_src_url', 'input_datetime_posixct',
    'download_datetime_posixct'
)
  
  # Drop specified columns if they exist in the DataFrame
  df_final <- data_cache$df[, !(names(data_cache$df) %in% columns_to_remove)]
  
  # Only save to disk if a path is specified.
  if (!is.null(output_csv_path)) {
    # Ensure the output directory exists
    output_dir <- dirname(output_csv_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    # Save the cleaned DataFrame to the specified CSV file path
    write.csv(df_final, output_csv_path, row.names = FALSE)
  
    if (verbosity > 1) cat(sprintf("Final data saved to %s\n", output_csv_path))
  }
  return(df_final)
}


#' This function downloads NetCDF data for a set of locations, times, and depths
#' specified in a data frame or CSV file (input_data_src) for the variables specified (var_names). 
#' This script will find the closest grid cells in the eReefs data, download the data
#' and save it as an extra column of data.
#'
#' This function was originally developed to assist downloading data from eReefs THREDDS
#' OpenDAP services, but it can also be used for other data sources that have a 
#' similar structure.
#'
#' It starts by calculating the grid cells of data to be downloaded. It then 
#' progressively downloads them one by one. Large downloads of many thousands of points 
#' will be slow (potentially several hours), however this script is robust. 
#' Progress is auto-saved regularly to a RDATA file so that when this function is 
#' restarted with the same inputs then it will read the cache and resume from where 
#' it was previously up to.
#' @param var_names string vector: a vector specifying the short names for 
#'        variables that you want from the netcdf file. Defaults to c('eta').
#' @param input_data_src string or data frame: path to the CSV file that specifies the locations of the
#'        points to download from. Each row should correspond to a single location, time and depth
#'        combination. OR a data frame corresponding to the points to download. Should have
#'        attributes 'lat', 'lon', and optionally 'datetime', 'depth' depending on the variables
#'        being downloaded. It will also accept variations with these names; it will automatically
#'        match names against a preset number of acceptable variations, or they can be
#'        manually specified with the lat_column, lon_column, etc. input parameters.
#'        This input can have extra data columns. They will be propagated to the output.
#'        The generated output will correspond to this input, plus the extra data that is
#'        downloaded from the input_data_src.
#' @param download_data_src_url_template string - URL of the data source to download from. This can be
#'        a link to a file or an OpenDAP URL end point. This can also be a template
#'        corresponding to a set of data end points. This is particularly common when 
#'        there is a set of files making up a time series, or a THREDDS data service
#'        that doesn't have an aggregation data service to combine all the files into
#'        a single time series. If '%(year)', '%(month)' or '%(day)' is included
#'        in the path then the dates from the input_data_src, adjusted for time zone differences,
#'        will be substituted into the URL. All these URLs are calculated and saved in
#'        cache prior to the data download commencing. 
#'        For example: Using GBR1 eReefs without the aggregate data services.
#'        Each individual daily files has an end point data URL like:
#'        https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_2024-01-17.nc
#'        We can use a template such as:
#'        https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_%(year)-%(month)-%(day).nc
#'        to specify how the URL for a given date should be calculated.
#' @param output_csv_path - string: Path to final download. Use this to check if the whole 
#'            download is complete and we don't need to reconnect to the data source.
#'            If set to NULL then no output is saved and the downloaded data is returned
#'            as a data frame from this function. If this is NULL then caching is 
#'            disabled, i.e. enable_disk_cache will be set to FALSE
#' @param timezone_input_data_src - string: Timezone of the date times in the input
#'        data source sample points. The format needs to be suitable for as.POSIXct tz
#'        See https://stat.ethz.ch/R-manual/R-devel/library/base/html/as.POSIXlt.html
#'        
#' @param timezone_download_data_src: Timezone of the opendap input data source. For example:
#'        'Etc/GMT-10' refers to Queensland time. While this seems counter 
#'        intuitive Etc/GMT-X refers to a time zone which is X hours ahead of UTC.
#'        If the strings specify the time zone, such as 1990-01-01 00:00:00 +10
#'        then this time zone is used internally.
#' @param verbosity integer - how much information to display along the way (0 to 2. Default is 1).
#' @param download_rate_limit_delay - numeric: number of seconds pause between requests
#'          sent to the server. Setting this to 0 will make the download as fast as 
#'          possible, at the expense of potentially heavily loading the server that the
#'          data is being loaded from. Having a small delay ensures that we don't 
#'          hog the server, or overload it. (default is 0.1 sec).
#' @param stop_on_error boolean: If TRUE then an error during the download will stop
#'          the application. If FALSE then a message will be reported and the variable
#'          saved as NA. Default is FALSE. Setting this to TRUE is useful for 
#'          running tests. Errors that occur during the startup before the download
#'          ignore this setting and always stop the application.
#' @param enable_disk_cache - boolean: If TRUE then the cache will be saved and read from disk.
#'          If FALSE then there is no effective cache as each run will forget download
#'          progress. This might still be useful if you are downloading a small number of 
#'          points and don't want to clutter the disk with cache files.
#' @return data frame corresponding to the input_data_src with the extra downloaded data.
#'          This includes columns that indicate the grid cell locations that the data was
#'          extracted from time_val, lat_val, lon_val, depth_val. The extracted columns will correspond
#'          to their variable names. The column names for the location and date time are the same
#'          as the input data.
netcdf_points_downloader <- function(var_names = c('eta'),
                                        input_data_src,
                                        download_data_src_url_template,
                                        output_csv_path,
                                        verbosity = 1,
                                        timezone_input_data_src = 'Etc/GMT-10',
                                        timezone_download_data_src = 'Etc/GMT-10',
                                        download_rate_limit_delay = 0.1,
                                        stop_on_error = FALSE,
                                        enable_disk_cache = TRUE){
    # if we are not saving the output then there is no progressive caching
    if (is.null(output_csv_path)) {
      if (enable_disk_cache) {
        stop("enable_disk_cache is TRUE but output_csv_path is NULL. The disk cache is saved based on output_csv_path.\n") 
      }
      enable_disk_cache <- FALSE
      data_cache_path <- NULL
    } else {
      # Handling file paths
      path_without_extension <- file_path_sans_ext(output_csv_path)
      # data_cache_path <- paste0(path_without_extension, "_cache.csv")
      data_cache_path <- paste0(path_without_extension, "_cache.rds")
    }
    
    # Print variables to download
    if (verbosity > 0) cat(sprintf("Variables to download (%s)\n", paste(var_names, collapse = ", ")))

    # Allow the user to provide a data frame specifying the points or
    # a path to a CSV file.
    # Check if input_data_src is a string (indicating a path to a CSV file)
    if (is.character(input_data_src) && length(input_data_src) == 1) {
      # It's a string. Attempt to read the CSV file.
      if (file.exists(input_data_src)) {
        input_df <- read.csv(input_data_src, stringsAsFactors = FALSE)
      } else {
        stop("Provided file path does not exist.", input_data_src)
      }
    } else if (is.data.frame(input_data_src)) {
      # It's a data frame. Use it directly.
      input_df <- input_data_src
    } else {
      stop("input_data_src must be a path to a CSV file or a data frame.")
    }

    # Setup the data_cache list that is used to track the download. We use the same
    # structure regardless of whether the enable_disk_cache is true or not. This is
    # to simplify the code flow. If the enable_disk_cache then this processing is
    # not really a cache because it will forget between successive runs.
    data_cache <- .setup_data_cache(input_df, data_cache_path, var_names, 
      download_data_src_url_template, output_csv_path, 
      timezone_input_data_src, timezone_download_data_src,
      verbosity, enable_disk_cache)

    if (all(data_cache$df$vars_download)) {
      if (verbosity > 0) cat(sprintf(
        "Download process is already complete. Delete cache (%s) to redownload. \n", 
        data_cache_path))
    } else {
      # For debugging:
      #ds_global <<- data_cache$ds
      if (verbosity > 1) cat(sprintf("Note: This script is restartable due to %s. It will resume download on restart.\n", data_cache_path))

      # Check that all the locations are valid
      data_cache <- .precalculate_indices(data_cache, timezone_download_data_src, verbosity)

      data_cache <- .download_data_with_precalculated_indices(
          data_cache, verbosity, download_rate_limit_delay, stop_on_error)
    }
    df_final <- .convert_cache_to_final_csv(data_cache, output_csv_path, verbosity)

    if (verbosity > 0) cat("Download is complete\n")
    return(df_final)
}
