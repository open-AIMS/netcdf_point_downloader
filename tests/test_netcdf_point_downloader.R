# Generate the warnings immediately so we can more easily determine which 
# test they correspond to.
options(warn = 1)

source('../netcdf_points_downloader.R')

capture_warnings_and_run <- function(expr) {
  warnings <- NULL
  result <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(result = result, warnings = warnings)
}

compare_vectors <- function(actual, expected, add_info = NULL, tolerance = 1e-5) {
  # Check lengths first
  if (length(actual) != length(expected)) {
    stop("Lengths of actual and expected vectors do not match.")
  }
  
  # Check for NA values alignment
  if (!all(is.na(actual) == is.na(expected))) {
    cat("Actual:\n")
    print(actual)
    cat("Expected:\n")
    print(expected)
    stop("NA values do not match between actual and expected vectors.")
  }
  
  # Replace NA with 0 for comparison, if necessary
  actual_clean <- ifelse(is.na(actual), 0, actual)
  expected_clean <- ifelse(is.na(expected), 0, expected)
  
  # Compare within tolerance
  differences <- abs(actual_clean - expected_clean)
  within_tolerance <- differences <= tolerance
  
  # If any values are not within tolerance, stop and throw an error
  # If any values are not within tolerance, stop and throw an error
  if (!all(within_tolerance)) {
    if (!is.null(add_info)) {
      cat("Additional Info:\n")
      print(add_info)
    }
    cat("Actual:\n")
    print(actual)
    cat("Expected:\n")
    print(expected)
    stop("Actual and expected values do not match within the specified tolerance.\n")
  }
  
  return(TRUE) # If no error is thrown, return TRUE
}

# Which tests to run. Use this to speed up the testing when developing new tests.
TESTS_TO_RUN <- c(
  TRUE, TRUE, TRUE, TRUE, TRUE,
  TRUE, TRUE, TRUE, TRUE) 

# ------------------  TEST 1 ---------------------
if (TESTS_TO_RUN[1]) {
  cat("Test 1: Extracting botz and eta, but no datetime specified\n")

  download_data_src_url_template <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr4_2.0.ncml'
  # Wrapping the function call in a tryCatch to catch the specific error
  result <- tryCatch({
    netcdf_points_downloader(var_names = c("botz", "eta"), 
                              input_data_src = "test-locations.csv", 
                              download_data_src_url_template, 
                              output_csv_path = "temp/test-1.csv", 
                              verbosity = 0)
    # If the function call does not throw an error, then the test has failed.
    stop("Test failed: Expected error was not thrown.")
  }, error = function(e) {
    # Check if the error message matches the expected pattern
    expected_message <- "Variables (eta) need a time dimension, but no time column found in input CSV."
    if(e$message == expected_message) {
      # Test passes if the error message contains the expected text
      cat("Test 1 passed: Expected error was caught.\n")
    } else {
      # The error was not what we expected, so the test fails
      stop(paste("Test failed: Unexpected error message:", e$message))
    }
  })
}



# ------------------  TEST 2 ---------------------
if (TESTS_TO_RUN[2]) {
  cat("Test 2: 12 points, from CSV, no date time, just botz.\n")
  output_csv_path = "temp/test-2.csv"
  ereefs_data <- netcdf_points_downloader(var_names = c("botz"), 
                              input_data_src = "test-locations.csv", 
                              download_data_src_url_template, 
                              output_csv_path, 
                              verbosity = 2,
                              stop_on_error = TRUE)
  # Was the data downloaded and saved
  if (!file.exists(output_csv_path)) {
    stop(paste("Test failed: output data not saved", output_csv_path))
  }

  # Does the output have NA in dry cells and values in wet cells
  # Find rows that should be wet for GBR 1
  wetcells <- ereefs_data[ereefs_data$WET_GBR4 == 1, ]

  # Fail if any of the botz column is NA
  if (any(is.na(wetcells$botz))) {
    stop("Test failed: NA values found in 'botz' column for wet cells.")
  }

  # Get the known dry or out of grid cells
  drycells <- ereefs_data[ereefs_data$WET_GBR4 == 0, ]

  # Fail if any of the botz column is no NA
  if (!all(is.na(drycells$botz))) {
    stop("Test failed: Non-NA values found in 'botz' column for dry cells.")
  }

  # Check for the presence of eReefs lat_val and lon_val columns
  if (!all(c("lat_val", "lon_val") %in% names(ereefs_data))) {
    stop("Test failed: eReefs lat_val and/or lon_val columns are missing.")
  }

  # Calculate the absolute differences
  lat_diff <- abs(ereefs_data$lat_val - ereefs_data$lat)
  lon_diff <- abs(ereefs_data$lon_val - ereefs_data$lon)

  # Identify points with differences greater than 0.05
  lat_issues <- which(lat_diff > 0.05)
  lon_issues <- which(lon_diff > 0.05)

  # Check botz values for those points
  if (any(!is.na(ereefs_data$botz[lat_issues])) || any(!is.na(ereefs_data$botz[lon_issues]))) {
    stop("Test failed: botz values should be NA for points with lat/lon differences greater than 0.05.")
  }

  cat("Test 2 passed.\n")
}




# ------------------  TEST 3 ---------------------
if (TESTS_TO_RUN[3]) {
  # Test 3 - Expect a helpful message
  cat("Test 3: request unknown variable against AIMS regular grid\n")
  download_data_src_url_template <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc"
  output_csv_path = "temp/test-3.csv"
  result <- tryCatch({
    netcdf_points_downloader(var_names = c("foobar"), 
                              input_data_src = "test-locations.csv", 
                              download_data_src_url_template, 
                              output_csv_path, 
                              verbosity = 0,
                              stop_on_error = TRUE)
    # If the function call does not throw an error, then the test has failed.
    stop("Test failed: Expected error was not thrown.")
  }, error = function(e) {
    # Check if the error message matches the expected pattern
    expected_message <- "The following variables are not available in dataset: foobar  available vars:  zc, mean_cur, salt, temp, u, v, mean_wspeed, eta, wspeed_u, wspeed_v"
    if(e$message == expected_message) {
      # Test passes if the error message contains the expected text
      cat("Test 3 passed: Expected error was caught.\n")
    } else {
      # The error was not what we expected, so the test fails
      stop(paste("Test failed: Unexpected error message:", e$message))
    }
  })
  cat("Test 3 passed.\n")
}




# ------------------  TEST 4 ---------------------
if (TESTS_TO_RUN[4]) {
  cat("Test 4: request multiple variables against AIMS regular grid\n")
  download_data_src_url_template <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc"
  output_csv_path = "temp/test-4.csv"
  input_df <- read.csv("test-locations.csv", stringsAsFactors = FALSE)
  input_df$datetime <- "2022/06/12 10:00"
  input_df$depth <- -2
  ereefs_data <- netcdf_points_downloader(var_names = c("temp", "eta"), 
                              input_data_src = input_df, 
                              download_data_src_url_template, 
                              output_csv_path, 
                              verbosity = 0,
                              stop_on_error = TRUE)
  temp_expected <- c(25.4402523, 25.31443214, 25.64244843, 26.2272644, NA, 26.54280281, 
    26.46940231, NA, NA, 25.45497894, 25.66849709, 24.61592484)

  eta_expected <- c(0.005879886, -0.001490339, 0.012014549, 0.018708244, NA, 0.1058442,
    0.111275263, NA, NA, -0.03240861, -0.041306648, 0.001989755)

  compare_vectors(ereefs_data$temp, temp_expected)
  compare_vectors(ereefs_data$eta, eta_expected)
  cat("Test 4 passed.\n")

}


# ------------------  TEST 5 ---------------------
if (TESTS_TO_RUN[5]) {
  cat("Test 5: Non eReefs data service - IMAS Seamap bathymetry\n")
  download_data_src_url_template <- "https://thredds.imas.utas.edu.au/thredds/dodsC/IMAS/SeamapAus_Bathymetry_National_50m/bathy_01_flt_clp.nc"
  #download_data_src_url_template <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc"
  #download_data_src_url_template = 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'
  #download_data_src_url_template <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr4_2.0.ncml'
  output_csv_path = "temp/test-5.csv"
  input_df <- read.csv("test-locations.csv", stringsAsFactors = FALSE)
  input_df$datetime <- "2022/06/12 10:00"
  input_df$depth <- -2

  ereefs_data <- netcdf_points_downloader(var_names = c("DEPTH"), 
                              input_data_src = input_df, 
                              download_data_src_url_template, 
                              output_csv_path, 
                              verbosity = 0,
                              stop_on_error = TRUE)
  DEPTH_expected <- c(-8.25, -3.5, -2, -55.5, -1162.5, -8.75, -10, NA, NA, -27.25, -33.5, 0)
  compare_vectors(ereefs_data$DEPTH, DEPTH_expected)
  cat("Test 5 passed.\n")
}




# ------------------  TEST 6 ---------------------
if (TESTS_TO_RUN[6]) {
  cat("Test 6: Out of time bounds behaviour.\n")
  # URLs and corresponding times
  # urls and times must have the same length
  urls <- c(
    "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc",
    'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml',
    'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr4_2.0.ncml'
  )
  times <- c("2000/06/12 10:00", "2000-06-12 10:00", "12-08-2000 10:00")

  # Iterate through the combinations of urls and times
  for (i in seq_along(urls)) {
    cat("Test 6-",i,": Out of temporal bounds.\n")
    
    # Define the expected warning start message dynamically
    expected_warning_start <- sprintf("Time (%s) for row 1 is out of range", times[i])
    
    # Prepare input data frame
    input_df <- data.frame(datetime = times[i], depth = c(-2), lon = c(147), lat = c(-19), stringsAsFactors = FALSE)

    # Run the function and capture warnings
    test_result <- capture_warnings_and_run({
      netcdf_points_downloader(var_names = c("temp", "eta"), 
                               input_data_src = input_df, 
                               download_data_src_url_template = urls[i], 
                               output_csv_path = NULL, #paste0("temp/test-6-", i, ".csv"), 
                               enable_disk_cache = FALSE,
                               verbosity = 0,
                               stop_on_error = TRUE)
    })
    
    # Check for expected warning
    warning_matched <- any(sapply(test_result$warnings, function(w) startsWith(w, expected_warning_start)))
    
    # Check if 'temp', 'eta', and 'time_val' are all NA in ereefs_data
    ereefs_data <- test_result$result
    data_checks_passed <- all(is.na(ereefs_data$temp[1]), is.na(ereefs_data$eta[1]), is.na(ereefs_data$time_val[1]))
    
    # Assert that both the warning was matched and the data checks passed
    if (warning_matched && data_checks_passed) {
      cat("Test 6-",i," passed: Expected warning captured and data checks out.\n")
    } else {
      print(test_result)
      stop("Test 6-",i," failed.\n")
    }
  }
}




# ------------------  TEST 7 ---------------------
if (TESTS_TO_RUN[7]) {
  # https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2_river_tracers/monthly.nc.html
  cat("Test 7: AIMS regular grid river tracer data\n")
  download_data_src_url_template <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr4_v2_river_tracers/monthly.nc"
  output_csv_path = "temp/test-7.csv"
  # Define start and end dates
  start_date <- as.Date("2019-01-01")
  end_date <- as.Date("2019-06-01")

  # Generate the sequence of dates, incrementing by month
  date_sequence <- format(seq.Date(from = start_date, to = end_date, by = "month"), '%Y-%m-%d %H:%M:%S')


  # Create a data frame with this date-time sequence
  input_df <- data.frame(datetime = date_sequence)
  input_df$depth <- -2

  # Cleveland bay
  input_df$lat <- -19.2
  input_df$long <- 146.9

  ereefs_data <- netcdf_points_downloader(var_names = c("bur", "hau"), # Burdekin, Haughton
                              input_data_src = input_df, 
                              download_data_src_url_template, 
                              output_csv_path, 
                              verbosity = 2,
                              stop_on_error = TRUE)
  bur_expected <- c(
      0.003952577,
      0.052879505,
      0.049450144,
      0.021398755,
      0.010989358,
      0.003535847
    )
  hau_expected <- c(
      0.000616914,
      0.013681924,
      0.002550242,
      0.001037028,
      0.000277827,
      4.84E-05
    )
  compare_vectors(ereefs_data$bur, bur_expected, ereefs_data)
  compare_vectors(ereefs_data$hau, hau_expected, ereefs_data)
      
  cat("Test 7 passed.\n")
}



# ------------------  TEST 8 ---------------------
if (TESTS_TO_RUN[8]) {
  # Test out of grid coordinates. In general these should be caught prior to
  # the requests to download the data.
  # There might be glitches in the out-of-bounds detection that might not align
  # with the server if the point is right near the boundary.
  # This is particularly complicated with the curvilinear grid where we determine
  # if a point is outside the grid if it is more than ~3 x the average distance 
  # between cells, if they were a regular grid. This means that for a squished
  # curvilinear grid the boundary detection might indicate that a point is inside
  # the bounds even though it is maybe 9 pixels from the edge. 
  # If a point gets through then when the request is made to the server it will
  # return a error, which will normally be converted to a warning and the result
  # recorded as an NA. The problem is that a network failure also results in
  # the same error message.
  cat("Test 8: Test out of bounds behaviour\n")

  ereefs_regular_url <- "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc"
  ereefs_curv_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'

  # Check if urls, lats, logs, in_bounds and temps must have the same length
  urls <- c(
    ereefs_regular_url,
    ereefs_regular_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url,
    ereefs_curv_url
  )

  lats <- c(
    -19, 
    25,
    -17.91410,
    -17.90083,
    -17.89512,
    -17.88315,
    -17.84900,
    25,         # Far out of bounds
    -28.6134,   # Just out of bounds on the southern side GBR1
    # Wet cells near PNG and the southern edge where cells are stretched, at intersection between cells.
    # If CURVILINEAR_GRID_SQUISH is too small then this will be considered out of bounds
    -7.8091,
    -8.0912,
    -28.5808
    )
  longs <- c(
    147, 
    25,
    # GBR 1
    148.38817, # Last wet pixel
    148.40066, # Last null pixel
    148.40461, # One pixel past null pixel
    148.41358, # Three pixels past null outer edge
    148.45108, # About 10 pixels past the edge
    25,        # Far out of bounds
    153.6781,  # Just out of bounds on the southern side GBR1
    144.8657,  # Wet cells at the edge
    145.9708,
    153.7137
    )
  in_bounds <- c(
    TRUE, 
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,     # Even though we are 3 pixels out of bounds, need to cater for CURVILINEAR_GRID_SQUISH
    FALSE,
    FALSE,    # Far out of bounds
    TRUE,     # Still just within the distance_threshold
    TRUE,     # Wet cells at the edge
    TRUE,
    TRUE)
  temps <- c(
    25.75602, 
    NA,
    26.1304,
    NA,
    NA,
    NA,     
    NA,
    NA,     # Far out of bounds
    22.43302,
    29.849951,  # Wet cells at the edge
    29.960915,
    22.62871)


  # Iterate through the combinations of urls and times
  for (i in seq_along(urls)) {
    cat("Test 8-",i,": Out of spatial bounds.\n")

    input_df <- data.frame(datetime = "2022/06/12 10:00", depth = -2, lat = lats[i], lon = longs[i])
    #print(input_df)
    # Run the function and capture warnings
      test_result <- capture_warnings_and_run({
        netcdf_points_downloader(var_names = c("temp"), 
                                 input_data_src = input_df, 
                                 download_data_src_url_template = urls[i], 
                                 output_csv_path = NULL,  
                                 enable_disk_cache = FALSE,
                                 verbosity = 0,
                                 stop_on_error = TRUE)
      })
    
    #print(test_result)
    # Check for expected warning
    warning_matched <- any(sapply(test_result$warnings,   
      function(w) startsWith(w, "Location is outside the acceptable distance")))
    vals <- test_result$result
    actual_vector <- c(vals$temp, vals$lat_val, vals$lon_val)
    if (in_bounds[i]) {
      expected_vector <- c(temps[i], lats[i], longs[i])
    } else {
      expected_vector <- c(temps[i], NA, NA)
    }
    compare_vectors(actual_vector, expected_vector, ereefs_data, tolerance = 3e-2)
    
    #print(warning_matched)
  }

  cat("Test 8 passed.\n")
}


# ------------------  TEST 9 ---------------------
if (TESTS_TO_RUN[9]) {
  cat("Test 9: Do individual eReefs opendap files match the aggregate services?\n")
  # In this test we download data from the individual eReefs end points and
  # compare it with the aggregate data services. This will test whether the
  # the download_data_src_url_template will work.

  aggregate_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr4_2.0.ncml'
  individual_url_template <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_%(year)-%(month)-%(day).nc'

  # Create a time series that will span multiple end points and have multiple
  # requests per end point. This is to test the reuse of the data source object
  # across requests.
  start_date <- as.POSIXct("2021-08-12 01:00", tz ="ETC/GMT-10")
  num_days <- 4
  num_intervals_per_day <- 4 # 24 hours / 6 hours

  # Calculate total number of intervals
  total_intervals <- num_days * num_intervals_per_day

  # Generate the series
  date_vector <- seq(start_date, by = "6 hours", length.out = total_intervals)
  input_df <- data.frame(datetime = date_vector, lat = -19, lon = 147, depth = -8)

  individual_data <- netcdf_points_downloader(var_names = c("eta"), 
                                 input_data_src = input_df, 
                                 download_data_src_url_template = individual_url_template, 
                                 output_csv_path = 'temp/test-9-individual.csv',  
                                 enable_disk_cache = TRUE,
                                 verbosity = 2,
                                 stop_on_error = TRUE)
                                 
  aggregate_data <- netcdf_points_downloader(var_names = c("eta"), 
                                 input_data_src = input_df, 
                                 download_data_src_url_template = aggregate_url, 
                                 output_csv_path = 'temp/test-9-aggregate.csv',  
                                 enable_disk_cache = TRUE,
                                 verbosity = 2,
                                 stop_on_error = TRUE)
                                 
  compare_vectors(individual_data$eta, aggregate_data$eta)
  
  DEBUGGING <- FALSE
  if (DEBUGGING) {
    # Both methods should have downloaded the same data.
    # Plot tideHeight
    plot(individual_data$datetime, individual_data$eta, type = 'l', col = 'blue', 
       xlab = 'DateTime', ylab = 'Tide (m LAT)', 
       main = 'Individual data and aggregate data service')
    # Add eta to the existing plot
    lines(aggregate_data$datetime, aggregate_data$eta, col = 'red')
  }
  cat("Test 9 passed.\n")
}