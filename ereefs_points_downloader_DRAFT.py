# This version of the tool is incomplete. The tool was originally
# developed in Python, then ported to R, then extended and refined.
# We have not back ported the changes to Python. I have left the 
# Python code here as a starting point for extending.
# There are a lot of functions that are marked as deprecated.
# These functions were developed, then refactoring was done, or
# a simplification was developed making these functions no longer
# needed. Unfortunately the development was not fully completed
# to the point where it was clear that this code would not longer
# be used.
# One key feature that was deprecated was having a grid lookup to
# speed up the calculation of the spatial indices. This did speed
# up this part of the calculation (about 4 - 10 times), but the 
# download was still the dominant processing time. I had trouble
# porting this feature to R and so decided to deprecate it for a
# slower but simplier method.

import pandas as pd
import numpy as np
import xarray as xr
import os
import time
from datetime import datetime, timedelta
import json

# Constants
#input_csv_path = 'original/TIDES_input-data-for-eReef_nGBR2023.csv'
input_csv_path = 'derived/Simulated_BROKEN_TIDES_input-data-for-eReef_nGBR2023.csv'
output_csv_path = 'output/TIDES_eReefs-GBR1-Hydro_v2_nGBR2023.csv'

input_csv_path = "tests/test_locations.csv"
#output_csv_path = "tests/test_locations_out_Python_GBR4.csv"
output_csv_path = "tests/test_locations_out_Python_GBR1.csv"

#opendap_url = "https://thredds.ereefs.aims.gov.au/thredds/dodsC/gbr1_2.0/daily.nc"
opendap_url = 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'
#opendap_url = 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr4_2.0.ncml'
variables = ['botz']

download_rate_limit_delay = 0.1

"""
This script downloads data from eReefs data services at the locations and times
specified in the input_csv_path file. The variables specified are downloaded from
OpenDAP service specified by opendap_url. 

To make the download process robust, this script creates a copy of the input data, 
then routinely saves processing and download progress to this file (data_cache). 
If the script is restarted then it resumes processing using the data_cache. If the
input location data (input_csv_path) or set of variables is changed then data_cache
is reset.

Once all the data is downloaded the final output file is created by stripping off 
the additional intermediate progress tracking attributes.

This script automatically detects a variety of different column names for latitude,
longitude, time and depth. It also adapts to different data time formats.
"""

    
def create_spatial_grid_flatten_deprecated(lat_array, long_array, grid_size=(1, 1)):
    """
    Create a grid corresponding to the bounding extents of the data.
    In each grid cell store the indicies of all the matching curvilinear
    cells. When we are trying to find the closed location, we only need to
    search all the elements in one grid cell.
    This speeds up the index look ups considerably.
    No grid: 2 m 40 sec
    50x50 Grid: 19 sec
    """
    lat_min, lat_max = np.nanmin(lat_array), np.nanmax(lat_array)
    long_min, long_max = np.nanmin(long_array), np.nanmax(long_array)

    lat_step = (lat_max - lat_min) / grid_size[0]
    long_step = (long_max - long_min) / grid_size[1]

    # Initialize an empty grid
    grid = [[[] for _ in range(grid_size[1])] for _ in range(grid_size[0])]

    # Populate the grid with indices
    for index, (lat, long) in enumerate(zip(lat_array.flatten(), long_array.flatten())):
        if np.isnan(lat) or np.isnan(long):
            continue  # Skip NaN values
        lat_idx = int((lat - lat_min) / lat_step)
        long_idx = int((long - long_min) / long_step)
        # Ensure indices fall within the grid
        lat_idx = min(max(lat_idx, 0), grid_size[0] - 1)
        long_idx = min(max(long_idx, 0), grid_size[1] - 1)
        grid[lat_idx][long_idx].append(index)

    return grid, (lat_min, long_min, lat_step, long_step)
    
def create_spatial_grid_deprecated(lat_array, long_array, grid_size=(100, 100)):
    """
    Create a low resolution regular grid covering the extent of the
    curvilinear grid (lat_array, long_array), where each cell of the
    regular grid has the indicies of the curvilinear grid that cover that
    region. This makes finding the closest curvilinear cell to specified
    location much faster as we can quickly look up in the small set of 
    cells that we need to perform the nearest neighbour checks on.
    The performance improvements plateaus at around a 50x50 grid
    Parameters:
    - lat_array: 2D vector of the centroids of the latitudes of each cell 
    in the curvilinear grid.
    - long_array: 2D vector of the centroids of the longitudes of each cell 
    in the curvilinear grid.
    Return:
    - grid: 2D grid where each cell contains a list of 
    """
    
    # The curvilinear grid can contain NaN cells as not all parts
    # of the lat_array and long_array are valid cells. For eReefs GBR1
    # there are cut outs. To work out the bounds of the valid grid
    # we need to ignore the NaN cells.
    lat_min, lat_max = np.nanmin(lat_array), np.nanmax(lat_array)
    long_min, long_max = np.nanmin(long_array), np.nanmax(long_array)

    lat_step = (lat_max - lat_min) / grid_size[0]
    long_step = (long_max - long_min) / grid_size[1]
    
    original_shape = lat_array.shape  # Store the original shape

    # Flatten the arrays once and store them
    flat_lat_array = lat_array.flatten()
    flat_long_array = long_array.flatten()

    # Initialize an empty grid
    grid = [[[] for _ in range(grid_size[1])] for _ in range(grid_size[0])]

    # Populate the grid with indices, using the pre-flattened arrays
    for index, (lat, long) in enumerate(zip(flat_lat_array, flat_long_array)):
        if np.isnan(lat) or np.isnan(long):
            continue  # Skip NaN values
        lat_idx = int((lat - lat_min) / lat_step)
        long_idx = int((long - long_min) / long_step)
        lat_idx = min(max(lat_idx, 0), grid_size[0] - 1)
        long_idx = min(max(long_idx, 0), grid_size[1] - 1)
        grid[lat_idx][long_idx].append(index)

    # Return the grid along with the flat arrays and grid parameters.
    # Do this 
    return grid, (lat_min, long_min, lat_step, long_step), flat_lat_array, flat_long_array, original_shape

    
def find_nearest_curvilinear_index_deprecated(flat_lat_array, flat_long_array, lat, long, spatial_grid, grid_params, original_shape):
    lat_min, long_min, lat_step, long_step = grid_params
    grid_size = (len(spatial_grid), len(spatial_grid[0]))

    # Find the grid cell of the query point
    lat_idx = int((lat - lat_min) / lat_step)
    long_idx = int((long - long_min) / long_step)
    lat_idx = min(max(lat_idx, 0), grid_size[0] - 1)
    long_idx = min(max(long_idx, 0), grid_size[1] - 1)

    # Get the indices in the cell
    indices = spatial_grid[lat_idx][long_idx]
    if not indices:
        raise ValueError("No points in the vicinity. Consider searching in neighboring cells.")

    # Use the pre-flattened arrays for the search
    subset_lat = flat_lat_array[indices]
    subset_long = flat_long_array[indices]
    difference_array = ((subset_lat - lat) ** 2) + ((subset_long - long) ** 2)
    nearest_index = indices[np.nanargmin(difference_array)]

    return np.unravel_index(nearest_index, original_shape)


def find_nearest_curvilinear_index_simple_deprecated(lat_array, long_array, lat, long):
    """
    Returns a tuple corresponding to the i, j index.
    Raises ValueError if input lat/long is out of the bounds of the input data or points to a dry cell.
    """
    # Find the vector distance squared between the point of interest and the curvilinear arrays. 
    difference_array = ((lat_array-lat)**2) + ((long_array-long)**2)
    # Find the closest point.
    index = np.unravel_index(np.nanargmin(difference_array), difference_array.shape)
    return index
    


def write_metadata_deprecated(metadata_path, df_metadata):
    """Writes metadata to a JSON file."""
    with open(metadata_path, 'w') as file:
        json.dump(df_metadata, file)

def read_metadata_deprecated(metadata_path):
    """Reads metadata from a JSON file, returns an empty dictionary if the file doesn't exist."""
    if os.path.exists(metadata_path):
        with open(metadata_path, 'r') as file:
            return json.load(file)
    return {}

def setup_data_cache_deprecated(input_csv_path, data_cache_path, ds, variables):
    metadata_path = data_cache_path + '.json'
    df_metadata = read_metadata(metadata_path)
    input_csv_mod_date = get_file_modification_date(input_csv_path)
    
    # Validate cache: Check modification date and variables
    if df_metadata:
        metadata_mod_date = df_metadata.get('modification_date', '')
        metadata_vars = df_metadata.get('variables', [])
        
        if metadata_mod_date != input_csv_mod_date or set(metadata_vars) != set(variables):
            print("Cache is invalid due to modification date or variables change. Setting up cache from scratch.")
            df_metadata = {}  # Invalidate metadata to trigger cache reset

    if df_metadata and os.path.exists(data_cache_path):
        print(f"Found existing data cache: {data_cache_path}")
        # Data cache is valid, read it
        df = pd.read_csv(data_cache_path)
    else:
        # Cache is invalid or does not exist, create from input CSV
        df = pd.read_csv(input_csv_path)
        
        # Initialise columns for recording progress
        for col in variables + ['time_idx', 'lat_idx', 'lon_idx']:
            df[col] = float('nan')
            
        for col in ['index_complete', 'vars_download']:
            df[col] = False
        
        # Update metadata
        df_metadata = {'modification_date': input_csv_mod_date, 'variables': variables}

        # Initialize additional processing
        identify_column_names_and_format(df, df_metadata, ds, variables)
        
        save_data_cache(df, df_metadata, data_cache_path)
        write_metadata(metadata_path, df_metadata)
        
    return df, df_metadata  # Return both DataFrame and metadata dictionary
    
def precalculate_indices_deprecated(df, df_metadata, ds, data_cache_path):
    """
    Precalculate the nearest indices for latitude, longitude, and time for each row in df,
    and record if cells are wet or dry based on 'eta'. Continue from the last unprocessed row.
    Saves the data cache every 1000 calculations and prints out absolute progress.
    
    Parameters:
    - df: The DataFrame to process.
    - ds: The data structure object representing the eReefs OpenDAP endpoint.
    - data_cache_path: Path to save the data cache periodically.
    """
    if df['index_complete'].all():
        print("All indices have already been calculated.")
        return df

    if 'wet_dry_status' not in df.columns:
        df['wet_dry_status'] = 'unprocessed'

    # If we have previously done some processing (a restart) indicate where we are 
    # resuming from. 
    already_processed = df['index_complete'].sum()
    if already_processed > 0:
        total_rows = len(df)
        print(f"Resuming index calculation from row {already_processed} of {total_rows}.")

    start_time = datetime.now()
    
    # Prepare a lookup array for time from the OpenDAP time dimension
    time_required = df_metadata['datetime_column'] is not None
    if time_required:
        time_array = pd.to_datetime(ds['time'].values)
    
    depth_required = df_metadata['depth_column'] is not None
    if depth_required:
        depth_array = ds['zc'].values  # Assuming 'zc' is the depth dimension name
    
    lat_array, long_array = ds['latitude'].values, ds['longitude'].values
    
    print("Downloading single time slice of eta to determine if cells are wet or dry (90MB) ...")
    eta_slice = ds['eta'].isel(time=0).values

    # Determine if the dataset represents a regular grid based on the shape of latitude and longitude arrays.
    regular_grid = lat_array.ndim == 1 and long_array.ndim == 1

    # -------- Create a fast lookup table for curvilinear grids ------
    if not regular_grid:
        print("Calculating coordinate lookup grid for curvilinear data")
        # Calculate a fast lookup grid (lut). This is a regular grid, making it
        # fast to find the appropriate cell for a given location. It
        # contains a list of all the curvilinear indices that should be
        # checked to find the closest point.
        
        lut_grid_size=(100, 100)
        
        # The curvilinear grid can contain NaN cells as not all parts
        # of the lat_array and long_array are valid cells. For eReefs GBR1
        # there are cut outs. To work out the bounds of the valid grid
        # we need to ignore the NaN cells.
        lat_min, lat_max = np.nanmin(lat_array), np.nanmax(lat_array)
        long_min, long_max = np.nanmin(long_array), np.nanmax(long_array)

        lat_step = (lat_max - lat_min) / lut_grid_size[0]
        long_step = (long_max - long_min) / lut_grid_size[1]
        
        original_shape = lat_array.shape  # Store the original shape

        # Flatten the arrays once and store them. Flatten operations is slow.
        flat_lat_array = lat_array.flatten()
        flat_long_array = long_array.flatten()

        # Initialize an empty lut_grid
        lut_grid = [[[] for _ in range(lut_grid_size[1])] for _ in range(lut_grid_size[0])]

        # Populate the grid with indices, using the pre-flattened arrays
        for index, (lat, long) in enumerate(zip(flat_lat_array, flat_long_array)):
            if np.isnan(lat) or np.isnan(long):
                continue  # Skip NaN values
            lut_lat_idx = int((lat - lat_min) / lat_step)
            lut_long_idx = int((long - long_min) / long_step)
            lut_lat_idx = min(max(lut_lat_idx, 0), lut_grid_size[0] - 1)
            lut_long_idx = min(max(lut_long_idx, 0), lut_grid_size[1] - 1)
            lut_grid[lut_lat_idx][lut_long_idx].append(index)

    print("Calculating indices for input locations")
    # Find the matching indices for each of the input locations. Saves these
    # in lat_idx, lon_idx, time_idx and depth_idx columnes.
    for index, row in df[df['index_complete'] == False].iterrows():
        input_lat = row[df_metadata['latitude_column']]
        input_lon = row[df_metadata['longitude_column']]
        
        if regular_grid:
            # Find the nearest index for the input location on a regular grid
            lat_idx = np.abs(lat_array - input_lat).argmin()
            lon_idx = np.abs(long_array - input_lon).argmin()
            
        else:
            # Process indices for curvilinear grid
            #grid, grid_params, flat_lat_array, flat_long_array, original_shape = create_spatial_grid(lat_array, long_array)

            #lat_idx, lon_idx = find_nearest_curvilinear_index(flat_lat_array, flat_long_array, input_lat, input_lon, grid, grid_params, original_shape)
            
            # ------- Find nearest curvilinear index -------

            # Find the grid cell of the query point
            lut_lat_idx = int((input_lat - lat_min) / lat_step)
            lut_long_idx = int((input_lon - long_min) / long_step)
            lut_lat_idx = min(max(lut_lat_idx, 0), lut_grid_size[0] - 1)
            lut_long_idx = min(max(lut_long_idx, 0), lut_grid_size[1] - 1)

            # Get the indices in the cell
            indices = lut_grid[lut_lat_idx][lut_long_idx]
            if not indices:
                raise ValueError("No points in the vicinity. Consider searching in neighboring cells.")

            # Use the pre-flattened arrays for the search
            subset_lat = flat_lat_array[indices]
            subset_long = flat_long_array[indices]
            difference_array = ((subset_lat - input_lat) ** 2) + ((subset_long - input_lon) ** 2)
            nearest_index = indices[np.nanargmin(difference_array)]
            lat_idx, lon_idx = np.unravel_index(nearest_index, original_shape)
        
        df.at[index, 'lat_idx'] = lat_idx
        df.at[index, 'lon_idx'] = lon_idx
        
        
        if time_required:
            observation_time = datetime.strptime(row[df_metadata['datetime_column']], df_metadata['datetime_format'])
            time_idx = np.abs(time_array - observation_time).argmin()
            df.at[index, 'time_idx'] = time_idx

        if depth_required:
            observation_depth = row[df_metadata['depth_column']]
            depth_idx = np.abs(depth_array - observation_depth).argmin()
            df.at[index, 'depth_idx'] = depth_idx

        df.at[index, 'index_complete'] = True
        df.at[index, 'wet_dry_status'] = 'dry' if np.isnan(eta_slice[lat_idx, lon_idx]) else 'wet'

        processed_count = index + 1
        if processed_count % 2000 == 0 or processed_count == len(df):
            save_data_cache(df, df_metadata, data_cache_path)
            print(f"Processed {processed_count} of {len(df)} rows. Continuing...")

    elapsed_time = datetime.now() - start_time
    print(f'Elapsed time: {elapsed_time}. Processing look up table index complete.')
    return df

def report_out_of_bounds_deprecated(df, data_cache_path):
    """
    Reports on the number of dry cells based on the 'wet_dry_status' column in the DataFrame.
    Then provides the path of the data cache CSV.

    Parameters:
    - df: The DataFrame to analyze. 
    - data_cache_path: The path to the data cache CSV file.
    """
    if 'wet_dry_status' in df.columns:
        # Count the number of dry cells
        dry_cells_count = (df['wet_dry_status'] == 'dry').sum()
        
        if dry_cells_count > 0:
            print(f"Number of out-of-bounds or dry cells: {dry_cells_count}")
            print(f"The specific problematic rows can be seen in: {data_cache_path}")
        else:
            print("No dry cells found.")
    else:
        print("The 'wet_dry_status' column does not exist in the DataFrame.")

def find_nearest_index(lat_array, long_array, lat, long):
    """
    Returns a tuple corresponding to the i, j index for both regular and curvilinear grids.
    Raises ValueError if input lat/long is out of the bounds of the input data or points to a dry cell.
    """
    if lat_array.ndim == 1 and long_array.ndim == 1:
        # Regular grid
        lat_diff = lat_array[:, None] - lat  # Column vector
        long_diff = long_array[None, :] - long  # Row vector
        difference_array = lat_diff**2 + long_diff**2
    elif lat_array.ndim == 2 and long_array.ndim == 2:
        # Curvilinear grid
        difference_array = (lat_array - lat)**2 + (long_array - long)**2
    else:
        raise ValueError("lat_array and long_array must both be either 1D or 2D arrays")
    
    # Find the closest point.
    index_flat = np.nanargmin(difference_array)
    index = np.unravel_index(index_flat, difference_array.shape)
    
    return index


def precalculate_indices(df, df_metadata, ds, data_cache_path):
    """
    Precalculate the nearest indices for latitude, longitude, and time for each row in df,
    and record if cells are wet or dry based on 'eta'. Continue from the last unprocessed row.
    Saves the data cache every 1000 calculations and prints out absolute progress.
    
    Parameters:
    - df: The DataFrame to process.
    - ds: The data structure object representing the eReefs OpenDAP endpoint.
    - data_cache_path: Path to save the data cache periodically.
    """
    if df['index_complete'].all():
        print("All indices have already been calculated.")
        return df

    # If we have previously done some processing (a restart) indicate where we are 
    # resuming from. 
    already_processed = df['index_complete'].sum()
    if already_processed > 0:
        total_rows = len(df)
        print(f"Resuming index calculation from row {already_processed} of {total_rows}.")

    start_time = datetime.now()
    
    # Prepare a lookup array for time from the OpenDAP time dimension
    time_required = df_metadata['datetime_column'] is not None
    if time_required:
        time_array = pd.to_datetime(ds['time'].values)
    
    depth_required = df_metadata['depth_column'] is not None
    if depth_required:
        depth_array = ds['zc'].values  # Assuming 'zc' is the depth dimension name
    
    lat_array, long_array = ds['latitude'].values, ds['longitude'].values

    # Determine if the dataset represents a regular grid based on the shape of latitude and longitude arrays.
    regular_grid = lat_array.ndim == 1 and long_array.ndim == 1


    print("Calculating indices for input locations")
    # Find the matching indices for each of the input locations. Saves these
    # in lat_idx, lon_idx, time_idx and depth_idx columnes.
    for index, row in df[df['index_complete'] == False].iterrows():
        input_lat = row[df_metadata['latitude_column']]
        input_lon = row[df_metadata['longitude_column']]
            
        lat_idx, lon_idx = find_nearest_index(lat_array, long_array, input_lat, input_lon)
        
        df.at[index, 'lat_idx'] = lat_idx
        df.at[index, 'lon_idx'] = lon_idx

        if time_required:
            observation_time = datetime.strptime(row[df_metadata['datetime_column']], df_metadata['datetime_format'])
            time_idx = np.abs(time_array - observation_time).argmin()
            df.at[index, 'time_idx'] = time_idx

        if depth_required:
            observation_depth = row[df_metadata['depth_column']]
            depth_idx = np.abs(depth_array - observation_depth).argmin()
            df.at[index, 'depth_idx'] = depth_idx

        df.at[index, 'index_complete'] = True

        processed_count = index + 1
        if processed_count % 200 == 0 or processed_count == len(df):
            save_data_cache(df, df_metadata, data_cache_path)
            print(f"Processed {processed_count} of {len(df)} rows. Continuing...")

    elapsed_time = datetime.now() - start_time
    print(f'Elapsed time: {elapsed_time.total_seconds():.2f} sec. Processing look up table index complete.')
    return df





    

def download_data_with_precalculated_indices(df, df_metadata, ds, data_cache_path, download_rate_limit_delay):
    """
    Download each row of data from the input CSV (represented by the data cache
    version of it in df) from the OpenDAP service (ds). Progress is saved
    every 10 rows, so the script is restartable.
    """
    report_rate = 10
    start_time = datetime.now()

    # Variables to download
    variables = df_metadata['variables']

    total_rows = len(df)
    
    # Determine the starting point for downloads
    already_downloaded = df['vars_download'].sum()

    # Report resuming or starting anew
    if already_downloaded > 0:
        print(f"Resuming download at row {already_downloaded} of {total_rows}.")
    else:
        print("Starting data download using precalculated indices...")
        
    time_required = df_metadata['datetime_column'] is not None
    depth_required = df_metadata['depth_column'] is not None
    
    # Download for all input rows that haven't already been processed
    for index, row in df[df['vars_download'] == False].iterrows():

        # Use precalculated indices
        lat_idx = int(row['lat_idx'])
        lon_idx = int(row['lon_idx'])
        if time_required:
            time_idx = int(row['time_idx'])
        if depth_required:
            depth_idx = int(row['depth_idx'])

        for var in variables:
            try:
                # Make OpenDAP request for each variable in row. Adapt
                # to the dimensions that the variable uses.
                if time_required and depth_required:
                    # Both time and depth indices are used
                    var_value = ds[var][time_idx, depth_idx, lat_idx, lon_idx].values.item()

                elif time_required and not depth_required:
                    # Only time index is used, assuming the variable does not have a depth dimension
                    var_value = ds[var][time_idx, lat_idx, lon_idx].values.item()

                elif not time_required and depth_required:
                    # Only depth index is used, assuming the variable does not have a time dimension
                    var_value = ds[var][depth_idx, lat_idx, lon_idx].values.item()

                else:
                    # Neither time nor depth indices are used, assuming the variable has neither dimension
                    var_value = ds[var][lat_idx, lon_idx].values.item()
                    
                df.at[index, var] = var_value  # Update DataFrame with the downloaded value
            except Exception as e:
                print(f"Error downloading {var} for row {index}: {e}")
                df.at[index, var] = float('nan')  # Set to NaN in case of error

        # Mark the row to indicate a download attempt has been made
        df.at[index, 'vars_download'] = True

        # Pace the downloads if necessary
        time.sleep(download_rate_limit_delay)

        # Calculate progress and estimate completion time every 10 rows
        if (index + 1 - already_downloaded) % report_rate == 0: 
            current_progress = index + 1
            
            # Estimate completion time
            elapsed_time = datetime.now() - start_time
            start_time = datetime.now()
            rows_left = total_rows - current_progress
            estimated_time_remaining = elapsed_time * rows_left / report_rate
            estimated_completion_time = datetime.now() + estimated_time_remaining

            print(f"Downloading: {current_progress} of {total_rows} rows complete. "
                f"Estimated completion time: {estimated_completion_time.strftime('%Y-%m-%d %H:%M:%S')}")

            # Save progress
            save_data_cache(df, df_metadata, data_cache_path)

    # Final save and report
    save_data_cache(df, df_metadata, data_cache_path)
    print(f"Data download completed. Results saved to {data_cache_path}. Total time: {datetime.now() - start_time}")


def get_file_modification_date(file_path):
    """Returns the modification date of a file."""
    modification_time = os.path.getmtime(file_path)
    return datetime.fromtimestamp(modification_time).strftime('%Y-%m-%d %H:%M:%S')



def setup_data_cache(input_csv_path, data_cache_path, ds, variables):

    df, df_metadata = read_data_cache(data_cache_path)

    input_csv_mod_date = get_file_modification_date(input_csv_path)
    
    # Check if cache is valid Check modification date and variables
    if df_metadata:
        if df_metadata['modification_date'] != input_csv_mod_date or set(df_metadata['variables']) != set(variables):
            print("Cache is invalid due to modification date or variables change.")
            df_metadata = {}  # Invalidate metadata to trigger cache reset

    # Invalid cache or first time, so setup the cache
    if not df_metadata:
        print("Setting up data cache: "+data_cache_path)
        # Cache is invalid or does not exist, base the cache on the input CSV
        df = pd.read_csv(input_csv_path)
        
        # Initialise columns for recording progress
        for col in variables + ['time_idx', 'lat_idx', 'lon_idx', 'depth_idx']:
            df[col] = float('nan')
            
        for col in ['index_complete', 'vars_download']:
            df[col] = False
        
        df_metadata['modification_date'] = input_csv_mod_date
        df_metadata['variables'] = variables

        # Work out the alignment between column names in the input CSV and internal names
        # This adds additional values to df_metadata
        identify_column_names_and_format(df, df_metadata, ds, variables)
        
        save_data_cache(df, df_metadata, data_cache_path)
        
    return df, df_metadata  # Return both DataFrame and metadata dictionary
    
    


def save_data_cache(df, df_metadata, data_cache_path):
    # Save DataFrame to CSV
    df.to_csv(data_cache_path, index=False)
    
    # Save metadata to JSON
    metadata_path = data_cache_path + '.json'
    with open(metadata_path, 'w') as file:
        json.dump(df_metadata, file)


def read_data_cache(data_cache_path):
    """
    Read the cache files if they exist. 
    
    Technically this function wastes loading the data cache, when it might be
    immediately invalidated when we check if the modification date of the
    input has changed. We do it this way for simplicity.
    """
    metadata_path = data_cache_path + '.json'
    if not os.path.exists(metadata_path) or not os.path.exists(data_cache_path):
        return None, {}
        
    with open(metadata_path, 'r') as file:
        df_metadata = json.load(file)

    df = pd.read_csv(data_cache_path)
    print(f"Found existing data cache: {data_cache_path}")
    
    return df, df_metadata

def convert_cache_to_final_csv(df, output_csv_path):
    """
    Converts the data cache DataFrame into the final saved CSV file by removing
    index attributes (lat_idx, lon_idx, time_idx) and attributes used for tracking
    progress ('vars_download' and 'index_complete'), and saves it to the specified path.

    Parameters:
    - df: The data cache DataFrame.
    - output_csv_path: Path where the final CSV file should be saved.
    """
    # Columns to remove
    columns_to_remove = ['lat_idx', 'lon_idx', 'time_idx', 'vars_download', 'index_complete', 'attempt_status']
    
    # Drop specified columns if they exist in the DataFrame
    df_final = df.drop(columns=[col for col in columns_to_remove if col in df.columns], errors='ignore')
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    
    # Save the cleaned DataFrame to the specified CSV file path
    df_final.to_csv(output_csv_path, index=False)
    
    print(f"Final data saved to {output_csv_path}")

def identify_datetime_format(df, datetime_column, datetime_formats):
    """
    Identifies the datetime format by checking a random sample of rows from a specified column.
    
    Parameters:
    - df: DataFrame containing the datetime column.
    - datetime_column: The name of the column containing datetime strings.
    - datetime_formats: A list of strings representing datetime formats to try.
    
    Returns:
    - The matching datetime format string.
    
    Raises:
    - ValueError: If no matching datetime format is found.
    """
    # Ensure the column exists in the DataFrame
    if datetime_column not in df.columns:
        raise ValueError(f"The column '{datetime_column}' does not exist in the DataFrame.")
    
    # Sample up to 10 non-na datetime values from the dataframe
    sample_datetimes = df[datetime_column].dropna().sample(min(len(df), 10))
    
    for fmt in datetime_formats:
        try:
            # Try to convert each sampled datetime with the current format
            converted_sample = pd.to_datetime(sample_datetimes, format=fmt, errors='coerce')
            # If no NA values are produced, consider it a match
            if not converted_sample.isna().any():
                print(f"Match found for datetime format: {fmt}")
                return fmt
        except (ValueError, TypeError):
            # If an exception is raised, try the next format
            continue
    
    # If no format matched, raise an error
    raise ValueError(f"No match found for datetime format. Must be one of {', '.join(datetime_formats)}")


def identify_column_names_and_format(df, df_metadata, ds, variables_to_download):
    """
    Enhances the script to be tolerant of different attribute names for latitude, longitude, datetime, and depth.
    Automatically determines the appropriate datetime format within a fixed set. 
    Checks required dimensions for variables to be downloaded and ensures necessary 
    columns are available in the input file. Also verifies that all variables to 
    download exist in the dataset.
    """
    possible_lat_names = ['lat', 'latitude', 'LAT', 'LATITUDE', 'Y', 'y']
    possible_lon_names = ['long', 'longitude', 'lon','LONG', 'LONGITUDE', 'LON', 'X', 'x']
    possible_datetime_names = ['datetime', 'date', 'Date', 'time', 'Time']
    possible_depth_names = ['depth', 'zc', 'Depth', 'depth_m']
    datetime_formats = [
        '%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M', 
        '%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M',
        '%Y/%m/%d %H:%M:%S', '%Y/%m/%d %H:%M',
        '%d-%m-%Y %H:%M:%S', '%d-%m-%Y %H:%M', 
        '%d/%m/%Y %H:%M:%S', '%d/%m/%Y %H:%M'
    ]
    
    # Find the column name used in the input from the possible list of
    # expected column names.
    def find_column(column_names):
        return next((name for name in column_names if name in df.columns), None)

    # Identify required columns
    df_metadata['latitude_column'] = find_column(possible_lat_names)
    df_metadata['longitude_column'] = find_column(possible_lon_names)
    df_metadata['datetime_column'] = find_column(possible_datetime_names)
    df_metadata['depth_column'] = find_column(possible_depth_names)
    
    # Verify that the variables to download exist in the dataset
    missing_vars = [var for var in variables_to_download if var not in ds.variables]
    if missing_vars:
        raise ValueError(f"The following variables are not available in the dataset: {', '.join(missing_vars)}")

    # Make sure that if we have variables that have time or depth dimensions
    # that we have an input datetime and/or depth columns. If not print an error.
    
    # Check if time or depth dimensions are required for any of the variables to be downloaded
    time_required = depth_required = False
    df_metadata['time_vars'] = []
    df_metadata['depth_vars'] = []
    for var in variables_to_download:
        if 'time' in ds[var].dims:
            time_required = True
            df_metadata['time_vars'].append(var)
        if any(d in ds[var].dims for d in ['zc', 'depth', 'k']):
            depth_required = True
            df_metadata['depth_vars'].append(var)
    
    # If time or depth required, then report an error if no matching column found.
    errors = []
    if time_required and not df_metadata.get('datetime_column'):
        errors.append("a time dimension, but no time column found in input CSV")
    if depth_required and not df_metadata.get('depth_column'):
        errors.append("a depth dimension, but no depth column found in input CSV")
    
    if errors:
        raise ValueError("Dimensions missing from input CSV: " + "; ".join(errors))
    
    
    # Make sure that the input CSV has columns for latitute and longitude
    errors = []
    if df_metadata["latitude_column"] is None:
        errors.append("latitude column")
    if df_metadata["longitude_column"] is None:
        errors.append("longitude column")

    if errors:
        error_message = "Match for "+', '.join(errors)+" in input CSV not found."
        raise ValueError(error_message)
    
    # Identify datetime format of the input file by iteratively trying to convert
    # the input date with a specific datetime format. If it fails then it is not
    # a match. If no match found then fail with error.
       
    if time_required:
        df_metadata['datetime_format'] = identify_datetime_format(df, df_metadata['datetime_column'], datetime_formats)


    return df, df_metadata


def check_variables_in_dataset(ds, variables):
    """
    Checks if the specified variables are available in the dataset.
    
    Parameters:
    - ds: xarray dataset object.
    - variables: list of strings, names of the variables to check.
    """
    missing_vars = [var for var in variables if var not in ds.variables]
    if missing_vars:
        available_vars = list(ds.variables)
        print(f"Error: The following variables are not available in the dataset: {', '.join(missing_vars)}")
        print("Available variables are:")
        for var in available_vars:
            print(var)
        raise ValueError("One or more requested variables are not available in the dataset.")

def read_lat_lon_time_depth_csv(input_csv_path, lat_column='auto', 
    lon_column='auto', depth_column='auto', datetime_column='auto', datetime_format='auto'):
    """
    This function reads in a CSV of point locations to be analysed. It returns a
    standardised 
    """

def get_ereefs_points_in_time(var_names,
                              lat_lon_time_depth_df,
                              input_file,
                              download_cache = True,
                              verbosity = 1,
                              date_format = 'auto'):
    """
    Extract point data for multiple 4D locations (latitude, longitude, depth, time) from
    eReefs models. Handles both curvilinear model data (raw eReefs data from NCI) and
    regridded regular grid data from AIMS aggregate data service. 
    - var_names: list of variables in the model data (for example: ["temp", "eta", "botz"]). 
        Not all variables need to have the same number of dimensions. 
    - lat_lon_time_depth_df: Data frame containing the coordinates of the data to be extracted.
        
    """

    print(f'Variables to download ({variables}) ')
    path_without_extension, _ = os.path.splitext(output_csv_path)
    data_cache_path = f'{path_without_extension}_cache.csv'
    # Setup or read data cache with specified variables
    print(f'Connecting to OpenDAP: {opendap_url}')
    # Open the dataset
    ds = xr.open_dataset(opendap_url)
    df, df_metadata = setup_data_cache(input_csv_path, data_cache_path, ds, variables)
    
    # Assuming df is your data cache DataFrame
    

    if df['vars_download'].all():
        print(f"Download process is already complete. Output data {output_csv_path}")
        return
 
    print("Download process is incomplete. ")
  
    print(f'Note: This script is restartable due to {data_cache_path}. It will resume processing and download on restart.')
        
    # Check that all the locations are valid
    df = precalculate_indices(df, df_metadata, ds, data_cache_path)
    
    download_data_with_precalculated_indices(df, df_metadata, ds, data_cache_path, download_rate_limit_delay)
    
    convert_cache_to_final_csv(df, output_csv_path)

def main():
    get_ereefs_points_in_time
# Run the script
if __name__ == "__main__":
    main()
