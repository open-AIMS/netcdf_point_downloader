# NetCDF Point Downloader
This utility library makes downloading data from NetCDF sources more straight forward, particularly data made available via THREDDS data services. This library was originally developed for downloading data from eReefs data services, but has been made flexible enough to download data from other data services. It has been tested on:
- [BARRA2](https://dapds00.nci.org.au/thredds/catalogs/ob53/catalog.html) - Australian regional atmospheric reanalysis - [Metadata](http://www.bom.gov.au/research/publications/researchreports/BRR-067.pdf)
- [eReefs CSIRO NCI raw model data](https://dapds00.nci.org.au/thredds/catalogs/fx3/catalog.html)
- [eReefs AIMS regridded aggregation data](https://thredds.ereefs.aims.gov.au/thredds/catalog/catalog.html)
- [IMOS GHRSST](https://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/SST/ghrsst/L3SM-1d/day/2020/catalog.html) - Satellite Sea Surface Temperature data
# NetCDF Points Downloader

The NetCDF Points Downloader is a R script designed for extracting data from NetCDF sources, such as eReefs, BARRA2 models, and IMOS GHRSST remote sensing satellite sea surface temperature. This tool facilitates the download of specific data points based on location (latitude and longitude), time, and depth, making it invaluable for researchers requiring precise data extraction from large collections of NetCDF files made available via a THREDDs data service.

## Features

- **Flexible Data Extraction**: Allows the extraction of data for specific variables from NetCDF files based on user-defined locations, times, and depths. These can be specified as a data frame or as a CSV.
- **Support for Various Data Sources**: Designed to work with multiple NetCDF data sources, including eReefs (curvilinear raw model data from NCI and regular gridded data from AIMS), BARRA2, and IMOS GHRSST. Other sources may work, but have not been tested.
- **Automatic Resumption**: Capable of auto-saving progress and resuming downloads after interruptions. It might be slow but it will get there. Once the data has been downloaded the data extraction pulls the data from the local copy. The progressive download is saved to an RData file that assists with debugging issues, and the final download is saved as a CSV that allows easy review of the data.
- **OpenDAP data access**: Data can be read directly from a THREDDs data service. Some services make their entire time series available as a single aggregate time series OpenDAP service. You can access all time slices from a single data end point in these cases.
- **Template URL Support**: Handles data sources with multiple end points through template URLs, adjusting for time zone differences. Long time series are often broken into multiple NetCDF files. When made available via THREDDs these each appear as different OpenDAP Data URLs. The set of OpenDAP end points can be described using a template, where the year, month and day of the input times are substituted into the URL to work out the end point URL. For example:
```
https://thredds.aodn.org.au/thredds/dodsC/IMOS/SRS/SST/ghrsst/L3SM-1d/day/%(year)/%(year)%(month)%(day)032000-ABOM-L3S_GHRSST-SSTskin-MultiSensor-1d_day.nc
```
- **Supports regular and curvilinear grid data**: Allows data download from the raw eReefs model data (curvilinear grid), providing access to the highest temporal resolution version of eReefs. It also supports regular rectangular grids used in most data sources including the temporal aggregation (daily, monthly, yearly) eReefs data services provided by AIMS.   

## How It Works

The `netcdf_points_downloader` script is driven from an input data frame (as a data structure or a CSV file) that contains latitude, longitude and optionally datetime and depth attributes. This information is used to download specified variables from a NetCDF data source, typically an OpenDAP data service provided by THREDDS data server. The download data source can be a single end point, or a set specified by a template, where the year, month and day are substituted to determine the URL of the data source for each point being downloaded. This makes this script particularly useful when dealing with data sources lacking an aggregate data service, requiring connections to individual files within a time series.

The script locates the nearest grid cells in the NetCDF data relative to the specified positions and downloads the data for the requested variables. It robustly saves progress, allowing downloads to resume after interruptions. If input data remains unchanged across restarts, the script utilizes the existing data cache to expedite the process. Any modification to the input data triggers a fresh download from the beginning.

To support the use of diverse datasets as inputs, the script automatically matches various potential names for latitude, longitude, datetime, and depth dimensions. For instance, it accepts multiple naming conventions for each dimension, such as 'lat', 'latitude', 'LAT', etc., for latitude. It applies similar automatic matching for names within the NetCDF data sources, with an option to override default matches through the `netcdf_points_downloader_globals` variables if needed.
The following are the default matching names:
- **latitude**: 'lat', 'latitude', 'LAT', 'LATITUDE', 'Y', 'y', 'i', 
- **longitude**: 'long', 'longitude', 'lon','LONG', 'LONGITUDE', 'LON', 'X', 'x', 'j'
- **date time**: 'datetime', 'dateTime', 'date', 'Date', 'time', 'Time'
- **depth**: 'depth', 'zc', 'Depth', 'depth_m', 'k'

The script also handles a variety of datetime formats automatically and can be customized by modifying the netcdf_points_downloader_globals$input_datetime_formats variable to include specific formats.

For time dimensions in download data sources, the script analyzes the 'units' metadata to determine the time unit ('days' or 'seconds') and the base date for time calculations. Notably, datasets may vary in their base date, and some (like BARRA2) may change this mid-series. These changes should be handled automatically.

Correct time zone specification is crucial, especially for data with a fine temporal scale. The script defaults to UTC+10 (Australia/Brisbane time) but can be adjusted. It supports R's time zone specifications. The allowable time zone strings varies across platforms, however most platforms support time zones of the form 'Etc/GMT+n' and 'Etc/GMT-n' (possibly also without prefix 'Etc/'), which assume a fixed offset from UTC, hence no DST (day light savings time). Contrary to some expectations, negative offsets are times ahead of (East of) UTC, positive offsets are times behind (West of) UTC. So UTC+10 should be specified as 'ETC/GMT-10'

The input data can include additional attributes, which are preserved and passed through to the output, provided there are no naming conflicts with the downloaded variables. (Note: I don't think the script detects and reports on these clashes, I suspect it will just overwrite any existing data).

The script extracts the data from the nearest data grid cell. To indicate the proximity of the downloaded data to the requested location, it calculates and includes lat_val, lon_val, time_val, and depth_val attributes in the output.
- **lat_val**: Latitude of the centroid of the grid cell that the data was extracted from.
- **lon_val**: Longitude of the centroid of the grid cell that the data was extracted from.
- **time_val**: Date time of the time slice that the data was extracted from.
- **depth_val**: Value of the depth (in the units of the dataset, typically metres) that the data was extracted from.


### Limitations

- **Speed**: The script is relatively slow, especially for large datasets. From a remote THREDDS service it can download approximately 10,000 - 20,000 points per hour. For extensive time series data, other tools might be more efficient. While the download can be slow it gets progressively saved and so a restart or interruption doesn't cause any problems, and the loading of the data is fast once the download is complete.
- **Requires R Environment**: Execution requires an R environment with necessary libraries (`ncdf4`).
- **Limited NetCDF support**: The NetCDF standard is very flexible and every new dataset that netcdf_points_downloader was tested on raised a new or varied way of expressing the data. It is highly likely that if you are using a new untested dataset that the code may not work. It is unlikely that the curvilinear data support would work on datasets other than eReefs as it required some customisations to correctly link the names of the curvilinear indices (i, j) to the 2D latitude and longitude specifications. There is an override for this link up `netcdf_points_downloader_globals$download_dim_map`, but this has not been tested on other curvilinear data.
- **Download variables might clash**: At the moment variables that are downloaded are joined to the input data frame using the name in the download source. This might result in a clash with other data in the input data frame. If this is a problem it can be worked around by downloading the data using an input data frame with just the location, time, depth attributes, then copy the downloaded variables from the resulting data frame to their final destination.


### Prerequisites

Ensure the `ncdf4` libraries are installed in your R environment.
Copy the `netcdf_points_downloader.R` into your project, or load it via source from where ever you put the file.

### Basic Usage

1. **Set Up**: Define the variables to download, the input data source (CSV file or data frame), and the NetCDF source URL.
2. **Execution**: Call `netcdf_points_downloader` with appropriate parameters.
3. **Output**: The script will output a data frame and/or CSV file, including the original input data plus the downloaded NetCDF data.

### Examples


The following example demonstrates how to use the script to compare weather station data with eReefs GBR1 Hydro 2.0 data:

```r
# First example:
source('netcdf_points_downloader.R')

if (!dir.exists('temp')) {
  dir.create('temp', recursive = TRUE)
}

# Create an example location to download
input_data <- data.frame(lat = c(-18.83162), lon = c(147.6345), time = c('2022-12-15 00:00:00'), depth = c(-4))

input_data_src <- 'temp/davies_ereefs_temp.csv'
write.csv(input_data, input_data_src)

# input_data_src <- input_data  # You can input a data frame or path to a CSV file

var_names <- c('salt')
download_data_src_url <- 'https://dapds00.nci.org.au/thredds/dodsC/fx3/model_data/gbr1_2.0.ncml'
output_csv_path <- 'temp/output.csv'
verbosity <- 2      # Tell us everything that is going on

downloaded_data <- netcdf_points_downloader(
    var_names = var_names,
    input_data_src = input_data_src,
    download_data_src_url_template = download_data_src_url,
    output_csv_path = output_csv_path,
    verbosity = verbosity,
    timezone_input_data_src = 'Etc/GMT-10',
    timezone_download_data_src = 'Etc/GMT-10'
)
print(downloaded_data)

```

This will fetch the specified `var_names` from the NetCDF data source and append them to the input data frame or save to a new CSV file as specified by `output_csv_path`.

There are additonal more extensive examples in the `examples` folder.

# Weather station example

![examples/weather-station-davies-rf/weather-st-davies-vs-ereefs-vs-barra2-wind.png](Plot of the wind speed from the weather station, eReefs and the BARRA2 models.)


## License

This script is provided under and Apache 2 license.




