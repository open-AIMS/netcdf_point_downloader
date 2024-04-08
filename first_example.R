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