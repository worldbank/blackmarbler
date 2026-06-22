# Download h5 files

Download h5 files from from [NASA Black Marble
data](https://blackmarble.gsfc.nasa.gov/) using `wget`. The
wget_h5_files() function requires the wget command line tool to be
installed on your system. If you do not have wget installed, please
install it from https://www.gnu.org/software/wget/.

## Usage

``` r
download_h5_files(
  roi_sf = NULL,
  product_id,
  date,
  h5_dir,
  bearer,
  download_method = "httr"
)
```

## Arguments

- roi_sf:

  Region of interest; sf polygon. Must be in the [WGS 84
  (epsg:4326)](https://epsg.io/4326) coordinate reference system. If
  `NULL`, all `h5` files for the inputted date(s) are downloaded.

- product_id:

  One of the following:

  - `"VNP46A1"`: Daily (raw)

  - `"VNP46A2"`: Daily (corrected)

  - `"VNP46A3"`: Monthly

  - `"VNP46A4"`: Annual

- date:

  Date(s) to download `h5` files.

  - For `product_id`s `"VNP46A1"` and `"VNP46A2"`, a date (eg,
    `"2021-10-03"`).

  - For `product_id` `"VNP46A3"`, a date or year-month (e.g.,
    `"2021-10-01"`, where the day will be ignored, or `"2021-10"`).

  - For `product_id` `"VNP46A4"`, year or date (e.g., `"2021-10-01"`,
    where the month and day will be ignored, or `2021`).

- h5_dir:

  Path to download `h5` files to.

- bearer:

  NASA bearer token. For instructions on how to create a token, see
  [here](https://github.com/worldbank/blackmarbler#bearer-token-).

## Value

`NULL`

## Author

Robert Marty <rmarty@worldbank.org>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define bearer token
bearer <- "BEARER-TOKEN-HERE"

# sf polygon of Ghana
library(geodata)
roi_sf <- gadm(country = "GHA", level=0, path = tempdir()) %>% st_as_sf()

# h5 files for Ghana for October 3, 2021
wget_h5_files(roi_sf = roi_sf,
              product_id = "VNP46A2",
              date = "2021-10-03",
              h5_dir = getwd(),        
              bearer = bearer)

# Make raster using h5_files
ken_202103_r <- bm_raster(roi_sf = roi_sf,
                          product_id = "VNP46A3",
                          date = "2021-03-01",
                          bearer = bearer,
                          h5_dir = getwd())

} # }
```
