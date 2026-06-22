# Extract and Aggregate Black Marble Data

Extract and aggregate nighttime lights data from [NASA Black Marble
data](https://blackmarble.gsfc.nasa.gov/)

## Usage

``` r
bm_extract(
  roi_sf,
  product_id,
  date,
  bearer,
  aggregation_fun = c("sum"),
  add_n_pixels = TRUE,
  variable = NULL,
  quality_flag_rm = NULL,
  check_all_tiles_exist = TRUE,
  interpol_na = FALSE,
  output_location_type = "memory",
  file_dir = NULL,
  file_prefix = NULL,
  file_skip_if_exists = TRUE,
  file_return_null = FALSE,
  h5_dir = NULL,
  download_method = "httr",
  quiet = FALSE,
  ...
)
```

## Arguments

- roi_sf:

  Region of interest; sf polygon. Must be in the [WGS 84
  (epsg:4326)](https://epsg.io/4326) coordinate reference system.

- product_id:

  One of the following:

  - `"VNP46A1"`: Daily (raw)

  - `"VNP46A2"`: Daily (corrected)

  - `"VNP46A3"`: Monthly

  - `"VNP46A4"`: Annual

- date:

  Date of raster data. Entering one date will produce a `SpatRaster`
  object. Entering multiple dates will produce a `SpatRaster` object
  with multiple bands; one band per date.

  - For `product_id`s `"VNP46A1"` and `"VNP46A2"`, a date (eg,
    `"2021-10-03"`).

  - For `product_id` `"VNP46A3"`, a date or year-month (e.g.,
    `"2021-10-01"`, where the day will be ignored, or `"2021-10"`).

  - For `product_id` `"VNP46A4"`, year or date (e.g., `"2021-10-01"`,
    where the month and day will be ignored, or `2021`).

- bearer:

  NASA bearer token. For instructions on how to create a token, see
  [here](https://github.com/worldbank/blackmarbler#bearer-token-).

- aggregation_fun:

  Function used to aggregate nighttime lights data to polygons; this
  values is passed to the `fun` argument in
  [exactextractr::exact_extract](https://github.com/isciences/exactextractr)
  (Default: `sum`).

- add_n_pixels:

  Whether to add a variable indicating the number of nighttime light
  pixels used to compute nighttime lights statistics (eg, number of
  pixels used to compute average of nighttime lights). When `TRUE`, it
  adds three values: `n_non_na_pixels` (the number of non-`NA` pixels
  used for computing nighttime light statistics); `n_pixels` (the total
  number of pixels); and `prop_non_na_pixels` the proportion of the two.
  (Default: `TRUE`).

- variable:

  Variable to used to create raster (default: `NULL`). If `NULL`, uses
  the following default variables:

  - For `product_id` `:VNP46A1"`, uses `DNB_At_Sensor_Radiance_500m`.

  - For `product_id` `"VNP46A2"`, uses
    `Gap_Filled_DNB_BRDF-Corrected_NTL`.

  - For `product_id`s `"VNP46A3"` and `"VNP46A4"`, uses
    `NearNadir_Composite_Snow_Free`. To see all variable choices, set
    `variable = ""` (this will create an error message that lists all
    valid variables). For additional information on variable choices,
    see
    [here](https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.2_April_2021.pdf);
    for `VNP46A1`, see Table 3; for `VNP46A2` see Table 6; for `VNP46A3`
    and `VNP46A4`, see Table 9.

- quality_flag_rm:

  Quality flag values to use to set values to `NA`. Each pixel has a
  quality flag value, where low quality values can be removed. Values
  are set to `NA` for each value in the `quality_flag_rm` vector. Note
  that `quality_flag_rm` does not apply for `VNP46A1`. (Default:
  `NULL`).

  For `VNP46A2` (daily data):

  - `0`: High-quality, Persistent nighttime lights

  - `1`: High-quality, Ephemeral nighttime Lights

  - `2`: Poor-quality, Outlier, potential cloud contamination, or other
    issues

  For `VNP46A3` and `VNP46A4` (monthly and annual data):

  - `0`: Good-quality, The number of observations used for the composite
    is larger than 3

  - `1`: Poor-quality, The number of observations used for the composite
    is less than or equal to 3

  - `2`: Gap filled NTL based on historical data

- check_all_tiles_exist:

  Check whether all Black Marble nighttime light tiles exist for the
  region of interest. Sometimes not all tiles are available, so the full
  region of interest may not be covered. If `TRUE`, skips cases where
  not all tiles are available. (Default: `TRUE`).

- interpol_na:

  When data for more than one date is downloaded, whether to interpolate
  `NA` values in rasters using the
  [`terra::approximate`](https://rspatial.github.io/terra/reference/approximate.html)
  function. Additional arguments for the
  [`terra::approximate`](https://rspatial.github.io/terra/reference/approximate.html)
  function can also be passed into `bm_extract` (eg, `method`, `rule`,
  `f`, `ties`, `z`, `NA_rule`). (Default: `FALSE`).

- output_location_type:

  Where to produce output; either `memory` or `file`. If `memory`,
  functions returns a dataframe in R. If `file`, function exports a
  `.csv` file and returns `NULL`.

- file_dir:

  (If `output_location_type = file`). The directory where data should be
  exported (default: `NULL`, so the working directory will be used)

- file_prefix:

  (If `output_location_type = file`). Prefix to add to the file to be
  saved. The file will be saved as the following:
  `[file_prefix][product_id]_t[date].csv`

- file_skip_if_exists:

  (If `output_location_type = file`). Whether the function should first
  check wither the file already exists, and to skip downloading or
  extracting data if the data for that date if the file already exists
  (default: `TRUE`).

- file_return_null:

  Whether to return `NULL` instead of a `dataframe`. When
  `output_location_type = 'file'`, the function will export data to the
  `file_dir` directory. When `file_return_null = FALSE`, the function
  will also return a `dataframe` of the queried data—so the data is
  available in R memory. Setting `file_return_null = TRUE`, data will be
  saved to `file_dir` but no data will be returned by the function to R
  memory (default: `FALSE`).

- h5_dir:

  Black Marble data are originally downloaded as `h5` files. If
  `h5_dir = NULL`, the function downloads to a temporary directory then
  deletes the directory. If `h5_dir` is set to a path, `h5` files are
  saved to that directory and not deleted. The function will then check
  if the needed `h5` file already exists in the directory; if it exists,
  the function will not re-download the `h5` file.

- download_method:

  Method to download data (h5 files) from NASA LAADS Archive: "`httr`"
  or "`wget`". If `httr`, uses the `httr2` R package to download data.
  If `wget`, uses the `wget` command line tool. `httr` is fully
  integrated in R, while `wget` requires the `wget` system command.
  `wget` can be more efficient and can help avoid network issues.
  (Default: `"httr"`).

- quiet:

  Suppress output that show downloading progress and other messages.
  (Default: `FALSE`).

- ...:

  Additional arguments for
  [`terra::approximate`](https://rspatial.github.io/terra/reference/approximate.html),
  if `interpol_na = TRUE`

## Value

Raster

## Author

Robert Marty <rmarty@worldbank.org>

## Examples

``` r
if (FALSE) { # \dontrun{
# Define bearer token
bearer <- "BEARER-TOKEN-HERE"

# sf polygon of Ghana
library(geodata)
roi_sf <- gadm(country = "GHA", level=1, path = tempdir()) %>% st_as_sf()

# Daily data: raster for October 3, 2021
ken_20210205_r <- bm_extract(roi_sf = roi_sf,
                            product_id = "VNP46A2",
                            date = "2021-10-03",
                            bearer = bearer)

# Monthly data: raster for March 2021
ken_202103_r <- bm_extract(roi_sf = roi_sf,
                          product_id = "VNP46A3",
                          date = "2021-03-01",
                          bearer = bearer)

# Annual data: raster for 2021
ken_2021_r <- bm_extract(roi_sf = roi_sf,
                        product_id = "VNP46A4",
                        date = 2021,
                        bearer = bearer)
} # }
```
