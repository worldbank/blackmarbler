#' Extract and Aggregate Black Marble Data
#'
#' Extract and aggregate nighttime lights data from [NASA Black Marble data](https://blackmarble.gsfc.nasa.gov/)

#' @param roi_sf Region of interest; sf polygon. Must be in the [WGS 84 (epsg:4326)](https://epsg.io/4326) coordinate reference system.
#' @param product_id One of the following:
#' * `"VNP46A1"`: Daily (raw)
#' * `"VNP46A2"`: Daily (corrected)
#' * `"VNP46A3"`: Monthly
#' * `"VNP46A4"`: Annual
#' @param date Date of raster data. Entering one date will produce a raster. Entering multiple dates will produce a raster stack.
#' * For `product_id`s `"VNP46A1"` and `"VNP46A2"`, a date (eg, `"2021-10-03"`).
#' * For `product_id` `"VNP46A3"`, a date or year-month (e.g., `"2021-10-01"`, where the day will be ignored, or `"2021-10"`).
#' * For `product_id` `"VNP46A4"`, year or date  (e.g., `"2021-10-01"`, where the month and day will be ignored, or `2021`).
#' @param bearer NASA bearer token. For instructions on how to create a token, see [here](https://github.com/worldbank/blackmarbler#bearer-token-).
#' @param aggregation_fun Function used to aggregate nighttime lights data to polygons; this values is passed to the `fun` argument in [exactextractr::exact_extract](https://github.com/isciences/exactextractr) (Default: `mean`).
#' @param add_n_pixels Whether to add a variable indicating the number of nighttime light pixels used to compute nighttime lights statistics (eg, number of pixels used to compute average of nighttime lights). When `TRUE`, it adds three values: `n_non_na_pixels` (the number of non-`NA` pixels used for computing nighttime light statistics); `n_pixels` (the total number of pixels); and `prop_non_na_pixels` the proportion of the two. (Default: `TRUE`).
#' @param variable Variable to used to create raster (default: `NULL`). If `NULL`, uses the following default variables:
#' * For `product_id` `:VNP46A1"`, uses `DNB_At_Sensor_Radiance_500m`.
#' * For `product_id` `"VNP46A2"`, uses `Gap_Filled_DNB_BRDF-Corrected_NTL`.
#' * For `product_id`s `"VNP46A3"` and `"VNP46A4"`, uses `NearNadir_Composite_Snow_Free`.
#' For information on other variable choices, see [here](https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.2_April_2021.pdf); for `VNP46A1`, see Table 3; for `VNP46A2` see Table 6; for `VNP46A3` and `VNP46A4`, see Table 9.
#' @param quality_flag_rm Quality flag values to use to set values to `NA`. Each pixel has a quality flag value, where low quality values can be removed. Values are set to `NA` for each value in ther `quality_flag_rm` vector. (Default: `NULL`).
#'
#'
#' For `VNP46A1` and `VNP46A2` (daily data):
#' - `0`: High-quality, Persistent nighttime lights
#' - `1`: High-quality, Ephemeral nighttime Lights
#' - `2`: Poor-quality, Outlier, potential cloud contamination, or other issues
#'
#'
#' For `VNP46A3` and `VNP46A4` (monthly and annual data):
#' - `0`: Good-quality, The number of observations used for the composite is larger than 3
#' - `1`: Poor-quality, The number of observations used for the composite is less than or equal to 3
#' - `2`: Gap filled NTL based on historical data
#' @param check_all_tiles_exist Check whether all Black Marble nighttime light tiles exist for the region of interest. Sometimes not all tiles are available, so the full region of interest may not be covered. If `TRUE`, skips cases where not all tiles are available. (Default: `TRUE`).
#' @param interpol_na When data for more than one date is downloaded, whether to interpolate `NA` values in rasters using the `raster::approxNA` function. Additional arguments for the `raster::approxNA` function can also be passed into `bm_extract` (eg, `method`, `rule`, `f`, `ties`, `z`, `NA_rule`). (Default: `FALSE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a dataframe in R. If `file`, function exports a `.csv` file and returns `NULL`.
#' @param file_dir (If `output_location_type = file`). The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix (If `output_location_type = file`). Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].csv`
#' @param file_skip_if_exists (If `output_location_type = file`). Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
#'
#' @param ... Additional arguments for `raster::approxNA`, if `interpol_na = TRUE`
#'
#' @return Raster
#'
#' @examples
#' \dontrun{
#' # Define bearer token
#' bearer <- "BEARER-TOKEN-HERE"
#'
#' # sf polygon of Ghana
#' library(geodata)
#' roi_sf <- gadm(country = "GHA", level = 1, path = tempdir()) %>% st_as_sf()
#'
#' # Daily data: raster for October 3, 2021
#' ken_20210205_r <- bm_extract(
#'   roi_sf = roi_sf,
#'   product_id = "VNP46A2",
#'   date = "2021-10-03",
#'   bearer = bearer
#' )
#'
#' # Monthly data: raster for March 2021
#' ken_202103_r <- bm_extract(
#'   roi_sf = roi_sf,
#'   product_id = "VNP46A3",
#'   date = "2021-03-01",
#'   bearer = bearer
#' )
#'
#' # Annual data: raster for 2021
#' ken_2021_r <- bm_extract(
#'   roi_sf = roi_sf,
#'   product_id = "VNP46A4",
#'   date = 2021,
#'   bearer = bearer
#' )
#' }
#'
#' @export
bm_raster <- function(roi_sf,
                      product_id,
                      date,
                      bearer,
                      variable = NULL,
               quality_flags_to_remove = NULL,
                      check_all_tiles_exist = TRUE,
                      interpol_na = FALSE,
                      output_location_type = "memory", # memory, file
                      file_dir = NULL,
                      file_prefix = NULL,
                      file_skip_if_exists = TRUE,
                      quiet = FALSE       ,
                      ...) {
  # Error Checks ---------------------------------------------------------------
  if (interpol_na && length(date) == 1) {
    stop("If interpol_na = TRUE, then must have more than one date")
  }

  if (interpol_na && output_location_type == "file") {
    interpol_na <- FALSE
    warning("interpol_na ignored. Interpolation only occurs when output_location_type = 'memory'")
  }

  # Assign Interpolation Variables ---------------------------------------------
  if (interpol_na) {
    method <- "linear"
    rule <- 1
    f <- 0
    ties <- mean
    z <- NULL
    NArule <- 1
  }

  # Define Temporary Directory -------------------------------------------------
  temp_dir <- tempfile(pattern = paste0("bm_raster_temp_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(temp_dir, showWarnings = FALSE)

  # Define NTL Variable --------------------------------------------------------
  blackmarble_variable <- define_blackmarble_variable(variable, product_id)

  # Download and Process Rasters -----------------------------------------------
  r_list <- lapply(date, function(date_i) {
    out <- tryCatch(
      {
        date_name_i <- define_raster_name(date_i, product_id)

        if (output_location_type == "file") {
          out_name <- paste0(file_prefix, product_id, "_", date_name_i, ".tif")
          out_path <- file.path(file_dir, out_name)

          if (file_skip_if_exists && file.exists(out_path)) {
            warning(paste0('"', out_path, '" already exists; skipping.\n'))
            return(NULL)
          }
        }

        r <- retrieve_and_process_nightlight_data(
          roi_sf = roi_sf,
          product_id = product_id,
          date = date_i,
          bearer = bearer,
          variable = blackmarble_variable,
          quality_flags_to_remove = quality_flags_to_remove,
          check_all_tiles_exist = check_all_tiles_exist,
          quiet = quiet,
          temp_dir = temp_dir
        )

        if (output_location_type == "file") {
          terra::writeRaster(r, out_path)
          return(NULL) # Saving as tif file, so output from function should be NULL
        } else {
          names(r) <- date_name_i
          return(r)
        }
      },
      error = function(e) {
        # add error message of error
        return(NULL)
      }
    )
  })

  # Clean output ---------------------------------------------------------------
  r_list <- r_list[!sapply(r_list, is.null)]

  r <- if (length(r_list) == 1) r_list[[1]] else terra::rast(r_list)

  # Interpolate ----------------------------------------------------------------
  if (interpol_na) {
    r <-
      terra::approximate(
        r,
        method = method,
        rule = rule,
        f = f,
        ties = ties,
        z = z,
        NArule = NArule
      )
  }

  unlink(temp_dir, recursive = TRUE)

  return(r)
}


#' Extracts raster data for specified region of interest
#'
#' @param roi_sf Spatial object representing the region of interest.
#' @param product_id Identifier for the product.
#' @param date Dates for which raster data is to be extracted.
#' @param bearer Authentication token for accessing the data.
#' @param aggregation_fun Aggregation function to summarize raster values (default is "mean").
#' @param add_n_pixels Logical indicating whether to add pixel count information (default is TRUE).
#' @param variable Variable to extract from the product (default is NULL).
#' @param quality_flags_to_remove Quality flag removal options (default is NULL).
#' @param check_all_tiles_exist Logical indicating whether to check if all tiles exist (default is TRUE).
#' @param interpol_na Logical indicating whether to interpolate missing values (default is FALSE).
#' @param output_location_type Output location type ("memory" or "file", default is "memory").
#' @param file_dir Directory for saving output files (used if output_location_type is "file", default is NULL).
#' @param file_prefix Prefix for output file names (default is NULL).
#' @param file_skip_if_exists Logical indicating whether to skip file creation if output file already exists (default is TRUE).
#' @param quiet Logical indicating whether to suppress progress messages (default is FALSE).
#' @return A data frame containing extracted raster data.
#' @export
bm_extract <- function(roi_sf,
                       product_id,
                       date,
                       bearer,
                       aggregation_fun = c("mean"),
                       add_n_pixels = TRUE,
                       variable = NULL,
                       quality_flags_to_remove = NULL,
                       check_all_tiles_exist = TRUE,
                       interpol_na = FALSE,
                       output_location_type = "memory",
                       file_dir = NULL,
                       file_prefix = NULL,
                       file_skip_if_exists = TRUE,
                       quiet = FALSE,
                       ...) {
  # Errors & Warnings ----------------------------------------------------------
  if (interpol_na & length(date) == 1) {
    stop("If interpol_na = TRUE, then must have more than one date")
  }

  if (interpol_na & output_location_type == "file") {
    warning("interpol_na ignored. Interpolation only occurs when output_location_type = 'memory'")
    interpol_na <- FALSE
  }

  # Define Tempdir -------------------------------------------------------------
  temp_dir <- tempdir()

  # NTL Variable ---------------------------------------------------------------
  variable <- define_blackmarble_variable(variable, product_id)

  if (interpol_na) {
    # Create raster
    bm_r <- bm_raster(
      roi_sf = roi_sf,
      product_id = product_id,
      date = date,
      bearer = bearer,
      variable = variable,
      quality_flags_to_remove = quality_flags_to_remove,
      check_all_tiles_exist = check_all_tiles_exist,
      interpol_na = FALSE,
      quiet = quiet,
      temp_dir = temp_dir
    )

    # Interpolate
    bm_r <- terra::approximate(bm_r,
      method = "linear",
      rule = 1,
      f = 0,
      ties = mean,
      z = NULL,
      NArule = 1
    )

    # Extract
    n_obs_df <- extract_and_process(bm_r, roi_sf, count_n_obs, quiet)
    ntl_df <- extract_and_process(bm_r, roi_sf, aggregation_fun, quiet)

    r <- bind_extracted_data(n_obs_df, ntl_df)
  } else {
    # Download data --------------------------------------------------------------
    r_list <- lapply(date, function(date_i) {
      tryCatch(
        {
          # Make name for raster based on date
          date_name_i <- define_raster_name(date_i, product_id)

          # If save to file
          if (output_location_type == "file") {
            out_name <- paste0(file_prefix, product_id, "_", date_name_i, ".Rds")
            out_path <- file.path(file_dir, out_name)

            if (file_skip_if_exists && file.exists(out_path)) {
              return(NULL)
            }

            r_agg <- extract_and_process_i(
              roi_sf, product_id, date_i, bearer, variable,
              quality_flags_to_remove, check_all_tiles_exist, quiet, temp_dir
            )

            export_extracted_data(r_agg, out_path)

            return(NULL) # Saving as file, so output from function should be NULL
          } else {
            r_out <- extract_and_process_i(
              roi_sf, product_id, date_i, bearer, variable,
              quality_flags_to_remove, check_all_tiles_exist, quiet, temp_dir
            )

            return(r_out)
          }
        },
        error = function(e) {
          return(NULL)
        }
      )
    })

    r <- bind_extracted_data_list(r_list)
  }

  unlink(temp_dir, recursive = TRUE)
  return(r)
}


# PENDING REPEAT DOCUMENTION OF PREVIOUS ----------------------------------
