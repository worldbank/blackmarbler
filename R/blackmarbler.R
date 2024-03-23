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
#' roi_sf <- gadm(country = "GHA", level=1, path = tempdir()) %>% st_as_sf()
#'
#' # Daily data: raster for October 3, 2021
#' ken_20210205_r <- bm_extract(roi_sf = roi_sf,
#'                             product_id = "VNP46A2",
#'                             date = "2021-10-03",
#'                             bearer = bearer)
#'
#' # Monthly data: raster for March 2021
#' ken_202103_r <- bm_extract(roi_sf = roi_sf,
#'                           product_id = "VNP46A3",
#'                           date = "2021-03-01",
#'                           bearer = bearer)
#'
#' # Annual data: raster for 2021
#' ken_2021_r <- bm_extract(roi_sf = roi_sf,
#'                         product_id = "VNP46A4",
#'                         date = 2021,
#'                         bearer = bearer)
#'}
#'
#' @export
bm_extract <- function(roi_sf,
                       product_id,
                       date,
                       bearer,
                       aggregation_fun = c("mean"),
                       add_n_pixels = TRUE,
                       variable = NULL,
                       quality_flag_rm = NULL,
                       check_all_tiles_exist = TRUE,
                       interpol_na = FALSE,
                       output_location_type = "memory", # memory, file
                       file_dir = NULL,
                       file_prefix = NULL,
                       file_skip_if_exists = TRUE,
                       quiet = FALSE,
                       ...){

  # Errors & Warnings ----------------------------------------------------------
  if( (interpol_na == T) & (length(date) == 1) ){
    stop("If interpol_na = TRUE, then must have more than one date")
  }

  if( (interpol_na == T) & (output_location_type == "file") ){
    interpol_na <- F
    warning("interpol_na ignored. Interpolation only occurs when output_location_type = 'memory'")
  }

  # Assign interpolation variables ---------------------------------------------
  if(interpol_na == T){
    if(!exists("method")) method <- "linear"
    if(!exists("rule"))   rule   <- 1
    if(!exists("f"))      f      <- 0
    if(!exists("ties"))   ties   <- mean
    if(!exists("z"))      z      <- NULL
    if(!exists("NArule")) NArule <- 1
  }

  # Define Tempdir -------------------------------------------------------------
  temp_main_dir = tempdir()

  current_time_millis = as.character(as.numeric(Sys.time())*1000) %>%
    str_replace_all("[:punct:]", "")
  temp_dir = file.path(temp_main_dir, paste0("bm_raster_temp_", current_time_millis))

  dir.create(temp_dir, showWarnings = F)

  # NTL Variable ---------------------------------------------------------------
  variable <- define_variable(variable, product_id)

  if(interpol_na == T){

    #### Create raster
    bm_r <- bm_raster(roi_sf = roi_sf,
                      product_id = product_id,
                      date = date,
                      bearer = bearer,
                      variable = variable,
                      quality_flag_rm = quality_flag_rm,
                      check_all_tiles_exist = check_all_tiles_exist,
                      interpol_na = F,
                      quiet = quiet,
                      temp_dir = temp_dir)

    bm_r <- raster::approxNA(bm_r,
                             method = method,
                             rule   = rule,
                             f      = f,
                             ties   = ties,
                             z      = z,
                             NArule = NArule)

    #### Extract
    roi_df <- roi_sf %>% st_drop_geometry()
    roi_df$date <- NULL

    n_obs_df <- exact_extract(bm_r, roi_sf, count_n_obs, progress = !quiet) %>%
      bind_cols(roi_df) %>%
      tidyr::pivot_longer(cols = -c(names(roi_df)),
                          names_to = c(".value", "date"),
                          names_sep = "\\.t") %>%
      dplyr::mutate(prop_non_na_pixels = .data$n_non_na_pixels / .data$n_pixels)

    ntl_df <- exact_extract(bm_r, roi_sf, aggregation_fun, progress = !quiet) %>%
      tidyr::pivot_longer(cols = everything(),
                          names_to = c(".value", "date"),
                          names_sep = "\\.t")

    names(ntl_df)[names(ntl_df) != "date"] <-
      paste0("ntl_", names(ntl_df)[names(ntl_df) != "date"])

    ntl_df$date <- NULL
    r <- bind_cols(n_obs_df, ntl_df)
    #r <- ntl_df %>%
    #  left_join(n_obs_df, by = "date")

    # Apply through each date, extract, then append
  } else{

    # Download data --------------------------------------------------------------
    r_list <- lapply(date, function(date_i){

      out <- tryCatch(
        {

          #### Make name for raster based on date
          date_name_i <- define_date_name(date_i, product_id)

          #### If save to file
          if(output_location_type == "file"){

            out_name <- paste0(file_prefix, product_id, "_", date_name_i, ".Rds")
            out_path <- file.path(file_dir, out_name)

            make_raster <- TRUE
            if(file_skip_if_exists & file.exists(out_path)) make_raster <- FALSE

            if(make_raster){

              #### Make raster
              r <- bm_raster_i(roi_sf = roi_sf,
                               product_id = product_id,
                               date = date_i,
                               bearer = bearer,
                               variable = variable,
                               quality_flag_rm = quality_flag_rm,
                               check_all_tiles_exist = check_all_tiles_exist,
                               quiet = quiet,
                               temp_dir = temp_dir)
              names(r) <- date_name_i

              #### Extract
              r_agg <- exact_extract(x = r, y = roi_sf, fun = aggregation_fun,
                                     progress = !quiet)
              roi_df <- roi_sf
              roi_df$geometry <- NULL

              if(length(aggregation_fun) > 1){
                names(r_agg) <- paste0("ntl_", names(r_agg))
                r_agg <- bind_cols(r_agg, roi_df)
              } else{
                roi_df[[paste0("ntl_", aggregation_fun)]] <- r_agg
                r_agg <- roi_df
              }

              if(add_n_pixels){

                r_n_obs <- exact_extract(r_out, roi_sf, function(values, coverage_fraction)
                  sum(!is.na(values)),
                  progress = !quiet)

                r_n_obs_poss <- exact_extract(r_out, roi_sf, function(values, coverage_fraction)
                  length(values),
                  progress = !quiet)

                r_agg$n_pixels           <- r_n_obs_poss
                r_agg$n_non_na_pixels    <- r_n_obs
                r_agg$prop_non_na_pixels <- r_agg$n_non_na_pixels / r_agg$n_pixels
              }

              r_agg$date <- date_i

              #### Export
              saveRDS(r_agg, out_path)

            } else{
              warning(paste0('"', out_path, '" already exists; skipping.\n'))
            }

            r_out <- NULL # Saving as file, so output from function should be NULL

          } else{
            r_out <- bm_raster_i(roi_sf = roi_sf,
                                 product_id = product_id,
                                 date = date_i,
                                 bearer = bearer,
                                 variable = variable,
                                 quality_flag_rm = quality_flag_rm,
                                 check_all_tiles_exist = check_all_tiles_exist,
                                 quiet = quiet,
                                 temp_dir = temp_dir)
            names(r_out) <- date_name_i

            if(add_n_pixels){

              r_n_obs <- exact_extract(r_out, roi_sf, function(values, coverage_fraction)
                sum(!is.na(values)),
                progress = !quiet)

              r_n_obs_poss <- exact_extract(r_out, roi_sf, function(values, coverage_fraction)
                length(values),
                progress = !quiet)

              roi_sf$n_pixels           <- r_n_obs_poss
              roi_sf$n_non_na_pixels    <- r_n_obs
              roi_sf$prop_non_na_pixels <- roi_sf$n_non_na_pixels / roi_sf$n_pixels
            }

            r_out <- exact_extract(x = r_out, y = roi_sf, fun = aggregation_fun,
                                   progress = !quiet)

            roi_df <- roi_sf
            roi_df$geometry <- NULL

            if(length(aggregation_fun) > 1){
              names(r_out) <- paste0("ntl_", names(r_out))
              r_out <- bind_cols(r_out, roi_df)
            } else{

              roi_df[[paste0("ntl_", aggregation_fun)]] <- r_out
              r_out <- roi_df
            }

            r_out$date <- date_i
          }

          return(r_out)

        },
        error=function(e) {
          return(NULL)
        }
      )

    })

    # Clean output ---------------------------------------------------------------
    # Remove NULLs
    r_list <- r_list[!sapply(r_list,is.null)]

    r <- r_list %>%
      bind_rows()

  }

  unlink(temp_dir, recursive = T)
  return(r)
}

#' Make Black Marble Raster
#'
#' Make a raster of nighttime lights from [NASA Black Marble data](https://blackmarble.gsfc.nasa.gov/)

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
#' @param interpol_na When data for more than one date is downloaded, whether to interpolate `NA` values using the `raster::approxNA` function. Additional arguments for the `raster::approxNA` function can also be passed into `bm_raster` (eg, `method`, `rule`, `f`, `ties`, `z`, `NA_rule`). (Default: `FALSE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a raster in R. If `file`, function exports a `.tif` file and returns `NULL`.
#' For `output_location_type = file`:
#' @param file_dir The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].tif`
#' @param file_skip_if_exists Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
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
#' roi_sf <- gadm(country = "GHA", level=0, path = tempdir()) %>% st_as_sf()
#'
#' # Daily data: raster for October 3, 2021
#' ken_20210205_r <- bm_raster(roi_sf = roi_sf,
#'                             product_id = "VNP46A2",
#'                             date = "2021-10-03",
#'                             bearer = bearer)
#'
#' # Monthly data: raster for March 2021
#' ken_202103_r <- bm_raster(roi_sf = roi_sf,
#'                           product_id = "VNP46A3",
#'                           date = "2021-03-01",
#'                           bearer = bearer)
#'
#' # Annual data: raster for 2021
#' ken_2021_r <- bm_raster(roi_sf = roi_sf,
#'                         product_id = "VNP46A4",
#'                         date = 2021,
#'                         bearer = bearer)
#'}
#'
#' @export
#'
#' @import readr
#' @import hdf5r
#' @import dplyr
#' @import sf
#' @import exactextractr
#' @import stringr
#' @import httr
#' @import lubridate
#' @rawNamespace import(tidyr, except = c(extract))
#' @rawNamespace import(purrr, except = c(flatten_df, values))
#' @rawNamespace import(raster, except = c(union, select, intersect, origin, tail, head, values))
#'
# @rawNamespace import(utils, except = c(stack, unstack))
bm_raster <- function(roi_sf,
                      product_id,
                      date,
                      bearer,
                      variable = NULL,
                      quality_flag_rm = NULL,
                      check_all_tiles_exist = TRUE,
                      interpol_na = FALSE,
                      output_location_type = "memory", # memory, file
                      file_dir = NULL,
                      file_prefix = NULL,
                      file_skip_if_exists = TRUE,
                      quiet = FALSE,
                      ...){

  # Errors & Warnings ----------------------------------------------------------
  if( (interpol_na == T) & (length(date) == 1) ){
    stop("If interpol_na = TRUE, then must have more than one date")
  }

  if( (interpol_na == T) & (output_location_type == "file") ){
    interpol_na <- F
    stop("interpol_na ignored. Interpolation only occurs when output_location_type = 'memory'")
  }

  # Assign interpolation variables ---------------------------------------------
  if(interpol_na == T){
    if(!exists("method")) method <- "linear"
    if(!exists("rule"))   rule   <- 1
    if(!exists("f"))      f      <- 0
    if(!exists("ties"))   ties   <- mean
    if(!exists("z"))      z      <- NULL
    if(!exists("NArule")) NArule <- 1
  }

  # Define Tempdir -------------------------------------------------------------
  temp_main_dir = tempdir()

  current_time_millis = as.character(as.numeric(Sys.time())*1000) %>%
    str_replace_all("[:punct:]", "")
  temp_dir = file.path(temp_main_dir, paste0("bm_raster_temp_", current_time_millis))

  dir.create(temp_dir, showWarnings = F)

  # NTL Variable ---------------------------------------------------------------
  variable <- define_variable(variable, product_id)

  # Download data --------------------------------------------------------------
  r_list <- lapply(date, function(date_i){

    out <- tryCatch(
      {

        #### Make name for raster based on date
        date_name_i <- define_date_name(date_i, product_id)

        #### If save as tif format
        if(output_location_type == "file"){
          out_name <- paste0(file_prefix, product_id, "_", date_name_i, ".tif")
          out_path <- file.path(file_dir, out_name)

          make_raster <- TRUE
          if(file_skip_if_exists & file.exists(out_path)) make_raster <- FALSE

          if(make_raster){

            r <- bm_raster_i(roi_sf = roi_sf,
                             product_id = product_id,
                             date = date_i,
                             bearer = bearer,
                             variable = variable,
                             quality_flag_rm = quality_flag_rm,
                             check_all_tiles_exist = check_all_tiles_exist,
                             quiet = quiet,
                             temp_dir = temp_dir)
            names(r) <- date_name_i

            writeRaster(r, out_path)

          } else{
            warning(paste0('"', out_path, '" already exists; skipping.\n'))
          }

          r_out <- NULL # Saving as tif file, so output from function should be NULL

        } else{
          r_out <- bm_raster_i(roi_sf = roi_sf,
                               product_id = product_id,
                               date = date_i,
                               bearer = bearer,
                               variable = variable,
                               quality_flag_rm = quality_flag_rm,
                               check_all_tiles_exist = check_all_tiles_exist,
                               quiet = quiet,
                               temp_dir = temp_dir)
          names(r_out) <- date_name_i

        }

        return(r_out)

      },
      error=function(e) {
        return(NULL)
      }
    )

  })

  # Clean output ---------------------------------------------------------------
  # Remove NULLs
  r_list <- r_list[!sapply(r_list,is.null)]

  if(length(r_list) == 1){
    r <- r_list[[1]]
  } else if (length(r_list) > 1){
    r <- raster::stack(r_list)
  } else{
    r <- NULL
  }

  # Interpolate ----------------------------------------------------------------
  if(interpol_na %in% T){
    r <- raster::approxNA(r,
                          method = method,
                          rule   = rule,
                          f      = f,
                          ties   = ties,
                          z      = z,
                          NArule = NArule)
  }

  unlink(temp_dir, recursive = T)

  return(r)
}
