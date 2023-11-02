# TODO
# 1. Cloud mask: https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.1_July_2020.pdf
#   file:///Users/robmarty/Downloads/Thesis_Zihao_Zheng.pdf

if(F){
  library(purrr)
  library(furrr)
  library(stringr)
  library(rhdf5)
  library(raster)
  library(dplyr)
  library(sf)
  library(lubridate)
  library(exactextractr)
}

#' Black Marble Tile Grid Shapefile
#' A dataset containing black marble grid tiles.
#'
#' @name bm_tiles_sf
#' @docType data
#' @usage data(bm_tiles_sf)
#' @format Sf polygon with 648 features
#' @author Robert Marty \email{rmarty@worldbank.org}
#' @references \url{https://blackmarble.gsfc.nasa.gov/}
#' @keywords datasets
NULL

month_start_day_to_month <- function(x){

  month <- NA

  if(x == "001") month <- "01"

  if(x == "032") month <- "02"

  if(x == "060") month <- "03"
  if(x == "061") month <- "03"

  if(x == "091") month <- "04"
  if(x == "092") month <- "04"

  if(x == "121") month <- "05"
  if(x == "122") month <- "05"

  if(x == "152") month <- "06"
  if(x == "153") month <- "06"

  if(x == "182") month <- "07"
  if(x == "183") month <- "07"

  if(x == "213") month <- "08"
  if(x == "214") month <- "08"

  if(x == "244") month <- "09"
  if(x == "245") month <- "09"

  if(x == "274") month <- "10"
  if(x == "275") month <- "10"

  if(x == "305") month <- "11"
  if(x == "306") month <- "11"

  if(x == "335") month <- "12"
  if(x == "336") month <- "12"

  return(month)
}

month_start_day_to_month <- Vectorize(month_start_day_to_month)

pad2 <- function(x){
  if(nchar(x) == 1) out <- paste0("0", x)
  if(nchar(x) == 2) out <- paste0(x)
  return(out)
}
pad2 <- Vectorize(pad2)

pad3 <- function(x){
  if(nchar(x) == 1) out <- paste0("00", x)
  if(nchar(x) == 2) out <- paste0("0", x)
  if(nchar(x) == 3) out <- paste0(x)
  return(out)
}
pad3 <- Vectorize(pad3)

file_to_raster <- function(f,
                           variable,
                           quality_flag_rm){
  # Converts h5 file to raster.
  # ARGS
  # --f: Filepath to h5 file

  ## Boundaries
  # Only works on later years
  #spInfo <- h5readAttributes(f,"/")

  #xMin<-spInfo$WestBoundingCoord
  #yMin<-spInfo$SouthBoundingCoord
  #yMax<-spInfo$NorthBoundingCoord
  #xMax<-spInfo$EastBoundingCoord

  ## Data
  h5_data <- H5Fopen(f)

  #### Daily
  if(f %>% str_detect("VNP46A1|VNP46A2")){

    tile_i <- f %>% str_extract("h\\d{2}v\\d{2}")

    bm_tiles_sf <- read_sf("https://raw.githubusercontent.com/ramarty/blackmarbler/main/data/blackmarbletiles.geojson")
    grid_i_sf <- bm_tiles_sf[bm_tiles_sf$TileID %in% tile_i,]

    grid_i_sf_box <- grid_i_sf %>%
      st_bbox()

    xMin <- min(grid_i_sf_box$xmin) %>% round()
    yMin <- min(grid_i_sf_box$ymin) %>% round()
    xMax <- max(grid_i_sf_box$xmax) %>% round()
    yMax <- max(grid_i_sf_box$ymax) %>% round()

    if(!(variable %in% names(h5_data$HDFEOS$GRIDS$VNP_Grid_DNB$`Data Fields`))){
      warning(paste0("'", variable, "'",
                     " not a valid variable option. Valid options include:\n",
                     paste(names(h5_data$HDFEOS$GRIDS$VNP_Grid_DNB$`Data Fields`), collapse = "\n")
      ))
    }

    out <- h5_data$HDFEOS$GRIDS$VNP_Grid_DNB$`Data Fields`[[variable]]
    qf  <- h5_data$HDFEOS$GRIDS$VNP_Grid_DNB$`Data Fields`$Mandatory_Quality_Flag

    if(length(quality_flag_rm) > 0){
      for(val in quality_flag_rm){ # out[qf %in% quality_flag_rm] doesn't work, so loop
        out[qf == val] <- NA
      }
    }

    #### Monthly/Annually
  } else{
    lat <- h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`$lat
    lon <- h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`$lon

    if(!(variable %in% names(h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`))){
      warning(paste0("'", variable, "'",
                     " not a valid variable option. Valid options include:\n",
                     paste(names(h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`), collapse = "\n")
      ))
    }

    out <- h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`[[variable]]

    if(length(quality_flag_rm) > 0){

      variable_short <- variable %>%
        str_replace_all("_Num", "") %>%
        str_replace_all("_Std", "")

      qf_name <- paste0(variable_short, "_Quality")
      if(qf_name %in% names(h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`)){

        qf <- h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`[[paste0(variable, "_Quality")]]

        for(val in quality_flag_rm){ # out[qf %in% quality_flag_rm] doesn't work, so loop
          out[qf == val] <- NA
        }
      }

    }

    if(class(out[1,1])[1] != "numeric"){
      out <- matrix(as.numeric(out),    # Convert to numeric matrix
                    ncol = ncol(out))
    }

    xMin <- min(lon) %>% round()
    yMin <- min(lat) %>% round()
    xMax <- max(lon) %>% round()
    yMax <- max(lat) %>% round()

  }

  ## Metadata
  nRows      <- nrow(out)
  nCols      <- ncol(out)
  res        <- nRows
  nodata_val <- NA
  myCrs      <- 4326

  ## Make Raster

  #transpose data to fix flipped row and column order
  #depending upon how your data are formatted you might not have to perform this
  out <- t(out)
  #land_water_mask <- t(land_water_mask)

  #assign data ignore values to NA
  out[out == nodata_val] <- NA

  #turn the out object into a raster
  outr <- raster(out,crs=myCrs)
  #land_water_mask_r <- raster(land_water_mask,crs=myCrs)

  #create extents class
  rasExt <- raster::extent(c(xMin,xMax,yMin,yMax))

  #assign the extents to the raster
  extent(outr) <- rasExt
  #extent(land_water_mask_r) <- rasExt

  #water to 0
  #outr[][outr[] %in% 65535] <- NA # This is a fill value; always exclude

  h5closeAll()

  return(outr)
}

read_bm_csv <- function(year,
                        day,
                        product_id){
  #print(paste0("Reading: ", product_id, "/", year, "/", day))
  df_out <- tryCatch(
    {
      df <- readr::read_csv(paste0("https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/",product_id,"/",year,"/",day,".csv"),
                            show_col_types = F)

      df$year <- year
      df$day <- day

      df
    },
    error = function(e){
      warning(paste0("Error with year: ", year, "; day: ", day))
      data.frame(NULL)
    }
  )

  Sys.sleep(0.1)

  return(df_out)
}

create_dataset_name_df <- function(product_id,
                                   all = TRUE,
                                   years = NULL,
                                   months = NULL,
                                   days = NULL){

  #### Prep dates
  if(product_id %in% c("VNP46A1", "VNP46A2")) months <- NULL
  if(product_id %in% c("VNP46A3"))            days <- NULL
  if(product_id %in% c("VNP46A4")){
    days <- NULL
    months <- NULL
  }

  #### Determine end year
  year_end <- Sys.Date() %>%
    substring(1,4) %>%
    as.numeric()

  #### Make parameter dataframe
  if(product_id %in% c("VNP46A1", "VNP46A2")){
    param_df <- cross_df(list(year = 2012:year_end,
                              day  = pad3(1:366)))
  }

  if(product_id == "VNP46A3"){
    param_df <- cross_df(list(year            = 2012:year_end,
                              day = c("001", "032", "061", "092", "122", "153", "183", "214", "245", "275", "306", "336",
                                      "060", "091", "121", "152", "182", "213", "244", "274", "305", "335")))
  }

  if(product_id == "VNP46A4"){
    param_df <- cross_df(list(year = 2012:year_end,
                              day  = "001"))
  }

  #### Add month if daily or monthly data
  if(product_id %in% c("VNP46A1", "VNP46A2", "VNP46A3")){

    param_df <- param_df %>%
      dplyr::mutate(month = day %>%
                      month_start_day_to_month() %>%
                      as.numeric())

  }

  #### Subset time period
  ## Year
  if(!is.null(years)){
    param_df <- param_df[param_df$year %in% years,]
  }

  ## Month
  if(product_id %in% c("VNP46A1", "VNP46A2", "VNP46A3")){

    if(!is.null(months)){
      param_df <- param_df[as.numeric(param_df$month) %in% as.numeric(months),]
    }

    if(!is.null(days)){
      param_df <- param_df[as.numeric(param_df$day) %in% as.numeric(days),]
    }

  }

  #### Create data
  files_df <- map2_dfr(param_df$year,
                       param_df$day,
                       read_bm_csv,
                       product_id)

  return(files_df)
}

download_raster <- function(file_name,
                            temp_dir,
                            variable,
                            bearer,
                            quality_flag_rm,
                            quiet){

  year       <- file_name %>% substring(10,13)
  day        <- file_name %>% substring(14,16)
  product_id <- file_name %>% substring(1,7)

  wget_command <- paste0('wget -e robots=off -m -np .html,.tmp -nH --cut-dirs=3 ',
                         '"https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/',product_id,'/', year, '/', day, '/', file_name,'"',
                         ' --header "Authorization: Bearer ',
                         bearer,
                         '" -P ',
                         temp_dir)

  if(quiet == FALSE) print(paste0("Downloading: ", file_name))
  tmp <- system(wget_command, intern = T, ignore.stdout = TRUE, ignore.stderr = TRUE)

  r <- file_to_raster(file.path(temp_dir, product_id, year, day, file_name),
                      variable,
                      quality_flag_rm)

  return(r)
}

define_variable <- function(variable, product_id){

  if(is.null(variable)){
    if(product_id == "VNP46A1") variable <- "DNB_At_Sensor_Radiance_500m"
    if(product_id == "VNP46A2") variable <- "Gap_Filled_DNB_BRDF-Corrected_NTL"
    if(product_id %in% c("VNP46A3", "VNP46A4")) variable <- "NearNadir_Composite_Snow_Free"
  }

  return(variable)
}

define_date_name <- function(date_i, product_id){

  #### Make name for raster based on date
  if(product_id %in% c("VNP46A1", "VNP46A2")){
    date_name_i <- paste0("t", date_i %>% str_replace_all("-", "_"))
  }

  if(product_id %in% c("VNP46A3")){
    date_name_i <- paste0("t", date_i %>% str_replace_all("-", "_") %>% substring(1,7))
  }

  if(product_id %in% c("VNP46A4")){
    date_name_i <- paste0("t", date_i %>% str_replace_all("-", "_") %>% substring(1,4))
  }

  return(date_name_i)
}

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
#' @param bearer NASA bearer token. For instructions on how to create a token, see [here](https://github.com/ramarty/blackmarbler#bearer-token-).
#' @param aggregation_fun Function used to aggregate nighttime lights data to polygons; this values is passed to the `fun` argument in [exactextractr::exact_extract](https://github.com/isciences/exactextractr) (Default: `mean`).
#' @param variable Variable to used to create raster (default: `NULL`). If `NULL`, uses the following default variables:
#' * For `product_id` `:VNP46A1"`, uses `DNB_At_Sensor_Radiance_500m`.
#' * For `product_id` `"VNP46A2"`, uses `Gap_Filled_DNB_BRDF-Corrected_NTL`.
#' * For `product_id`s `"VNP46A3"` and `"VNP46A4"`, uses `NearNadir_Composite_Snow_Free`.
#' For information on other variable choices, see [here](https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.2_April_2021.pdf); for `VNP46A1`, see Table 3; for `VNP46A2` see Table 6; for `VNP46A3` and `VNP46A4`, see Table 9.
#' @param quality_flag_rm Quality flag values to use to set values to `NA`. Each pixel has a quality flag value, where low quality values can be removed. Values are set to `NA` for each value in ther `quality_flag_rm` vector. (Default: `c(255)`).
#'
#'
#' For `VNP46A1` and `VNP46A2` (daily data):
#' - `0`: High-quality, Persistent nighttime lights
#' - `1`: High-quality, Ephemeral nighttime Lights
#' - `2`: Poor-quality, Outlier, potential cloud contamination, or other issues
#' - `255`: No retrieval, Fill value (masked out on ingestion)
#'
#'
#' For `VNP46A3` and `VNP46A4` (monthly and annual data):
#' - `0`: Good-quality, The number of observations used for the composite is larger than 3
#' - `1`: Poor-quality, The number of observations used for the composite is less than or equal to 3
#' - `2`: Gap filled NTL based on historical data
#' - `255`: Fill value
#' @param check_all_tiles_exist Check whether all Black Marble nighttime light tiles exist for the region of interest. Sometimes not all tiles are available, so the full region of interest may not be covered. If `TRUE`, skips cases where not all tiles are available. (Default: `TRUE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a raster in R. If `file`, function exports a `.tif` file and returns `NULL`.
#'
#' For `output_location_type = file`:
#' @param file_dir The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].tif`
#' @param file_skip_if_exists Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
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
                       variable = NULL,
                       quality_flag_rm = 255,
                       check_all_tiles_exist = TRUE,
                       output_location_type = "memory", # memory, file
                       file_dir = NULL,
                       file_prefix = NULL,
                       file_skip_if_exists = TRUE,
                       quiet = FALSE){

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
            r_agg <- exact_extract(x = r, y = roi_sf, fun = aggregation_fun)
            roi_df <- roi_sf
            roi_df$geometry <- NULL

            if(length(aggregation_fun) > 1){
              names(r_agg) <- paste0("ntl_", names(r_agg))
              r_agg <- bind_cols(r_agg, roi_df)
            } else{
              roi_df[[paste0("ntl_", aggregation_fun)]] <- r_agg
              r_agg <- roi_df
            }
            r_agg$date <- date_i

            #### Export
            saveRDS(r_agg, out_path)

          } else{
            cat(paste0('"', out_path, '" already exists; skipping.\n'))
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

          r_out <- exact_extract(x = r_out, y = roi_sf, fun = aggregation_fun)
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
#' @param bearer NASA bearer token. For instructions on how to create a token, see [here](https://github.com/ramarty/blackmarbler#bearer-token-).
#' @param variable Variable to used to create raster (default: `NULL`). If `NULL`, uses the following default variables:
#' * For `product_id` `:VNP46A1"`, uses `DNB_At_Sensor_Radiance_500m`.
#' * For `product_id` `"VNP46A2"`, uses `Gap_Filled_DNB_BRDF-Corrected_NTL`.
#' * For `product_id`s `"VNP46A3"` and `"VNP46A4"`, uses `NearNadir_Composite_Snow_Free`.
#' For information on other variable choices, see [here](https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.2_April_2021.pdf); for `VNP46A1`, see Table 3; for `VNP46A2` see Table 6; for `VNP46A3` and `VNP46A4`, see Table 9.
#' @param quality_flag_rm Quality flag values to use to set values to `NA`. Each pixel has a quality flag value, where low quality values can be removed. Values are set to `NA` for each value in ther `quality_flag_rm` vector. (Default: `c(255)`).
#'
#'
#' For `VNP46A1` and `VNP46A2` (daily data):
#' - `0`: High-quality, Persistent nighttime lights
#' - `1`: High-quality, Ephemeral nighttime Lights
#' - `2`: Poor-quality, Outlier, potential cloud contamination, or other issues
#' - `255`: No retrieval, Fill value (masked out on ingestion)
#'
#'
#' For `VNP46A3` and `VNP46A4` (monthly and annual data):
#' - `0`: Good-quality, The number of observations used for the composite is larger than 3
#' - `1`: Poor-quality, The number of observations used for the composite is less than or equal to 3
#' - `2`: Gap filled NTL based on historical data
#' - `255`: Fill value
#' @param check_all_tiles_exist Check whether all Black Marble nighttime light tiles exist for the region of interest. Sometimes not all tiles are available, so the full region of interest may not be covered. If `TRUE`, skips cases where not all tiles are available. (Default: `TRUE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a raster in R. If `file`, function exports a `.tif` file and returns `NULL`.
#'
#' For `output_location_type = file`:
#' @param file_dir The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].tif`
#' @param file_skip_if_exists Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
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
#' @import purrr
#' @import furrr
#' @import stringr
#' @import rhdf5
#' @import dplyr
#' @import sf
#' @import lubridate
#' @import readr
#' @import exactextractr
#' @import purrr
#' @rawNamespace import(raster, except = c(union, select, intersect, origin, tail, head))

# @rawNamespace import(utils, except = c(stack, unstack))
bm_raster <- function(roi_sf,
                      product_id,
                      date,
                      bearer,
                      variable = NULL,
                      quality_flag_rm = 255,
                      check_all_tiles_exist = TRUE,
                      output_location_type = "memory", # memory, file
                      file_dir = NULL,
                      file_prefix = NULL,
                      file_skip_if_exists = TRUE,
                      quiet = FALSE){

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
            cat(paste0('"', out_path, '" already exists; skipping.\n'))
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

  unlink(temp_dir, recursive = T)

  return(r)
}

bm_raster_i <- function(roi_sf,
                        product_id,
                        date,
                        bearer,
                        variable,
                        quality_flag_rm,
                        check_all_tiles_exist,
                        quiet,
                        temp_dir){

  # Checks ---------------------------------------------------------------------
  if(!("sf" %in% class(roi_sf))){
    stop("roi must be an sf object")
  }

  # Black marble grid ----------------------------------------------------------
  bm_tiles_sf <- read_sf("https://raw.githubusercontent.com/ramarty/blackmarbler/main/data/blackmarbletiles.geojson")

  # Prep dates -----------------------------------------------------------------
  ## For monthly, allow both yyyy-mm and yyyy-mm-dd (where -dd is ignored)
  if(product_id == "VNP46A3"){

    if(nchar(date) %in% 7){
      date <- paste0(date, "-01")
    }

  }

  ## For year, allow both yyyy and yyyy-mm-dd (where -mm-dd is ignored)
  if(product_id == "VNP46A4"){

    if(nchar(date) %in% 4){
      date <- paste0(date, "-01-01")
    }

  }

  # Grab tile dataframe --------------------------------------------------------
  #product_id <- "VNP46A4"
  #date <- "2021-10-15"

  year  <- date %>% year()
  month <- date %>% month()
  day   <- date %>% yday()

  bm_files_df <- create_dataset_name_df(product_id = product_id,
                                        all = T,
                                        years = year,
                                        months = month,
                                        days = day)

  # Intersecting tiles ---------------------------------------------------------
  # Remove grid along edges, which causes st_intersects to fail
  bm_tiles_sf <- bm_tiles_sf[!(bm_tiles_sf$TileID %>% str_detect("h00")),]
  bm_tiles_sf <- bm_tiles_sf[!(bm_tiles_sf$TileID %>% str_detect("v00")),]

  #inter <- st_intersects(bm_tiles_sf, roi_1row_sf, sparse = F) %>% as.vector()
  # inter <- st_intersects(bm_tiles_sf, roi_sf, sparse = F) %>%
  #   apply(1, sum)


  inter <- tryCatch(
    {
      inter <- st_intersects(bm_tiles_sf, roi_sf, sparse = F) %>%
        apply(1, sum)

      inter
    },
    error = function(e){
      warning("Issue with `roi_sf` intersecting with blackmarble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
      stop("Issue with `roi_sf` intersecting with blackmarble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
      #stop(st_intersects(bm_tiles_sf, roi_sf, sparse = F))
    }
  )


  grid_use_sf <- bm_tiles_sf[inter>0,]

  # Make Raster ----------------------------------------------------------------
  tile_ids_rx <- grid_use_sf$TileID %>% paste(collapse = "|")
  bm_files_df <- bm_files_df[bm_files_df$name %>% str_detect(tile_ids_rx),]

  if( (nrow(bm_files_df) < nrow(grid_use_sf)) & check_all_tiles_exist){
    warning("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
    stop("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
  }

  unlink(file.path(temp_dir, product_id), recursive = T)

  r_list <- lapply(bm_files_df$name, function(name_i){
    #if(quiet == FALSE) print(paste0("Downloading ", nrow(bm_files_df), " tiles."))
    download_raster(name_i, temp_dir, variable, bearer, quality_flag_rm, quiet)
  })

  if(length(r_list) == 1){
    r <- r_list[[1]]
  } else{

    ## Mosaic rasters together
    names(r_list)    <- NULL
    r_list$fun       <- max

    r <- do.call(raster::mosaic, r_list)

  }

  ## Crop
  r <- r %>% crop(roi_sf)

  unlink(file.path(temp_dir, product_id), recursive = T)

  return(r)
}


