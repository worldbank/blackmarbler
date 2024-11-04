# BlackMarblerR

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

pad3 <- function(x){
  if(nchar(x) == 1) out <- paste0("00", x)
  if(nchar(x) == 2) out <- paste0("0", x)
  if(nchar(x) == 3) out <- paste0(x)
  return(out)
}
pad3 <- Vectorize(pad3)

remove_fill_value <- function(x, variable){
  # Remove fill values
  
  # https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
  # * Table 3 (page 12)
  # * Table 6 (page 16)
  # * Table 9 (page 18)
  
  #### 255
  if(variable %in% c(
    "Granule",
    "Mandatory_Quality_Flag",
    "Latest_High_Quality_Retrieval",
    "Snow_Flag",
    "DNB_Platform",
    "Land_Water_Mask",
    "AllAngle_Composite_Snow_Covered_Quality",
    "AllAngle_Composite_Snow_Free_Quality",
    "NearNadir_Composite_Snow_Covered_Quality",
    "NearNadir_Composite_Snow_Free_Quality",
    "OffNadir_Composite_Snow_Covered_Quality",
    "OffNadir_Composite_Snow_Free_Quality"
  )){
    x[][x[] == 255] <- NA
  }
  
  #### -999.9 
  if(variable %in% c("UTC_Time")){
    x[][x[] == -999.9] <- NA
  }
  
  #### -32768 
  if(variable %in% c("Sensor_Azimuth",
                     "Sensor_Zenith",
                     "Solar_Azimuth",
                     "Solar_Zenith",
                     "Lunar_Azimuth",
                     "Lunar_Zenith",
                     "Glint_Angle",
                     "Moon_Illumination_Fraction",
                     "Moon_Phase_Angle")){
    x[][x[] == -32768] <- NA
  }
  
  
  #### 65535
  if(variable %in% c(
    "DNB_At_Sensor_Radiance_500m",
    "BrightnessTemperature_M12",
    "BrightnessTemperature_M13",
    "BrightnessTemperature_M15",
    "BrightnessTemperature_M16",
    "QF_Cloud_Mask",
    "QF_DNB",
    "QF_VIIRS_M10",
    "QF_VIIRS_M11",
    "QF_VIIRS_M12",
    "QF_VIIRS_M13",
    "QF_VIIRS_M15",
    "QF_VIIRS_M16",
    "Radiance_M10",
    "Radiance_M11",
    "QF_Cloud_Mask",
    "DNB_BRDF-Corrected_NTL",
    "DNB_Lunar_Irradiance",
    "Gap_Filled_DNB_BRDF-Corrected_NTL",
    "AllAngle_Composite_Snow_Covered",
    "AllAngle_Composite_Snow_Covered_Num",
    "AllAngle_Composite_Snow_Free",
    "AllAngle_Composite_Snow_Free_Num",
    "NearNadir_Composite_Snow_Covered",
    "NearNadir_Composite_Snow_Covered_Num",
    "NearNadir_Composite_Snow_Free",
    "NearNadir_Composite_Snow_Free_Num",
    "OffNadir_Composite_Snow_Covered",
    "OffNadir_Composite_Snow_Covered_Num",
    "OffNadir_Composite_Snow_Free",
    "OffNadir_Composite_Snow_Free_Num",
    "AllAngle_Composite_Snow_Covered_Std",
    "AllAngle_Composite_Snow_Free_Std",
    "NearNadir_Composite_Snow_Covered_Std",
    "NearNadir_Composite_Snow_Free_Std",
    "OffNadir_Composite_Snow_Covered_Std",
    "OffNadir_Composite_Snow_Free_Std"
  )){
    x[][x[] == 65535] <- NA
  }
  
  return(x)
}

apply_scaling_factor <- function(x, variable){
  # Apply scaling factor to variables according to Black Marble user guide
  
  # https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
  # * Table 3 (page 12)
  # * Table 6 (page 16)
  # * Table 9 (page 18)
  
  if(variable %in% c(
    
    # VNP46A1
    "DNB_At_Sensor_Radiance",
    
    # VNP46A2
    "DNB_BRDF-Corrected_NTL",
    "Gap_Filled_DNB_BRDF-Corrected_NTL",
    "DNB_Lunar_Irradiance",
    
    # VNP46A3/4
    "AllAngle_Composite_Snow_Covered",
    "AllAngle_Composite_Snow_Covered_Std",
    "AllAngle_Composite_Snow_Free",
    "AllAngle_Composite_Snow_Free_Std",
    "NearNadir_Composite_Snow_Covered",
    "NearNadir_Composite_Snow_Covered_Std",
    "NearNadir_Composite_Snow_Free",
    "NearNadir_Composite_Snow_Free_Std",
    "OffNadir_Composite_Snow_Covered",
    "OffNadir_Composite_Snow_Covered_Std",
    "OffNadir_Composite_Snow_Free",
    "OffNadir_Composite_Snow_Free_Std")
  ){
    
    x <- x * 0.1
    
  }
  
  return(x)
}

file_to_raster <- function(h5_file,
                           variable,
                           quality_flag_rm){
  # Converts h5 file to raster.
  # ARGS
  # --f: Filepath to h5 file
  
  ## Data
  h5_data <- h5file(h5_file, "r+")
  
  #### Daily
  if(h5_file %>% str_detect("VNP46A1|VNP46A2")){
    
    tile_i <- h5_file %>% stringr::str_extract("h\\d{2}v\\d{2}")
    
    bm_tiles_sf <- read_sf("https://raw.githubusercontent.com/worldbank/blackmarbler/main/data/blackmarbletiles.geojson")
    grid_i_sf <- bm_tiles_sf[bm_tiles_sf$TileID %in% tile_i,]
    
    grid_i_sf_box <- grid_i_sf %>%
      st_bbox()
    
    xMin <- min(grid_i_sf_box$xmin) %>% round()
    yMin <- min(grid_i_sf_box$ymin) %>% round()
    xMax <- max(grid_i_sf_box$xmax) %>% round()
    yMax <- max(grid_i_sf_box$ymax) %>% round()
    
    var_names <- h5_data[["HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields"]]$names
    
    if(!(variable %in% var_names)){
      warning(paste0("'", variable, "'",
                     " not a valid variable option. Valid options include:\n",
                     paste(var_names, collapse = "\n")
      ))
    }
    
    out <- h5_data[[paste0("HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/", variable)]][,]
    qf  <- h5_data[["HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/Mandatory_Quality_Flag"]][,]
    
    if(length(quality_flag_rm) > 0){
      if(variable %in% c("DNB_BRDF-Corrected_NTL",
                         "Gap_Filled_DNB_BRDF-Corrected_NTL",
                         "Latest_High_Quality_Retrieval")){
        
        for(val in quality_flag_rm){ # out[qf %in% quality_flag_rm] doesn't work, so loop
          out[qf == val] <- NA
        }
      }
    }
    
    # # Above doesn't fully capture
    # if(variable %in% "Latest_High_Quality_Retrieval"){
    #   out[out == 255] <- NA
    # }
    
    #### Monthly/Annually
  } else{
    
    lat <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lat"]][]
    lon <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lon"]][]
    
    var_names <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"]]$names
    
    if(!(variable %in% var_names)){
      warning(paste0("'", variable, "'",
                     " not a valid variable option. Valid options include:\n",
                     paste(var_names, collapse = "\n")
      ))
    }
    
    out <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", variable)]][,]
    
    if(length(quality_flag_rm) > 0){
      
      variable_short <- variable %>%
        str_replace_all("_Num", "") %>%
        str_replace_all("_Std", "")
      
      qf_name <- paste0(variable_short, "_Quality")
      
      if(qf_name %in% var_names){
        
        qf <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", qf_name)]][,]
        
        for(val in quality_flag_rm){ # out[qf %in% quality_flag_rm] doesn't work, so loop
          out[qf == val] <- NA
        }
        
      }
      
    }
    
    if(class(out[1,1])[1] != "numeric"){
      out <- matrix(as.numeric(out),  # Convert to numeric matrix
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
  #nodata_val <- NA
  myCrs      <- "EPSG:4326"
  
  ## Make Raster
  
  #transpose data to fix flipped row and column order
  #depending upon how your data are formatted you might not have to perform this
  out <- t(out)
  
  #assign data ignore values to NA
  #out[out == nodata_val] <- NA
  
  #turn the out object into a raster
  outr <- terra::rast(out,
                      crs = myCrs,
                      extent = c(xMin,xMax,yMin,yMax))
  
  #set fill values to NA
  outr <- remove_fill_value(outr, variable)
  
  #apply scaling factor
  outr <- apply_scaling_factor(outr, variable)
  
  #h5closeAll()
  h5_data$close_all()
  
  return(outr)
}

read_bm_csv <- function(year,
                        day,
                        product_id){
  
  
  # 
  df_out <- tryCatch(
    {
      df <- readr::read_csv(paste0("https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/",product_id,"/",year,"/",day,".csv"),
                            show_col_types = F)
      
      
      df$year <- year
      df$day <- day
      
      df
    },
    error = function(e){
      #warning(paste0("Error with year: ", year, "; day: ", day))
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
    param_df <- cross_df(list(year = 2012:year_end,
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
  # files_df <- purrr::map2_dfr(param_df$year,
  #                             param_df$day,
  #                             read_bm_csv,
  #                             product_id)
  
  files_df <- purrr::map2(param_df$year,
                          param_df$day,
                          read_bm_csv,
                          product_id) %>%
    bind_rows()
  
  
  return(files_df)
}

download_raster <- function(file_name,
                            temp_dir,
                            variable,
                            bearer,
                            quality_flag_rm,
                            h5_dir,
                            quiet){
  
  year       <- file_name %>% substring(10,13)
  day        <- file_name %>% substring(14,16)
  product_id <- file_name %>% substring(1,7)
  
  url <- paste0('https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/',
                product_id, '/', year, '/', day, '/', file_name)
  
  headers <- c('Authorization' = paste('Bearer', bearer))
  
  if(is.null(h5_dir)){
    download_path <- file.path(temp_dir, file_name)
  } else{
    download_path <- file.path(h5_dir, file_name)
  }
  
  if(!file.exists(download_path)){
    
    if(quiet == FALSE) message(paste0("Processing: ", file_name))
    
    if(quiet == TRUE){
      
      response <- httr::GET(url, 
                            httr::timeout(60),
                            httr::add_headers(headers), 
                            httr::write_disk(download_path, overwrite = TRUE))
      
    } else{
      response <- httr::GET(url, 
                            httr::timeout(60),
                            httr::add_headers(headers), 
                            httr::write_disk(download_path, overwrite = TRUE),
                            httr::progress())
      
    }
    
    
    
    if(response$status_code != 200){
      message("Error in downloading data")
      message(response)
    }
    
    if(response$all_headers[[1]]$status != 200){
      message("**Error in downloading data; bearer token likely invalid.** Try regenerating the bearer token; please see this link for instructions to obtain a bearer token: https://github.com/worldbank/blackmarbler?tab=readme-ov-file#bearer-token-")
    }
    
  }
  
  r <- file_to_raster(download_path,
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


count_n_obs <- function(values, coverage_fraction) {
  ## Function to count observations, for exact_extract
  
  orig_vars <- names(values)
  
  values %>%
    dplyr::mutate(across(orig_vars, ~ as.numeric(!is.na(.)) )) %>%
    dplyr::summarise(across(orig_vars, sum, .names = "n_non_na_pixels.{.col}"),
                     across(orig_vars, ~length(.), .names = "n_pixels.{.col}"))
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
#' @param date Date of raster data. Entering one date will produce a `SpatRaster` object. Entering multiple dates will produce a `SpatRaster` object with multiple bands; one band per date.
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
#' @param interpol_na When data for more than one date is downloaded, whether to interpolate `NA` values in rasters using the `terra::approximate` function. Additional arguments for the `terra::approximate` function can also be passed into `bm_extract` (eg, `method`, `rule`, `f`, `ties`, `z`, `NA_rule`). (Default: `FALSE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a dataframe in R. If `file`, function exports a `.csv` file and returns `NULL`.
#' @param file_dir (If `output_location_type = file`). The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix (If `output_location_type = file`). Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].csv`
#' @param file_skip_if_exists (If `output_location_type = file`). Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param file_return_null Whether to return `NULL` instead of a `dataframe`. When `output_location_type = 'file'`, the function will export data to the `file_dir` directory. When `file_return_null = FALSE`, the function will also return a `dataframe` of the queried data---so the data is available in R memory. Setting `file_return_null = TRUE`, data will be saved to `file_dir` but no data will be returned by the function to R memory (default: `FALSE`).
#' @param h5_dir Black Marble data are originally downloaded as `h5` files. If `h5_dir = NULL`, the function downloads to a temporary directory then deletes the directory. If `h5_dir` is set to a path, `h5` files are saved to that directory and not deleted. The function will then check if the needed `h5` file already exists in the directory; if it exists, the function will not re-download the `h5` file.
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
#'
#' @param ... Additional arguments for `terra::approximate`, if `interpol_na = TRUE`
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
                       file_return_null = FALSE,
                       h5_dir = NULL,
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
  
  if(class(roi_sf)[1] == "SpatVector") roi_sf <- roi_sf %>% st_as_sf()
  if(!("sf" %in% class(roi_sf))){
    stop("roi must be an sf object")
  }
  
  # Required parameters used in try statement, so error not generated when used, 
  # so use them here
  roi_sf <- roi_sf
  product_id <- product_id
  date <- date
  bearer <- bearer
  
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
  
  # Filename root --------------------------------------------------------------
  # Define outside of lapply, as use this later to aggregate rasters
  if(output_location_type == "file"){
    out_name_begin <- paste0(file_prefix, 
                             product_id, "_", 
                             variable, "_",
                             "qflag", 
                             quality_flag_rm %>% paste0(collapse="_"), "_",
                             aggregation_fun %>% paste0(collapse="_"))
  }
  
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
                      h5_dir = h5_dir,
                      quiet = quiet,
                      temp_dir = temp_dir)
    
    bm_r <- terra::approximate(bm_r,
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
    
    # Apply through each date, extract, then append
  } else{
    
    # Download data --------------------------------------------------------------
    r_list <- lapply(date, function(date_i){
      
      #out <- tryCatch(
      #  {
      
      #### Make name for raster based on date
      date_name_i <- define_date_name(date_i, product_id)
      
      #### If save to file
      if(output_location_type == "file"){
        
        out_name_end <- paste0("_", date_name_i, ".Rds")
        out_name <- paste0(out_name_begin, out_name_end)
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
                           h5_dir = h5_dir,
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
            
            r_n_obs <- exact_extract(r, roi_sf, function(values, coverage_fraction)
              sum(!is.na(values)),
              progress = !quiet)
            
            r_n_obs_poss <- exact_extract(r, roi_sf, function(values, coverage_fraction)
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
                             h5_dir = h5_dir,
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
      
      ## HERE
      #  },
      #  error=function(e) {
      #    return(NULL)
      #  }
      #)
      
    })
    
    # Clean output ---------------------------------------------------------------
    # Remove NULLs
    r_list <- r_list[!sapply(r_list,is.null)]
    
    r <- r_list %>%
      bind_rows()
    
  }
  
  # Output dataframe when output_location_type = "file" ------------------------
  if(output_location_type == "file"){
    if(!file_return_null){
      
      print(out_name_begin)
      
      ## Output path
      date_names <- define_date_name(date, product_id)
      
      out_name_end <- paste0("_",
                             date_names,
                             ".Rds")
      out_name <- paste0(out_name_begin, out_name_end)
      
      r <- file.path(file_dir, out_name) %>%
        map_df(readRDS)
      
      # r <- file_dir %>%
      #   list.files(full.names = T,
      #              pattern = paste0("*.Rds")) %>%
      #   str_subset(out_name_begin) %>%
      #   map_df(readRDS)
    } else{
      r <- NULL
    }
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
#' @param date Date of raster data. Entering one date will produce a `SpatRaster` object. Entering multiple dates will produce a `SpatRaster` object with multiple bands; one band per date.
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
#' @param interpol_na When data for more than one date is downloaded, whether to interpolate `NA` values using the `terra::approximate` function. Additional arguments for the `terra::approximate` function can also be passed into `bm_raster` (eg, `method`, `rule`, `f`, `ties`, `z`, `NA_rule`). (Default: `FALSE`).
#' @param output_location_type Where to produce output; either `memory` or `file`. If `memory`, functions returns a raster in R. If `file`, function exports a `.tif` file and returns `NULL`.
#' For `output_location_type = file`:
#' @param file_dir The directory where data should be exported (default: `NULL`, so the working directory will be used)
#' @param file_prefix Prefix to add to the file to be saved. The file will be saved as the following: `[file_prefix][product_id]_t[date].tif`
#' @param file_skip_if_exists Whether the function should first check wither the file already exists, and to skip downloading or extracting data if the data for that date if the file already exists (default: `TRUE`).
#' @param file_return_null Whether to return `NULL` instead of a `SpatRaster`. When `output_location_type = 'file'`, the function will export data to the `file_dir` directory. When `file_return_null = FALSE`, the function will also return a `SpatRaster` of the queried data---so the data is available in R memory. Setting `file_return_null = TRUE`, data will be saved to `file_dir` but no data will be returned by the function to R memory (default: `FALSE`).
#' @param h5_dir Black Marble data are originally downloaded as `h5` files. If `h5_dir = NULL`, the function downloads to a temporary directory then deletes the directory. If `h5_dir` is set to a path, `h5` files are saved to that directory and not deleted. The function will then check if the needed `h5` file already exists in the directory; if it exists, the function will not re-download the `h5` file.
#' @param quiet Suppress output that show downloading progress and other messages. (Default: `FALSE`).
#' @param ... Additional arguments for `terra::approximate`, if `interpol_na = TRUE`
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
#' @rawNamespace import(terra, except = c(intersect, values, origin, union))
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
                      file_return_null = FALSE,
                      h5_dir = NULL,
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
  
  if(class(roi_sf)[1] == "SpatVector") roi_sf <- roi_sf %>% st_as_sf()
  if(!("sf" %in% class(roi_sf))){
    stop("roi must be an sf object")
  }
  
  # Required parameters used in try statement, so error not generated when used, 
  # so use them here
  roi_sf     <- roi_sf
  product_id <- product_id
  date       <- date
  bearer     <- bearer
  
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
  
  # Filename root --------------------------------------------------------------
  # Define outside of lapply, as use this later to aggregate rasters
  if(output_location_type == "file"){
    out_name_begin <- paste0(file_prefix, 
                             product_id, "_", 
                             variable, "_",
                             "qflag", 
                             quality_flag_rm %>% paste0(collapse="_"))
  }
  
  # Download data --------------------------------------------------------------
  r_list <- lapply(date, function(date_i){
    
    out <- tryCatch(
      {
        
        #### Make name for raster based on date
        date_name_i <- define_date_name(date_i, product_id)
        
        #### If save as tif format
        if(output_location_type == "file"){
          
          ## Output path
          out_name_end <- paste0("_",
                                 date_name_i,
                                 ".tif")
          out_name <- paste0(out_name_begin, out_name_end)
          
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
                             h5_dir = h5_dir,
                             quiet = quiet,
                             temp_dir = temp_dir)
            names(r) <- date_name_i
            
            writeRaster(r, out_path)
            
          } else{
            message(paste0('"', out_path, '" already exists; skipping.\n'))
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
                               h5_dir = h5_dir,
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
    r <- terra::rast(r_list)
  } else{
    r <- NULL
  }
  
  # Interpolate ----------------------------------------------------------------
  if(interpol_na %in% T){
    r <- terra::approximate(r,
                            method = method,
                            rule   = rule,
                            f      = f,
                            ties   = ties,
                            z      = z,
                            NArule = NArule)
  }
  
  unlink(temp_dir, recursive = T)
  
  # Output raster when output_location_type = "file" ---------------------------
  if(output_location_type == "file"){
    if(!file_return_null){
      
      ## Output path
      date_names <- define_date_name(date, product_id)
      
      out_name_end <- paste0("_",
                             date_names,
                             ".tif")
      out_name <- paste0(out_name_begin, out_name_end)
      
      r <- file.path(file_dir, out_name) %>%
        rast()

    } else{
      r <- NULL
    }
  }
  
  return(r)
}

bm_raster_i <- function(roi_sf,
                        product_id,
                        date,
                        bearer,
                        variable,
                        quality_flag_rm,
                        check_all_tiles_exist,
                        h5_dir,
                        quiet,
                        temp_dir){
  
  # Black marble grid ----------------------------------------------------------
  bm_tiles_sf <- read_sf("https://raw.githubusercontent.com/worldbank/blackmarbler/main/data/blackmarbletiles.geojson")
  
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
  
  
  inter <- tryCatch(
    {
      inter <- st_intersects(bm_tiles_sf, roi_sf, sparse = F) %>%
        apply(1, sum)
      
      inter
    },
    error = function(e){
      warning("Issue with `roi_sf` intersecting with blackmarble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
      stop("Issue with `roi_sf` intersecting with blackmarble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
    }
  )
  
  grid_use_sf <- bm_tiles_sf[inter > 0,]
  
  # Make Raster ----------------------------------------------------------------
  tile_ids_rx <- grid_use_sf$TileID %>% paste(collapse = "|")
  bm_files_df <- bm_files_df[bm_files_df$name %>% str_detect(tile_ids_rx),]
  
  if( (nrow(bm_files_df) < nrow(grid_use_sf)) & check_all_tiles_exist){
    warning("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
    stop("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
  }
  
  unlink(file.path(temp_dir, product_id), recursive = T)
  
  if(quiet == F){
    message(paste0("Processing ", nrow(bm_files_df), " nighttime light tiles"))
  }
  
  
  r_list <- lapply(bm_files_df$name, function(name_i){
    download_raster(name_i, temp_dir, variable, bearer, quality_flag_rm, h5_dir, quiet)
  })
  
  
  
  if(length(r_list) == 1){
    r <- r_list[[1]]
  } else{
    
    r <- do.call(terra::mosaic, c(r_list, fun = "max"))
  }
  
  ## Crop
  r <- r %>% terra::crop(roi_sf)
  
  unlink(file.path(temp_dir, product_id), recursive = T)
  
  return(r)
}


