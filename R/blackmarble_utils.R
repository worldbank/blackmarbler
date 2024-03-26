# box::use(
#   terra[rast, crop, mosaic, approximate, ext],
#   sf[read_sf, st_intersects, st_drop_geometry],
#   purrr[list_rbind, map2],
#   hdf5r[h5file],
#   stringr[stringr::str_replace_all, str_detect, str_extract],
#   readr[read_csv],
#   dplyr[mutate, across, summarise, bind_cols, bind_rows],
#   httr2[req_headers, request, req_perform, req_user_agent, resp_status, req_progress],
#   lubridate[year, month, yday],
#   exactextractr[exact_extract]
# )

#' Map Black Marble Tiles
#'
#' Function to generate an interactive map displaying black marble tiles using Leaflet.
#'
#' This function requires the \code{leaflet} package for generating interactive maps.
#'
#' @details
#' This function generates an interactive map displaying black marble tiles obtained from NASA's Black Marble project using the Leaflet library. The map includes black marble tiles overlaid on a base map, with a marker indicating the center coordinates (latitude and longitude).
#'
#' @seealso
#' \url{https://leafletjs.com/}
#'
#' @return An interactive Leaflet map object.
#'
#' @examples
#' library(leaflet)
#' # Generate the map
#' map_black_marble_tiles()
map_black_marble_tiles <- function() {
  # Define the center coordinates and zoom level
  center <- c(19.432608, -99.133209) # Mexico City coordinates

  zoom <- 2 # Adjust the zoom level as needed

  map <- leaflet::leaflet() |>
    leaflet::addProviderTiles("Stadia.AlidadeSmooth") |>
    leaflet::addPolygons(
      data = black_marble_tiles_sf,
      color = "#35B779",
      weight = 3,
      fillOpacity = 0,
      opacity = 1
    ) |>
    leaflet::addProviderTiles("NASAGIBS.ViirsEarthAtNight2012", options = leaflet::providerTileOptions(opacity = 0.5)) |>
    leaflet::addMarkers(lng = center[2], lat = center[1], popup = paste("Latitude:", center[1], "<br>Longitude:", center[2])) |> # Add marker with lat lon info
    leaflet::setView(lng = center[2], lat = center[1], zoom = zoom)

  return(map)
}


#' Translate Julian Dates to Regular Month Representation
#'
#' This function translates Julian dates to regular month representations.
#'
#' @param julian_date A character string representing the day of the year in Julian format (e.g., "001" for January 1st).
#' @return A number representing the month corresponding to the given Julian date (e.g., 1 for January).
#' @examples
#' julian_to_month("001")
#' # [1] 1
#'
#' julian_to_month("032")
#' # [1] 2
#' @references
#' For more information on the Julian day system, see: https://en.wikipedia.org/wiki/Julian_day
#' @export
julian_to_month <- function(julian_date) {
  months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  j_day <- as.integer(julian_date)
  month_index <- findInterval(j_day, c(0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366))
  return(as.numeric(months[month_index]))
}

#' Remove Artifact Values from Satellite Data
#'
#' This function removes artifact values from satellite data, replacing them with NA.
#' Artifact values are commonly used in satellite data to represent fill values,
#' invalid measurements, or missing data.
#'
#' The artifact values and their corresponding variables are based on the data format
#' described in the Black Marble User Guide (https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf).
#'  # * Table 3 (page 12)
# * Table 6 (page 16)
# * Table 9 (page 18)
#' @param data A matrix or data frame containing satellite data.
#' @param variable A character string specifying the variable for which artifact values should be removed.
#' @return The input data with artifact values replaced by NA.
#' @examples
#' # Remove artifact values from UTC_Time variable in satellite data
#' cleaned_data <- remove_fill_value_from_satellite_data(satellite_data, "UTC_Time")
#' @references
#' For more information on artifact values in satellite data, refer to the Black Marble User Guide.
#'
#' @export
remove_fill_value_from_satellite_data <- function(data, variable) {
  artifact_values_mapping <- list(
    list(255, c(
      "Granule", "Mandatory_Quality_Flag", "Latest_High_Quality_Retrieval",
      "Snow_Flag", "DNB_Platform", "Land_Water_Mask",
      "AllAngle_Composite_Snow_Covered_Quality", "AllAngle_Composite_Snow_Free_Quality",
      "NearNadir_Composite_Snow_Covered_Quality", "NearNadir_Composite_Snow_Free_Quality",
      "OffNadir_Composite_Snow_Covered_Quality", "OffNadir_Composite_Snow_Free_Quality"
    )),
    list(-999.9, "UTC_Time"),
    list(-32768, c(
      "Sensor_Azimuth", "Sensor_Zenith", "Solar_Azimuth", "Solar_Zenith",
      "Lunar_Azimuth", "Lunar_Zenith", "Glint_Angle", "Moon_Illumination_Fraction",
      "Moon_Phase_Angle"
    )),
    list(65535, c(
      "DNB_At_Sensor_Radiance_500m", "BrightnessTemperature_M12", "BrightnessTemperature_M13",
      "BrightnessTemperature_M15", "BrightnessTemperature_M16", "QF_Cloud_Mask", "QF_DNB",
      "QF_VIIRS_M10", "QF_VIIRS_M11", "QF_VIIRS_M12", "QF_VIIRS_M13", "QF_VIIRS_M15", "QF_VIIRS_M16",
      "Radiance_M10", "Radiance_M11", "DNB_BRDF-Corrected_NTL", "DNB_Lunar_Irradiance",
      "Gap_Filled_DNB_BRDF-Corrected_NTL", "AllAngle_Composite_Snow_Covered",
      "AllAngle_Composite_Snow_Covered_Num", "AllAngle_Composite_Snow_Free",
      "AllAngle_Composite_Snow_Free_Num", "NearNadir_Composite_Snow_Covered",
      "NearNadir_Composite_Snow_Covered_Num", "NearNadir_Composite_Snow_Free",
      "NearNadir_Composite_Snow_Free_Num", "OffNadir_Composite_Snow_Covered",
      "OffNadir_Composite_Snow_Covered_Num", "OffNadir_Composite_Snow_Free",
      "OffNadir_Composite_Snow_Free_Num", "AllAngle_Composite_Snow_Covered_Std",
      "AllAngle_Composite_Snow_Free_Std", "NearNadir_Composite_Snow_Covered_Std",
      "NearNadir_Composite_Snow_Free_Std", "OffNadir_Composite_Snow_Covered_Std",
      "OffNadir_Composite_Snow_Free_Std"
    ))
  )

  mapping_found <- FALSE


  for (mapping in artifact_values_mapping) {
    variables <- mapping[[2]]
    if (variable %in% variables || variable %in% unlist(variables)) {
      value <- mapping[[1]]
      data[data == value] <- NA
      mapping_found <- TRUE
      break # exit loop once artifact value is found
    }
  }

  if (!mapping_found) {
    warning(paste("Variable", variable, "not found in artifact values mapping. No action taken."))
  }

  return(data)
}

#' Apply Scaling Factor to VIIRS Data
#'
#' Apply scaling factor to variables according to Black Marble user guide.
#' The scaling factor is 0.1 for specific VIIRS variables.
#'
#' @param x A numeric vector or matrix representing the VIIRS data.
#' @param variable A character string specifying the variable name.
#' @param quiet Logical, if TRUE, suppresses warnings when the variable is not found in scaling variables.
#'
#' @return A numeric vector or matrix with the scaling factor applied to the specified variables.
#'
#' @details
#' This function applies a scaling factor of 0.1 to specific VIIRS variables
#' according to the Black Marble user guide.
#'
#' The following VIIRS variables are affected:
#' \itemize{
#'   \item DNB_At_Sensor_Radiance (VNP46A1)
#'   \item DNB_BRDF-Corrected_NTL (VNP46A2)
#'   \item Gap_Filled_DNB_BRDF-Corrected_NTL (VNP46A2)
#'   \item DNB_Lunar_Irradiance (VNP46A2)
#'   \item AllAngle_Composite_Snow_Covered (VNP46A3/4)
#'   \item AllAngle_Composite_Snow_Covered_Std (VNP46A3/4)
#'   \item AllAngle_Composite_Snow_Free (VNP46A3/4)
#'   \item AllAngle_Composite_Snow_Free_Std (VNP46A3/4)
#'   \item NearNadir_Composite_Snow_Covered (VNP46A3/4)
#'   \item NearNadir_Composite_Snow_Covered_Std (VNP46A3/4)
#'   \item NearNadir_Composite_Snow_Free (VNP46A3/4)
#'   \item NearNadir_Composite_Snow_Free_Std (VNP46A3/4)
#'   \item OffNadir_Composite_Snow_Covered (VNP46A3/4)
#'   \item OffNadir_Composite_Snow_Covered_Std (VNP46A3/4)
#'   \item OffNadir_Composite_Snow_Free (VNP46A3/4)
#'   \item OffNadir_Composite_Snow_Free_Std (VNP46A3/4)
#' }
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
apply_scaling_factor_to_viirs_data <- function(x, variable, quiet = TRUE) {
  # Apply scaling factor to VIIRS data
  scaling_variables <- c(
    "DNB_At_Sensor_Radiance",
    "DNB_BRDF-Corrected_NTL",
    "Gap_Filled_DNB_BRDF-Corrected_NTL",
    "DNB_Lunar_Irradiance",
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
    "OffNadir_Composite_Snow_Free_Std"
  )

  if (variable %in% scaling_variables) {
    x <- x * 0.1
  } else {
    if (!quiet) {
      warning(paste("Variable", variable, "not found in scaling variables. No scaling applied."))
    }
  }

  return(x)
}

#' Download VIIRS Satellite Image in HDF5 Format
#'
#' Downloads a VIIRS satellite image in HDF5 format from NASA's LADSWeb and saves it to a temporary directory or persistent location.
#'
#' @param file_name A character string representing the name of the file to download.
#' @param temp_dir A character string specifying the temporary directory where the file will be saved.
#' @param bearer A character string containing the authorization token for accessing NASA's LADSWeb.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data (optional).
#' @param quiet Logical; indicating whether to suppress progress messages (default: FALSE).
#'
#'
#' @details This function downloads a VIIRS satellite image in HDF5 format from NASA's LADSWeb
#' based on the provided file name. It then processes the downloaded data, removing quality flag values
#' if specified. The function also provides an option to suppress progress messages.
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
download_h5_viirs_sat_image <- function(file_name,
                                        download_path,
                                        bearer,
                                        quality_flags_to_remove = numeric(),
                                        quiet = FALSE) {
  # Extract file metadata
  year <- substr(file_name, 10, 13)
  day <- substr(file_name, 14, 16)
  product_id <- substr(file_name, 1, 7)


  # Construct download URL
  url <- paste0(
    "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/",
    product_id, "/", year, "/", day, "/", file_name
  )

  # Define url request
  # This is because regular httr2 setting doesn't return an integer in libcurl version, generating an API error related to versioning

  versions <- c(
    httr2 = utils::packageVersion("httr2"),
    `r-curl` = utils::packageVersion("curl"),
    libcurl = sub("-DEV", "", curl::curl_version()$version)
  )
  string <- paste0(names(versions), "/", versions, collapse = " ")

  request <- httr2::request(url) |>
    httr2::req_headers(
      "Authorization" = paste("Bearer", bearer)
    ) |>
    httr2::req_user_agent(string)

  # Display processing message if not quiet
  if (!quiet) message(paste0("Processing: ", file_name))

  # Perform the download
  if (quiet) {
    response <- request |>
      httr2::req_perform(
        path = download_path
      )
  } else {
    response <- request |>
      httr2::req_progress(type = "down") |>
      httr2::req_perform(
        path = download_path
      )
  }

  # Check for successful download
  if (httr2::resp_status(response) != 200) {
    message("Error in downloading data")
    message(response |>
      httr2::resp_status_desc())
  }
}
#' Extracts bounding box coordinates based on tile ID from the file path
#'
#' This code snippet extracts the bounding box coordinates (minimum longitude, minimum latitude, maximum longitude, maximum latitude)
#' based on the tile ID extracted from the file path. It uses the `file_path` to extract the tile ID, then retrieves the corresponding
#' grid from the `black_marble_tiles_sf` dataset. Finally, it calculates the bounding box of the grid and rounds the coordinates.
#'
#' @param file_path A character string specifying the file path.
#' @param black_marble_tiles_sf The spatial dataset containing the grid tiles information.
#'
#' @return A list containing the minimum and maximum longitude and latitude coordinates of the bounding box.
#'
#' @export
extract_bounding_box <- function(file_path, black_marble_tiles_sf) {
  tile_i <- file_path |>
    stringr::str_extract("h\\d{2}v\\d{2}")

  grid_i_sf <- black_marble_tiles_sf[black_marble_tiles_sf$TileID %in% tile_i, ]

  grid_i_sf_box <- grid_i_sf |>
    sf::st_bbox()

  min_lon <- min(grid_i_sf_box$xmin) |> round()
  min_lat <- min(grid_i_sf_box$ymin) |> round()
  max_lon <- max(grid_i_sf_box$xmax) |> round()
  max_lat <- max(grid_i_sf_box$ymax) |> round()

  return(list(min_lon = min_lon, min_lat = min_lat, max_lon = max_lon, max_lat = max_lat))
}

#' Extract Daily Data from HDF5 File
#'
#' This function extracts daily data from an HDF5 file based on the provided file path and variable name.
#'
#' @param file_path A character string specifying the file path.
#' @param h5_data The HDF5 file data.
#' @param variable_name A character string specifying the variable name.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data.
#'
#' @return A list containing the extracted data and metadata.
#'
#' @details This function extracts daily data from an HDF5 file for the specified variable.
#' It retrieves the data and applies quality flag removal if specified.
#' Additionally, it calculates the bounding box coordinates based on the tile ID extracted from the file path.
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
extract_daily_data <- function(file_path,
                               h5_data,
                               variable_name, quality_flags_to_remove) {
  allowed_variables_for_daily_data <- c(
    "DNB_At_Sensor_Radiance",
    "DNB_BRDF-Corrected_NTL",
    "Gap_Filled_DNB_BRDF-Corrected_NTL",
    "DNB_Lunar_Irradiance",
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
    "OffNadir_Composite_Snow_Free_Std"
  )

  if (!(variable_name %in% allowed_variables_for_daily_data)) {
    stop("Variable name must be one of the specified values.")
  }

  out <- h5_data[[paste0("HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/", variable_name)]][, ]
  qf <- h5_data[["HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/Mandatory_Quality_Flag"]][, ]

  # Check if there are quality flags to remove and if the variable is applicable
  if (length(quality_flags_to_remove) > 0 && variable_name %in% c(
    "DNB_BRDF-Corrected_NTL",
    "Gap_Filled_DNB_BRDF-Corrected_NTL",
    "Latest_High_Quality_Retrieval"
  )) {
    # Iterate over each quality flag value to remove
    for (flag_value in quality_flags_to_remove) {
      # Set values to NA where quality flag matches the specified value
      out[qf == flag_value] <- NA
    }
  }

  # Call the extract_bounding_box function to get the bounding box coordinates
  bounding_box <- extract_bounding_box(file_path, black_marble_tiles_sf)

  # Extract the coordinates from the bounding box list
  min_lon <- bounding_box$min_lon
  min_lat <- bounding_box$min_lat
  max_lon <- bounding_box$max_lon
  max_lat <- bounding_box$max_lat

  # Return the extracted data and bounding box coordinates
  return(list(data = out, min_lon = min_lon, max_lon = max_lon, min_lat = min_lat, max_lat = max_lat))
}

#' Extract Monthly Data from HDF5 File
extract_monthly_data <- function(h5_data, variable_name, quality_flags_to_remove) {
  # Extracting monthly data logic from the original function
  lat <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lat"]][]
  lon <- h5_data[["HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/lon"]][]

  out <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", variable_name)]][, ]

  # Check if there are quality flags to remove
  if (length(quality_flags_to_remove) > 0) {
    # Extract the base variable name without "_Num" or "_Std"
    variable_short <- variable_name |>
      stringr::str_remove_all("_Num") |>
      stringr::str_remove_all("_Std")

    # Construct the quality flag variable name
    qf_name <- paste0(variable_short, "_Quality")

    # Check if the quality flag variable exists in the data
    if (qf_name %in% names(h5_data)) {
      # Extract the quality flag data
      qf <- h5_data[[paste0("HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/", qf_name)]][, ]

      # Set values to NA where quality flag matches the specified value
      for (val in quality_flags_to_remove) {
        out[qf == val] <- NA
      }
    }
  }


  # Check if the data type of the first element is not numeric
  if (class(out[1, 1]) != "numeric") {
    # Convert the entire matrix to numeric
    out <- as.matrix(as.numeric(out))
  }


  # Extract lon/lat if available
  # if not i think the rest fails. so we need to check if it is available
  if ("lon" %in% colnames(out) && "lat" %in% colnames(out)) {
    min_lon <- min(out$lon)
    max_lon <- max(out$lon)
    min_lat <- min(out$lat)
    max_lat <- max(out$lat)
  } else {
    min_lon <- NULL
    max_lon <- NULL
    min_lat <- NULL
    max_lat <- NULL
  }

  return(list(data = out, min_lon = min_lon, max_lon = max_lon, min_lat = min_lat, max_lat = max_lat))
}


#' Extract Data and Metadata from HDF5 File
#'
#' This function extracts data and metadata from an HDF5 file.
#'
#' @param h5_data The HDF5 file data.
#' @param file_path A character string specifying the file path.
#' @param variable_name A character string specifying the variable name.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data.
#'
#' @return A list containing the extracted data and metadata.
#'
#' @details This function extracts data and metadata from an HDF5 file based on the specified file path and variable name.
#' Depending on the file path pattern (daily or monthly/annually), it calls different functions to extract the data.
#' Quality flag values specified in \code{quality_flags_to_remove} are removed from the data.
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
extract_data_and_metadata_from_hdf5 <- function(h5_data,
                                                download_path,
                                                variable_name,
                                                quality_flags_to_remove) {
  if (grepl("VNP46A1|VNP46A2", download_path, ignore.case = TRUE)) {
    print("im in daily_result")
    # Extract data for daily files
    daily_result <- extract_daily_data(
      download_path,
      h5_data,
      variable_name,
      quality_flags_to_remove
    )

    data <- daily_result$data
    min_lon <- daily_result$min_lon
    max_lon <- daily_result$max_lon
    min_lat <- daily_result$min_lat
    max_lat <- daily_result$max_lat
  } else {
    print("im in monthly_result")
    # Extract data for monthly/annually files
    monthly_result <- extract_monthly_data(h5_data, variable_name, quality_flags_to_remove)
    data <- monthly_result$data
    min_lon <- monthly_result$min_lon
    max_lon <- monthly_result$max_lon
    min_lat <- monthly_result$min_lat
    max_lat <- monthly_result$max_lat
  }

  # Construct metadata
  metadata <-
    list(
      nRows = nrow(data),
      nCols = ncol(data),
      res = nrow(data),
      nodata_val = NA,
      myCrs = "epsg:4326", # default crs for VIIRS data
      min_lon = min_lon,
      max_lon = max_lon,
      min_lat = min_lat,
      max_lat = max_lat
    )

  return(list(data = data, metadata = metadata))
}
#' Create Raster Object from Data and Metadata using terra
#'
#' This function creates a raster object from the provided data and metadata using the terra package.
#'
#' @param data The data matrix for the raster.
#' @param metadata A list containing metadata information including nodata value, bounding box coordinates, and coordinate reference system (CRS).
#'
#' @return A raster object created from the input data and metadata.
#'
#' @details This function transposes the input data, assigns nodata values to NA, creates an extent class from bounding box coordinates, and finally creates a raster object with the given data, extent, and CRS.
#'
#' @export
create_raster_from_data_metadata <- function(data, metadata) {
  #transpose data to fix flipped row and column order
  #depending upon how your data are formatted you might not have to perform this
 # transposed_data <- t(data)

  # Assign nodata values to NA
  # transposed_data[transposed_data == metadata$nodata_val] <- NA
  data[data == metadata$nodata_val] <- NA

  # Create extent class
  rasExt <- terra::ext(c(
    metadata$min_lon,
    metadata$max_lon,
    metadata$min_lat,
    metadata$max_lat
  ))

  # Create raster object
  my_raster <- terra::rast(data,
    extent = rasExt,
    crs = metadata$myCrs
  )

  return(my_raster)
}
#' Clean Raster Data
#'
#' This function cleans raster data by removing fill values and applying scaling factors.
#'
#' @param raster_obj The raster object to be cleaned.
#' @param variable_name A character string specifying the variable name.
#'
#' @return A cleaned raster object.
#'
#' @details This function removes fill values and applies scaling factors to the specified variable in the raster object.
#'
#' @seealso \code{\link{remove_fill_value_from_satellite_data}}, \code{\link{apply_scaling_factor_to_viirs_data}}
#'
#' @export
clean_raster_data <- function(raster_obj, variable_name) {
  # Remove fill values
  raster_obj <- remove_fill_value_from_satellite_data(raster_obj, variable_name)

  # Apply scaling factor
  raster_obj <- apply_scaling_factor_to_viirs_data(raster_obj, variable_name)

  return(raster_obj)
}

#' Convert HDF5 File to Raster
#'
#' Converts an HDF5 file to a raster object.
#'
#' @param file_path A character string representing the filepath to the HDF5 file.
#' @param variable_name A character string specifying the variable name to extract from the HDF5 file.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data (optional).
#'
#' @return A raster object containing the extracted variable data from the HDF5 file.
#'
#' @details
#' This function converts an HDF5 file to a raster object. It extracts the specified variable
#' from the HDF5 file and optionally removes specific quality flag values from the data.
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
convert_h5_to_raster <- function(download_path,
                                 variable_name,
                                 quality_flags_to_remove = numeric()) {
  # Load HDF5 file
  h5_data <- hdf5r::h5file(download_path, "r+")

  # Extract data and metadata
  # print("extract_data_and_metadata_from_hdf5")

  result_metadata_list <- extract_data_and_metadata_from_hdf5(
    h5_data,
    download_path,
    variable_name,
    quality_flags_to_remove
  )

  data <- result_metadata_list$data

  metadata <- result_metadata_list$metadata

  print("create_raster")
  # Convert data to raster

  raster_obj <- create_raster_from_data_metadata(data, metadata)

  print("clean_raster_data")

  # Clean raster data
  clean_raster_obj <- clean_raster_data(raster_obj, variable_name)

  # Close HDF5 file
  h5_data$close_all()

  return(clean_raster_obj)
}
#' Download and Convert Raster Data
#'
#' Downloads raster data from NASA's LADSWeb and converts it to a raster object.
#'
#' @param file_name A character string representing the name of the file to download.
#' @param temp_dir A character string specifying the temporary directory where the file will be saved.
#' @param variable A character string specifying the variable to extract from the raster data.
#' @param bearer A character string containing the authorization token for accessing NASA's LADSWeb.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data (optional).
#' @param quiet Logical; indicating whether to suppress progress messages (default: FALSE).
#'
#' @return A raster object containing the downloaded and processed raster data.
#'
#' @details This function downloads raster data from NASA's LADSWeb based on the provided file name.
#' It then converts the downloaded data to a raster object, extracting the specified variable and removing
#' quality flag values if specified. The function also provides an option to suppress progress messages.
#'
#' @references
#' Black Marble User Guide - https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf
#'
#' @export
download_and_convert_raster <- function(file_name,
                                        temp_dir,
                                        variable,
                                        bearer,
                                        quality_flags_to_remove = numeric(),
                                        quiet = FALSE) {

  # Define download path
  download_path <- file.path(temp_dir, file_name)

  # Download VIIRS satellite image in HDF5 format
  download_h5_viirs_sat_image(
    file_name,
    download_path,
    bearer,
    quality_flags_to_remove,
    quiet
  )

  # Convert downloaded file to raster
  raster_data <- convert_h5_to_raster(
    download_path,
    variable,
    quality_flags_to_remove
  )

  return(raster_data)
}

#' Read Black Marble CSV Data
#'
#' Reads Black Marble CSV data for a specific year and day.
#'
#' @param year The year of the data.
#' @param day The day of the year (1-365 or 1-366 for leap years).
#' @param product_id The product ID specifying the type of data to retrieve.
#'
#' @return A data frame containing the Black Marble CSV data for the specified year and day.
#'
#' @details This function reads Black Marble CSV data from the NASA LADS website for the given \code{year},
#' \code{day}, and \code{product_id}. It constructs the URL based on the provided parameters and attempts
#' to read the CSV file using \code{readr::read_csv}. If successful, it adds columns for \code{year} and \code{day}
#' to the data frame. If an error occurs during the reading process, it returns an empty data frame and issues a warning.
#' Additionally, to avoid overloading the server with rapid requests, the function includes a small delay (0.1 seconds) between
#' consecutive requests using \code{Sys.sleep(0.1)}.
#'
#' @examples
#' # Read Black Marble CSV data for the year 2023, day 150, and product ID "VIIRS_SNPP_CorrectedReflectance_BandM3"
#' data <- read_bm_csv(2023, 150, "VIIRS_SNPP_CorrectedReflectance_BandM3")
#'
#' @export
read_black_marble_csv <- function(year, day, product_id) {
  df_out <- tryCatch(
    {
      df <- readr::read_csv(paste0("https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/", product_id, "/", year, "/", day, ".csv"),
        show_col_types = FALSE
      )

      df$year <- year
      df$day <- day

      df
    },
    error = function(e) {
      warning(paste0("Error with year: ", year, "; day: ", day))
      data.frame(NULL)
    }
  )

  Sys.sleep(0.1) # Adding a small delay to avoid overloading the server

  return(df_out)
}

#' Create Black Marble Dataset DataFrame
#'
#' Creates a data frame containing Black Marble dataset filenames based on the specified parameters.
#'
#' @param product_id The product ID specifying the type of Black Marble data.
#' @param all Logical; indicating whether to create filenames for all available data or only for specific years, months, or days (default: TRUE).
#' @param years A numeric vector specifying the years for which to create filenames (optional).
#' @param months A numeric vector specifying the months for which to create filenames (optional).
#' @param days A numeric vector specifying the days for which to create filenames (optional).
#'
#' @return A data frame containing the filenames of Black Marble datasets based on the specified parameters.
#'
#' @details This function generates a data frame with filenames of Black Marble datasets. It allows filtering
#' by year, month, and day based on the provided parameters. Depending on the \code{product_id}, it generates
#' filenames for daily, monthly, or yearly data. The generated filenames are based on the specified product ID,
#' year, and day or month. The resulting data frame includes columns for the year and day or month, depending on the
#' type of data.
#'
#' @examples
#' # Generate filenames for all available Black Marble data
#' all_data <- create_black_marble_dataset_df("VNP46A1")
#'
#' # Generate filenames for Black Marble data for specific years and months
#' specific_data <- create_black_marble_dataset_df("VNP46A2", years = c(2018, 2019), months = 1:6)
#'
#' @export
create_black_marble_dataset_df <- function(product_id,
                                           all = TRUE,
                                           years = NULL,
                                           months = NULL,
                                           days = NULL) {
  # Define product-specific parameters and conditions
  product_params <- list(
    "VNP46A1" = list(months = NULL, days = 1:366, add_month = TRUE),
    "VNP46A2" = list(months = NULL, days = 1:366, add_month = TRUE),
    "VNP46A3" = list(months = NULL, days = c(
      "001", "032", "061", "092", "122", "153", "183", "214", "245", "275", "306", "336",
      "060", "091", "121", "152", "182", "213", "244", "274", "305", "335"
    ), add_month = TRUE),
    "VNP46A4" = list(months = NULL, days = "001", add_month = FALSE)
  )

  # Retrieve product-specific parameters
  params <- product_params[[product_id]]

  # Determine end year #what is this for?
  year_end <- as.numeric(format(Sys.Date(), "%Y"))

  # Generate parameter dataframe
  param_df <- tidyr::expand_grid(
    years = 2012:year_end,
    days = params$days # sprintf("%03d", params$days) days already come in the correct format
  )


  # Add month if required
  if (params$add_month) {
    param_df <- param_df |>
      dplyr::mutate(
        months = days |>
          purrr::map_int(julian_to_month)
        )
  }

  # Subset time period
  if (!is.null(years)) {
    param_df <- param_df[param_df$years %in% years, ]
  }

  if (!is.null(months)) {
    param_df <- param_df[as.numeric(param_df$months) %in% as.numeric(months), ]
  }

  if (!is.null(days)) {
    param_df <- param_df[as.numeric(param_df$days) %in% as.numeric(days), ]
  }

  # Create data
  files_df <- purrr::map2(
    param_df$years,
    param_df$days,
    read_black_marble_csv,
    product_id
  ) |>
    purrr::list_rbind()

  return(files_df)
}
#' Define Black Marble Variable
#'
#' Defines the variable based on the Black Marble product ID if it is NULL.
#'
#' @param variable A character string specifying the variable to define.
#' @param product_id A character string representing the product ID.
#'
#' @return A character string representing the defined variable.
#'
#' @export
define_blackmarble_variable <- function(variable, product_id) {
  if (is.null(variable)) {
    variable <- switch(product_id,
      "VNP46A1" = "DNB_At_Sensor_Radiance_500m",
      "VNP46A2" = "Gap_Filled_DNB_BRDF-Corrected_NTL",
      "VNP46A3" = "NearNadir_Composite_Snow_Free",
      "VNP46A4" = "NearNadir_Composite_Snow_Free",
      variable
    )
  }

  return(variable)
}



#' Define Raster Name by Date
#'
#' Generates a name for the raster based on the given date and product ID.
#'
#' @param date_string A character string representing the date in the format "YYYY-MM-DD".
#' @param product_id A character string representing the product ID (e.g., "VNP46A1", "VNP46A2").
#'
#' @return A character string representing the generated raster name.
#'
#' @export
define_raster_name <- function(date_string, product_id) {
  raster_name <- switch(product_id,
    "VNP46A1",
    "VNP46A2" = paste0("t", stringr::str_replace_all(date_string, "-", "_")),
    "VNP46A3" = paste0("t", stringr::str_replace_all(date_string, "-", "_") |> substring(1, 7)),
    "VNP46A4" = paste0("t", stringr::str_replace_all(date_string, "-", "_") |> substring(1, 4))
  )

  return(raster_name)
}

#' Count Observations for Exact Extract
#'
#' Counts observations for each variable, considering coverage fraction, for use in exact_extract.
#'
#' @param values A data frame containing the values.
#' @param coverage_fraction A numeric value representing the coverage fraction.
#'
#' @return A data frame with the count of non-NA pixels and total pixels for each variable.
#'
#' @export
count_n_obs <- function(values) {
  # coverage_fraction was a param but not used perhaps the other functions needs it
  orig_vars <- names(values)

  values |>
    dplyr::mutate(
      dplyr::across(orig_vars, ~ as.numeric(!is.na(.)))
    ) |>
    dplyr::summarise(
      dplyr::across(orig_vars, sum, .names = "n_non_na_pixels.{.col}"),
      dplyr::across(orig_vars, ~ length(.), .names = "n_pixels.{.col}")
    )
}

process_tiles <-
  function(bm_files_df,
           grid_use_sf,
           check_all_tiles_exist,
           temp_dir,
           product_id,
           variable,
           bearer,
           quality_flags_to_remove,
           quiet) {

    tile_ids_rx <- grid_use_sf$TileID |>
      paste(collapse = "|")

    selected_bm_files_df <- bm_files_df[bm_files_df$name |>
      stringr::str_detect(tile_ids_rx), ]


    if ((nrow(selected_bm_files_df) < nrow(grid_use_sf)) && check_all_tiles_exist) {
      message("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
      stop("Not all satellite imagery tiles for this location exist, so skipping. To ignore this error and process anyway, set check_all_tiles_exist = FALSE")
    }

    unlink(file.path(temp_dir, product_id), recursive = TRUE)

    if (!quiet) {
      message(paste0("Processing ", nrow(selected_bm_files_df), " nighttime light tiles"))
    }

    r_list <- lapply(selected_bm_files_df$name, function(name_i) {
      download_and_convert_raster(
        name_i,
        temp_dir,
        variable,
        bearer,
        quality_flags_to_remove,
        quiet
      )
    })

    if (length(r_list) == 1) {
      return(r_list[[1]])
    } else {
      names(r_list) <- NULL
      r_list$fun <- max
      return(do.call(terra::mosaic, r_list))
    }
  }
#' Intersect black marble tiles with region of interest
#'
#' This function intersects black marble tiles with a region of interest, removing grid along edges to avoid issues with \code{st_intersects}.
#'
#' @param black_marble_tiles_sf Spatial object representing black marble tiles.
#' @param roi_sf Spatial object representing the region of interest.
#' @return Spatial object containing black marble tiles intersecting with the region of interest.
#' @export
#' @import sf
#' @importFrom stringr str_detect
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(sf)
#'
#' # Create example data
#' black_marble_tiles_sf <- st_read("path_to_black_marble_tiles_shapefile")
#' roi_sf <- st_read("path_to_roi_shapefile")
#'
#' # Intersect black marble tiles with region of interest
#' intersected_tiles <- intersect_bm_tiles(black_marble_tiles_sf, roi_sf)
#' }
intersect_bm_tiles <- function(black_marble_tiles_sf, roi_sf) {
  # Remove grid along edges, which causes st_intersects to fail
  black_marble_tiles_sf <- black_marble_tiles_sf[!(stringr::str_detect(black_marble_tiles_sf$TileID, "h00") | stringr::str_detect(black_marble_tiles_sf$TileID, "v00")), ]

  inter <- tryCatch(
    {
      inter <- sf::st_intersects(black_marble_tiles_sf, roi_sf, sparse = FALSE) |>
        apply(1, sum)

      inter
    },
    error = function(e) {
      warning("Issue with `roi_sf` intersecting with black marble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
      stop("Issue with `roi_sf` intersecting with black marble tiles; try buffering by a width of 0: eg, st_buffer(roi_sf, 0)")
    }
  )

  grid_use_sf <- black_marble_tiles_sf[inter > 0, ]
  return(grid_use_sf)
}


#' Retrieve and Process Nightlight Data
#'
#' Retrieves and processes nighttime light raster data from the Black Marble dataset.
#'
#' @param roi_sf An sf object representing the region of interest.
#' @param product_id A character string specifying the product ID.
#' @param date A character string representing the date.
#' @param bearer A character string representing the authorization bearer token.
#' @param variable A character string specifying the variable to extract.
#' @param quality_flags_to_remove A numeric vector containing quality flag values to be removed from the data (optional).
#' @param check_all_tiles_exist A logical value indicating whether to check if all satellite imagery tiles for the location exist (default is TRUE).
#' @param quiet A logical value indicating whether to suppress progress messages (default is FALSE).
#' @param temp_dir A character string representing the temporary directory path.
#'
#' @return A raster object containing the processed nighttime light data.
#'
#' @details
#' This function retrieves and processes nighttime light raster data from the Black Marble dataset.
#' It downloads the data for the specified region of interest and date, removes any quality-flagged
#' pixels as specified, and mosaics the tiles together if necessary. It then crops the raster to the
#' region of interest.
#'
#' @export
retrieve_and_process_nightlight_data <- function(roi_sf,
                                                 product_id,
                                                 date,
                                                 bearer,
                                                 variable,
                                                 quality_flags_to_remove,
                                                 check_all_tiles_exist = TRUE,
                                                 quiet = FALSE,
                                                 temp_dir) {
  # Checks ---------------------------------------------------------------------
  if (!("sf" %in% class(roi_sf))) {
    stop("roi must be an sf object")
  }


  # Prep dates -----------------------------------------------------------------
  date <- switch(product_id,
    "VNP46A3" = ifelse(nchar(date) == 7, paste0(date, "-01"), date),
    "VNP46A4" = ifelse(nchar(date) == 4, paste0(date, "-01-01"), date),
    date
  )

  # Grab tile dataframe --------------------------------------------------------
  year <- date |> lubridate::year()
  month <- date |> lubridate::month()
  day <- date |> lubridate::yday()

  bm_files_df <- create_black_marble_dataset_df(
    product_id = product_id,
    all = T,
    years = year,
    months = month,
    days = day
  )

  # Intersecting tiles ---------------------------------------------------------
  # Remove grid along edges, which causes st_intersects to fail
  # Intersect black marble tiles with region of interest
  intersected_tiles <- intersect_bm_tiles(
    black_marble_tiles_sf,
    roi_sf
  )
  # Make Raster ----------------------------------------------------------------
  print("processing tiles...")

  raster <- process_tiles(
    bm_files_df,
    intersected_tiles,
    check_all_tiles_exist,
    temp_dir,
    product_id,
    variable,
    bearer,
    quality_flags_to_remove,
    quiet
  )

  ## Crop
  raster <- raster |>
    terra::crop(roi_sf)

  unlink(file.path(temp_dir, product_id), recursive = T)

  return(raster)
}

# Helper functions ------------------------------------------------------------

#' Extract and process raster data
extract_and_process <- function(raster, roi_sf, fun, add_n_pixels = TRUE, quiet) {
  extracted_data <- exactextractr::exact_extract(raster, roi_sf, fun, progress = !quiet)
  roi_df <- sf::st_drop_geometry(roi_sf)
  roi_df$date <- NULL

  if (add_n_pixels) {
    # Compute additional pixel information if add_n_pixels is TRUE
    roi_df$n_pixels <- exactextractr::exact_extract(raster, roi_sf, function(values, coverage_fraction) {
      sum(!is.na(values))
    },
    progress = !quiet
    )
    roi_df$n_non_na_pixels <- exactextractr::exact_extract(raster, roi_sf, function(values, coverage_fraction) {
      length(values)
    },
    progress = !quiet
    )
    roi_df$prop_non_na_pixels <- roi_df$n_non_na_pixels / roi_df$n_pixels
  }

  if (length(fun) > 1) {
    names(extracted_data) <- paste0("ntl_", names(extracted_data))
    extracted_data <- dplyr::bind_cols(extracted_data, roi_df)
  } else {
    roi_df[[paste0("ntl_", fun)]] <- extracted_data
    extracted_data <- roi_df
  }

  return(extracted_data)
}

#' Extract and process raster data for individual dates
extract_and_process_i <- function(roi_sf, product_id, date_i, bearer, variable,
                                  quality_flags_to_remove, check_all_tiles_exist, add_n_pixels = TRUE, quiet, temp_dir) {

  bm_r <- retrieve_and_process_nightlight_data(
    roi_sf = roi_sf,
    product_id = product_id,
    date = date_i,
    bearer = bearer,
    variable = variable,
    quality_flags_to_remove = quality_flags_to_remove,
    check_all_tiles_exist = check_all_tiles_exist,
    quiet = quiet,
    temp_dir = temp_dir
  )

  r_agg <- exactextractr::exact_extract(
    x = bm_r, y = roi_sf, fun = aggregation_fun,
    progress = !quiet
  )

  if (add_n_pixels) {
    # Compute additional pixel information if add_n_pixels is TRUE
    roi_sf$n_pixels <- exactextractr::exact_extract(bm_r, roi_sf, function(values, coverage_fraction) {
      sum(!is.na(values))
    },
    progress = !quiet
    )
    roi_sf$n_non_na_pixels <- exactextractr::exact_extract(bm_r, roi_sf, function(values, coverage_fraction) {
      length(values)
    },
    progress = !quiet
    )
    roi_sf$prop_non_na_pixels <- roi_sf$n_non_na_pixels / roi_sf$n_pixels
  }

  return(r_agg)
}

#' Bind extracted data
bind_extracted_data <- function(n_obs_df, ntl_df) {
  names(ntl_df)[names(ntl_df) != "date"] <- paste0("ntl_", names(ntl_df)[names(ntl_df) != "date"])
  ntl_df$date <- NULL

  r <- dplyr::bind_cols(n_obs_df, ntl_df)
  return(r)
}

#' Bind extracted data list
bind_extracted_data_list <- function(r_list) {
  r_list <- r_list[!sapply(r_list, is.null)]
  r <- dplyr::bind_rows(r_list)
  return(r)
}
