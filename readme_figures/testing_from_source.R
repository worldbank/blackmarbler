# Testing from Source

#### Packages
library(readr)
library(hdf5r)
library(dplyr)
library(purrr)
library(lubridate)
library(tidyr)
#library(raster)
library(sf)
library(exactextractr)
library(stringr)
library(httr2)
#library(httr)
library(terra)

source("~/Documents/Github/blackmarbler/R/blackmarbler.R")

#bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token

roi_sf <- data.frame(lat = -1.943889, lon = 30.059444, id = 1) |>
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) |>
  st_buffer(dist = 20000)

r_20210205 <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A3",
                        date = "2021-02-05",
                        bearer = bearer,
                        quiet = T)

library(leaflet)
r <- r_20210205
# Convert raster to a color palette
palette <- colorNumeric(palette = "viridis", domain = terra::values(r), na.color = "transparent")

# Create Leaflet map and add raster
leaflet() %>%
  addTiles() %>%  # Add base map
  addRasterImage(r, colors = palette, opacity = 0.8) %>%  # Overlay raster
  addLegend(pal = palette, values = terra::values(r), title = "Raster Values")


url1 <- "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A3/2021/032/VNP46A3.A2021032.h21v09.001.2021144024822.h5"
response <- request(url1) %>%
  req_headers('Authorization' = paste('Bearer', bearer)) %>%
  req_timeout(60) %>%
  req_perform() 

# Assuming `response` is your HTTP response object
raw_data <- resp_body_raw(response)  # Get the raw data from the response

# Create a raw connection
raw_conn <- rawConnection(raw_data, open = "rb")

result = base::unserialize(memDecompress(raw_conn))


# Use `rast()` to load the data from the raw connection
h5_data <- rast(raw_data)
h5_data <- rast(raw_conn)
h5_data <- read_stars(raw_data)
h5_data <- read_stars(raw_conn)

library(stars)



close(raw_conn)


h5_data <- rast(resp_body_raw(response))

writeBin(resp_body_raw(response), "~/Desktop/VNP46A3.A2021036.h21v09.001.2021104201220.h5")

h5_file <- "~/Desktop/VNP46A3.A2021036.h21v09.001.2021104201220.h5"
variable <- "AllAngle_Composite_Snow_Free"
quality_flag_rm <- NULL
library(terra)

file_to_raster <- function(h5_file,
                           variable,
                           quality_flag_rm){
  # Converts h5 file to raster.
  # ARGS
  # --f: Filepath to h5 file
  
  # Load data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  h5_data <- terra::rast(h5_file)
  
  # Get extent - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tile_i <- h5_file %>% stringr::str_extract("h\\d{2}v\\d{2}")
  
  bm_tiles_sf <- read_sf("https://raw.githubusercontent.com/worldbank/blackmarbler/main/data/blackmarbletiles.geojson")
  grid_i_sf <- bm_tiles_sf[bm_tiles_sf$TileID %in% tile_i,]
  
  grid_i_sf_box <- grid_i_sf %>%
    st_bbox()
  
  xMin <- min(grid_i_sf_box$xmin) %>% round()
  yMin <- min(grid_i_sf_box$ymin) %>% round()
  xMax <- max(grid_i_sf_box$xmax) %>% round()
  yMax <- max(grid_i_sf_box$ymax) %>% round()
  
  # Load raster of select variable - - - - - - - - - - - - - - - - - - - - - - -
  var_names <- names(h5_data)
  
  if(!(variable %in% var_names)){
    warning(paste0("'", variable, "'",
                   " not a valid variable option. Valid options include:\n",
                   paste(var_names, collapse = "\n")
    ))
  }
  
  out <- h5_data[[variable]]
  
  # Filter by quality - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(length(quality_flag_rm) > 0){
    #### Get quality raster
    if(h5_file %>% str_detect("VNP46A2")){
      qf <- h5_data$Mandatory_Quality_Flag
      
    } else{
      
      variable_short <- variable %>%
        str_replace_all("_Num", "") %>%
        str_replace_all("_Std", "")
      
      qf_name <- paste0(variable_short, "_Quality")
      
      qf <- h5_data[[qf_name]]
    }
    
    #### Filter
    for(val in quality_flag_rm){ # out[qf %in% quality_flag_rm] doesn't work, so loop
      out[qf == val] <- NA
    }
  }
  
  # Add CRS and Extent - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  crs(out) <- "EPSG:4326"
  ext(out) <- c(xMin,xMax,yMin,yMax)
  
  # Additional filtering - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  #set fill values to NA
  out <- remove_fill_value(out, variable)
  
  #apply scaling factor
  out <- apply_scaling_factor(out, variable)
  
  return(out)
}
