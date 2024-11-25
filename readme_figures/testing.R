# Testing

library(readr)
library(hdf5r)
library(dplyr)
library(purrr)
library(lubridate)
library(tidyr)
library(raster)
library(sf)
library(exactextractr)
library(stringr)
library(httr2)
library(httr)

source("~/Documents/Github/blackmarbler/R/blackmarbler.R")

bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token

roi_sf <- data.frame(lat = -1.943889, lon = 30.059444, id = 1) |>
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) |>
  st_buffer(dist = 20000)

r_20210205 <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A2",
                        date = "2021-02-05",
                        bearer = bearer,
                        quiet = T)

r_20210205 <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A1",
                        date = "2021-02-05",
                        variable = "Sensor_Zenith",
                        bearer = bearer)


bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token
bearer <- "BEARER HERE"
#library(blackmarbler)
library(geodata)

username <- readline("Username")
password <- readline("Password")
bearer <- get_nasa_token(username, password)

out <- h5_data[[paste0("HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/", variable)]][,]
qf  <- h5_data[["HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/Mandatory_Quality_Flag"]][,]

h5_data[["HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/QF_Cloud_Mask"]][,] %>% table()

h5_file <- "~/Desktop/VNP46A1.A2021036.h20v09.001.2021037082509.h5"

h5_data <- h5file(h5_file, "r+")

# Setup ------------------------------------------------------------------------
source("~/Documents/Github/blackmarbler/R/blackmarbler.R")
library(blackmarbler)
library(geodata)

bearer <- read.csv("~/Desktop/bearer_bm.csv")$token

roi_sf <- gadm(country = "CHE", level=1, path = tempdir()) 

# bm_raster: Basic -------------------------------------------------------------
r1 <- bm_raster(roi_sf = roi_sf,
                product_id = "VNP46A2",
                date = "2021-10-01", 
                bearer = bearer)

r2 <- bm_raster(roi_sf = roi_sf,
                product_id = "VNP46A3",
                date = "2021-10-01", 
                bearer = bearer)

r3 <- bm_raster(roi_sf = roi_sf,
                product_id = "VNP46A4",
                date = "2021", 
                bearer = bearer)






dir.create("~/Desktop/h5_tmp")
r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = "2021-10-01", 
               bearer = bearer,
               h5_dir = "~/Desktop/h5_tmp")

dir.create("~/Desktop/ntl_tmp")
r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = "2021-10-01", 
               bearer = bearer,
               h5_dir = "~/Desktop/h5_tmp",
               output_location_type = "file",
               file_dir = "~/Desktop/ntl_tmp",
               file_return_null = T)

r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A4",
               date = 2018:2020, 
               bearer = bearer,
               h5_dir = "~/Desktop/h5_tmp",
               interpol_na = T)

exact_extract(r, st_as_sf(roi_sf), fun = "mean")

df <- bm_extract(roi_sf = roi_sf,
                 product_id = "VNP46A4",
                 date = 2018:2020, 
                 bearer = bearer,
                 h5_dir = "~/Desktop/h5_tmp",
                 output_location_type = "file",
                 file_dir = "~/Desktop/ntl_tmp",
                 file_return_null = F,
                 aggregation_fun = c("mean", "sum"))

bm_r <- terra::approximate(r,
                           method = method,
                           rule = rule,
                           f = f,
                           ties = ties,
                           z = z,
                           NArule = NArule)

r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A4",
               date = 2015:2020, 
               bearer = bearer,
               h5_dir = "~/Desktop/h5_tmp",
               output_location_type = "file",
               file_dir = "~/Desktop/ntl_tmp")

# bm_extract: Basic ------------------------------------------------------------
ntl_df <- bm_extract(roi_sf = roi_sf,
                     product_id = "VNP46A3",
                     date = "2021-10-01", 
                     bearer = bearer)

ntl_df <- bm_extract(roi_sf = roi_sf,
                     product_id = "VNP46A4",
                     date = 2020, 
                     bearer = bearer)

r_202110 <- bm_raster(roi_sf = roi_sf,
                      product_id = "VNP46A3",
                      variable = "NearNadir_Composite_Snow_Free",
                      date = "2021-10-01", 
                      bearer = bearer,
                      h5_dir = "~/Desktop/h5_test",
                      quiet = T,
                      output_location_type = "file",
                      file_dir = "~/Desktop")

r_202110 <- bm_extract(roi_sf = roi_sf,
                       product_id = "VNP46A3",
                       variable = "NearNadir_Composite_Snow_Free",
                       date = "2021-10-01", 
                       bearer = bearer,
                       h5_dir = "~/Desktop/h5_test",
                       output_location_type = "file",
                       file_dir = "~/Desktop")

a <- terra::extract(r_202110, roi_sf, fun = sum, exact = T)$t2021_10
b <- exact_extract(r_202110, roi_sf, fun = "sum")

plot(r_202110)
plot(roi_sf,add=T)

e_202110 <- bm_raster(roi_sf = roi_sf,
                      product_id = "VNP46A3",
                      date = c("2021-10-01", "2021-11-01"), 
                      bearer = bearer,
                      output_location_type = "file",
                      file_dir = "~/Desktop/test1")

e_202110 <- bm_extract(roi_sf = roi_sf,
                       product_id = "VNP46A3",
                       date = c("2021-10-01", "2021-11-01"), 
                       bearer = bearer)


