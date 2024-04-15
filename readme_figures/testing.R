# Testing

# Setup ------------------------------------------------------------------------
library(geodata)
library(sf)
library(terra)
library(ggplot2)

library(readr)
library(hdf5r)
library(dplyr)
library(purrr)
library(lubridate)
library(tidyr)
library(sf)
library(exactextractr)
library(stringr)
library(httr)

bearer <- read.csv("~/Desktop/bearer_bm.csv")$token

roi_sf <- gadm(country = "CHE", level=1, path = tempdir()) |> st_as_sf()

roi_sf = roi_sf
product_id = "VNP46A3"
date = "2018-04"
bearer = bearer
variable = "AllAngle_Composite_Snow_Free"
quality_flag_rm = NULL
check_all_tiles_exist = TRUE
interpol_na = FALSE
output_location_type = "memory"
file_dir = NULL
file_prefix = NULL
file_skip_if_exists = TRUE
quiet = FALSE

r_202110 <- bm_raster(roi_sf = roi_sf,
                      product_id = "VNP46A3",
                      date = "2021-10-01", 
                      bearer = bearer)

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


