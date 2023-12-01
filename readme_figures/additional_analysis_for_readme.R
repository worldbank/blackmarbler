# Additional Analysis for Readme

#### Setup
# Load packages
library(blackmarbler)
library(geodata)
library(sf)
library(raster)
library(ggplot2)

#### Define NASA bearer token
bearer <- "BEARER-TOKEN-HERE"
bearer <- read.csv("~/Desktop/bearer_bm.csv")$token

### ROI
# Define region of interest (roi). The roi must be (1) an sf polygon and (2)
# in the WGS84 (epsg:4326) coordinate reference system. Here, we use the 
# getData function to load a polygon of Ghana
roi_sf <- gadm(country = "NER", level=0, path = tempdir()) |> st_as_sf()

### Daily data: raster for February 5, 2021
r <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A4",
                        date = "2021",
                        bearer = bearer,
               variable = "NearNadir_Composite_Snow_Free")

r_cf <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A2",
               date = "2023-07-23",
               bearer = bearer,
               variable = "Gap_Filled_DNB_BRDF-Corrected_NTL")
r_cf[] <- as.numeric(r_cf[] == 0)
plot(r_cf)


r <- bm_raster(roi_sf = roi_sf,
                  product_id = "VNP46A2",
                  date = "2023-07-23",
                  bearer = bearer,
                  variable = "Gap_Filled_DNB_BRDF-Corrected_NTL",
               quality_flag_rm = c(255))

r[] <- log(r[] + 1)
plot(r)
