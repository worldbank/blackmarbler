# Testing

# devtools::install_github("worldbank/blackmarbler", auth_token = "")

library(blackmarbler)
library(sf)
library(terra)

# Using package ----------------------------------------------------------------
bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token

roi_sf <- data.frame(lat = -1.943889, lon = 30.059444, id = 1) |>
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) |>
  st_buffer(dist = 20000)

r_20210205 <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A2",
                        date = seq.Date(ymd("2024-01-29"),
                                        ymd("2024-02-01"),
                                        1),
                        bearer = bearer)



