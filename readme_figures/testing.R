# Testing

library(blackmarbler)
library(sf)
library(geodata)
library(exactextractr)

# Using package ----------------------------------------------------------------
bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token

roi_sf <- data.frame(lat = -1.943889, lon = 30.059444, id = 1) |>
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) |>
  st_buffer(dist = 20000)

sudan_sf <- gadm("SDN", level = 2, path = tempdir()) %>%
  st_as_sf()

r24 <- bm_raster(roi_sf = sudan_sf,
                        product_id = "VNP46A4",
                        date = 2024,
                        bearer = bearer)
r19 <- bm_raster(roi_sf = sudan_sf,
                 product_id = "VNP46A4",
                 date = 2019,
                 bearer = bearer)
r23 <- bm_raster(roi_sf = sudan_sf,
                 product_id = "VNP46A4",
                 date = 2023,
                 bearer = bearer)
sudan_sf$ntl24 <- exact_extract(r24, sudan_sf, "sum")
sudan_sf$ntl23 <- exact_extract(r23, sudan_sf, "sum")
sudan_sf$ntl19 <- exact_extract(r19, sudan_sf, "sum")

r23 <- r23 %>% mask(sudan_sf[sudan_sf$NAME_2 == "Al Fasher",]) %>% crop(sudan_sf[sudan_sf$NAME_2 == "Al Fasher",])
table(r23[] > 0)
r23[][r23[] > 0] %>% summary()

r <- log(r)
r_ll <- r
pal <- colorNumeric(
  palette = "viridis",
  domain  = terra::values(r_ll),
  na.color = "transparent"
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(
    r_ll,
    colors  = pal,
    opacity = 0.8
  ) %>%
  addLegend(
    pal     = pal,
    values  = terra::values(r_ll),
    title   = "Raster value"
  )



Al-Fasher