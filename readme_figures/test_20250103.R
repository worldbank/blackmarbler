# Test Package

# INSTALL DEV VERSION FROM GITHUB
# devtools::install_github("worldbank/blackmarbler", auth_token = "")

library(blackmarbler)
library(sf)
library(terra)
library(httr2)

# ADD BEARER TOKEN HERE
bearer <- read.csv("~/Dropbox/bearer_bm.csv")$token

# CHECK 1 ----------------------------------------------------------------------
# With sessionInfo, under "other attached packages", make sure the blackmarbler
# version is: blackmarbler_0.2.4
sessionInfo()

# CHECK 2 ----------------------------------------------------------------------
# What is the response? In particular, what is the:
# - Status
# - Content-Type
# - Body (how many bytes?)
url1 <- "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A3/2021/032/VNP46A3.A2021032.h21v09.001.2021144024822.h5"
response <- request(url1) %>%
  req_headers('Authorization' = paste('Bearer', bearer)) %>%
  req_timeout(60) %>%
  req_perform() 

response

length(response$body) < 10000

# CHECK 3 ----------------------------------------------------------------------
roi_sf <- data.frame(lat = -1.943889, lon = 30.059444, id = 1) |>
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) |>
  st_buffer(dist = 20000)

r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A2",
               date = "2018-09-01",
               bearer = bearer)

plot(r)


