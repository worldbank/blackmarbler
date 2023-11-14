
#remove.packages("blackmarbler")
#devtools::install_github("worldbank/blackmarbler")
library(purrr)
#library(furrr)
library(stringr)
#library(rhdf5)
library(raster)
library(dplyr)
library(sf)
library(lubridate)
library(geodata)
library(exactextractr)

#library(blackmarbler)
library(geodata)
library(sf)
library(dplyr)
library(httr)

library(hdf5r)
source("~/Documents/Github/blackmarbler/R/blackmarbler.R")

roi <- gadm(country = "GHA", 
            level=0, 
            path = tempdir()) |> 
  st_as_sf()

bearer <- read.csv("~/Desktop/bearer_bm.csv") %>%
  pull(token)

r <- bm_raster(roi_sf = roi,
                product_id = "VNP46A2",
                date = "2020-01-01",
                bearer = bearer)




#bearer <- "BEARER-TOKEN"





r <- r |> mask(roi_sf) 

r_df <- rasterToPoints(r, spatial = TRUE) |> as.data.frame()
names(r_df) <- c("value", "x", "y")

## Transform NTL
r_df$value[r_df$value <= 1] <- 0

r_df$value_adj <- log(r_df$value+1)

##### Map 
p <- ggplot() +
  geom_raster(data = r_df, 
              aes(x = x, y = y, 
                  fill = value_adj)) +
  scale_fill_gradient2(low = "black",
                       mid = "yellow",
                       high = "red",
                       midpoint = 4.5) +
  coord_quickmap() + 
  theme_void() +
  theme(legend.position = "none")





bearer <- read.csv("~/Desktop/bearer_bm.csv") %>%
  pull(token)

gha_0_sf <- gadm(country = "GHA", level=0, path = tempdir()) %>% st_as_sf()

r_2021 <- bm_raster(roi_sf = gha_0_sf,
                    product_id = "VNP46A4",
                    date = 2018,
                    bearer = bearer)

r_2021_c <- bm_raster(roi_sf = gha_0_sf,
                    product_id = "VNP46A4",
                    date = 2020,
                    variable = "NearNadir_Composite_Snow_Free",
                    bearer = bearer)

r_2021_c_n <- bm_raster(roi_sf = gha_0_sf,
                      product_id = "VNP46A3",
                      date = "2021-10",
                      variable = "NearNadir_Coamposite_Snow_Free_Num",
                      bearer = bearer)

h5_data <- H5Fopen("~/Desktop/testaa.h5")

variable <- "NearNadir_Composite_Snow_Free_Num"
out <- h5_data$HDFEOS$GRIDS$VIIRS_Grid_DNB_2d$`Data Fields`[[variable]]

summary(as.vector(out))

headers <- c('Authorization' = paste('Bearer', bearer))
response <- GET(urla, add_headers(headers), write_disk("~/Desktop/aaaa.h5", overwrite = TRUE))
response

response <- GET(urla, add_headers(headers), write_disk(download_path, overwrite = TRUE))




library(rnaturalearth)
library(rgeos)
c_df <- ne_countries(type = "countries", scale = "medium", returnclass = "sf")
c_df <- c_df[c_df$continent != "Antarctica",]
c_df <- c_df[c_df$continent != "Seven seas (open ocean)",]
c_df <- c_df %>% gBuffer(width = 0, byid = T)
c_df <- c_df %>% st_as_sf()
c_df <- st_buffer(c_df, 0)

bm_tiles_sf <- bm_tiles_sf[!(bm_tiles_sf$TileID %>% str_detect("h00")),]
bm_tiles_sf <- bm_tiles_sf[!(bm_tiles_sf$TileID %>% str_detect("v00")),]

st_intersects(bm_tiles_sf, c_df, sparse = F)

r_ntl <- bm_raster(roi_sf = gha_1_sf,
                   product_id = "VNP46A4",
                   date = 2021,
                   bearer = bearer,
                   quiet = F)

r_ntl <- bm_raster(roi_sf = gha_1_sf,
                   product_id = "VNP46A2",
                   date = "2021-01-01",
                   bearer = bearer)

r_ntl <- bm_raster(roi_sf = gha_1_sf,
                   product_id = "VNP46A4",
                   date = 2021,
                   bearer = bearer)

ntl_df <- bm_extract(roi_sf = gha_1_sf,
                     product_id = "VNP46A3",
                     date = "2021-01-01",
                     bearer = bearer)

r_ntl <- bm_raster(roi_sf = gha_1_sf,
                   product_id = "VNP46A3",
                   date = ymd("2020-01-01"),
                   bearer = bearer,
                   variable = "NearNadir_Composite_Snow_Free")

r_ntl2 <- bm_raster(roi_sf = gha_1_sf,
                    product_id = "VNP46A3",
                    date = ymd("2020-01-01"),
                    bearer = bearer,
                    variable = "NearNadir_Composite_Snow_Free_Quality")


roi_sf = gha_1_sf
product_id = "VNP46A3"
date = ymd("2020-01-01")
bearer = bearer
variable = "NearNadir_Composite_Snow_Free_Quality"

output_location_type = "r_memory"
file_dir = NULL
file_prefix = NULL
file_skip_if_exists = TRUE


r <- bm_raster(roi_sf = gha_1_sf,
               product_id = "VNP46A4",
               date = 2019:2020,
               bearer = bearer,
               output_type = "aggregation",
               output_location_type = "file",
               aggregation_fun = c("mean", "median"),
               file_dir = "~/Desktop/bmtest",
               file_prefix = NULL,
               file_skip_if_exists = TRUE)

r <- bm_raster(roi_sf = gha_1_sf,
               product_id = "VNP46A4",
               date = 2019:2020,
               bearer = bearer,
               output_type = "aggregation",
               output_location_type = "file",
               aggregation_fun = c("mean", "median"),
               file_dir = "~/Desktop/bmtest",
               file_prefix = NULL,
               file_skip_if_exists = TRUE)


r <- bm_raster(roi_sf = gha_1_sf,
               product_id = "VNP46A4",
               date = 2019:2020,
               bearer = bearer,
               output_type = "aggregation",
               output_location_type = "file",
               aggregation_fun = c("mean", "median"),
               file_dir = "~/Desktop/bmtest",
               file_prefix = NULL,
               file_skip_if_exists = TRUE)


