#### Setup
# Load packages
#library(blackmarbler)
library(geodata)
library(sf)
library(raster)
library(ggplot2)

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
library(httr)

#devtools::install_github("worldbank/blackmarbler@new-version")

#### Define NASA bearer token
bearer <- "BEARER-TOKEN-HERE"
bearer <- read.csv("~/Desktop/bearer_bm.csv")$token

### ROI
# Define region of interest (roi). The roi must be (1) an sf polygon and (2)
# in the WGS84 (epsg:4326) coordinate reference system. Here, we use the 
# getData function to load a polygon of Ghana
#roi_sf <- gadm(country = "GHA", level=1, path = tempdir()) |> st_as_sf()

#########
roi_sf <- gadm(country = "CHE", level=1, path = tempdir()) |> st_as_sf()

########
r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = c("2021-01-01"),
               bearer = bearer,
               quiet = T)

r <- bm_extract(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = c("2021-01-01"),
               bearer = bearer,
               quiet = F)

r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = c("2021-01-01",
                        "2021-02-01",
                        "2021-03-01"),
               bearer = bearer,
               interpol_na = T,
               method = "linear",
               quality_flag_rm = c(255,2))

bm_r <- raster::approxNA(r)

ro_df <- bm_extract(roi_sf = roi_sf,
                   product_id = "VNP46A3",
                   date = c("2021-01-01",
                            "2021-02-01",
                            "2021-03-01"),
                   bearer = bearer,
                   quality_flag_rm = c(255,2))

r_df <- bm_extract(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = c("2021-01-01",
                        "2021-02-01",
                        "2021-03-01"),
               bearer = bearer,
               interpol_na = T,
               method = "linear",
               quality_flag_rm = c(255,2))

r_df <- r_df %>%
  arrange(date, NAME_1)

ro_df <- ro_df %>%
  arrange(date, NAME_1)

count_n_obs <- function(values, coverage_fraction) {
  
  orig_vars <- names(values)
  
  values %>%
    dplyr::mutate(across(orig_vars, ~ as.numeric(!is.na(.)) )) %>%
    dplyr::summarise(across(orig_vars, sum, .names = "n_non_na_pixels.{.col}"),
                     across(orig_vars, ~length(.), .names = "n_pixels.{.col}"))
}

roi_df <- roi_sf %>% st_drop_geometry()
roi_df$date <- NULL

n_obs_df <- exact_extract(r, roi_sf, count_n_obs) %>%
  bind_cols(roi_df) %>%
  tidyr::pivot_longer(cols = -c(names(roi_df)),
               names_to = c(".value", "date"),
               names_sep = "\\.t") %>%
  dplyr::mutate(prop_non_na_pixels = n_non_na_pixels / n_pixels)


df <- exact_extract(r, roi_sf, "mean")

df %>%
  pivot_longer(cols = everything(),
               names_to = c(".value", "date"),
               names_sep = "\\.t")





r_n_obs <- exact_extract(r, roi_sf, function(values, coverage_fraction)
  sum(!is.na(values)))

r_n_obs_poss <- exact_extract(r, roi_sf, function(values, coverage_fraction)
  length(values))

roi_sf$n_pixels           <- r_n_obs_poss
roi_sf$n_non_na_pixels    <- r_n_obs
roi_sf$prop_non_na_pixels <- roi_sf$n_non_na_pixels / roi_sf$n_pixels 

plot(r)

#########

### Daily data: raster for February 5, 2021
r_20210205 <- bm_raster(roi_sf = roi_sf,
                        product_id = "VNP46A2",
                        date = "2021-02-05",
                        bearer = bearer)

### Monthly data: raster for October 2021
r_202110 <- bm_raster(roi_sf = roi_sf,
                      product_id = "VNP46A3",
                      date = "2021-10-01", # The day is ignored
                      bearer = bearer)

### Annual data: raster for 2021
r_2021 <- bm_raster(roi_sf = roi_sf,
                    product_id = "VNP46A4",
                    date = 2021,
                    bearer = bearer)

#### Daily data in March 2021
r_daily <- bm_raster(roi_sf = roi_sf,
                     product_id = "VNP46A3",
                     date = seq.Date(from = ymd("2021-03-01"), to = ymd("2021-03-31"), by = "day"),
                     bearer = bearer)

#### Monthly aggregated data in 2021 and 2022
r_monthly <- bm_raster(roi_sf = roi_sf,
                       product_id = "VNP46A3",
                       date = seq.Date(from = ymd("2021-01-01"), to = ymd("2022-12-01"), by = "month"),
                       bearer = bearer)

#### Yearly aggregated data in 2012 and 2021
r_annual <- bm_raster(roi_sf = roi_sf,
                      product_id = "VNP46A4",
                      date = 2012:2021,
                      bearer = bearer)

#### Make raster
r <- bm_raster(roi_sf = roi_sf,
               product_id = "VNP46A3",
               date = "2021-10-01",
               bearer = bearer)

#### Prep data
r <- r |> mask(roi_sf) 

r_df <- rasterToPoints(r, spatial = TRUE) |> as.data.frame()
names(r_df) <- c("value", "x", "y")

## Remove very low values of NTL; can be considered noise 
r_df$value[r_df$value <= 2] <- 0

## Distribution is skewed, so log
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
  labs(title = "NTL, October 2021") +
  coord_quickmap() + 
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

#### Extract annual data
ntl_df <- bm_extract(roi_sf = roi_sf,
                     product_id = "VNP46A4",
                     date = 2012:2022,
                     bearer = bearer)

#### Trends over time
ntl_df |>
  ggplot() +
  geom_col(aes(x = date,
               y = ntl_mean),
           fill = "darkorange") +
  facet_wrap(~NAME_1) +
  labs(x = NULL,
       y = "NTL Luminosity",
       title = "Ghana Admin Level 1: Annual Average Nighttime Lights") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) 

