# Figures for readme

if(F){
  library(blackmarbler)
  library(readr)
  library(geodata)
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(raster)
  
  # Setup ------------------------------------------------------------------------
  bearer <- read_csv("~/Desktop/bearer_bm.csv") %>%
    pull(token)
  
  roi_sf <- gadm(country = "GHA", level=1, path = tempdir()) %>% st_as_sf()
  product_id <- "VNP46A3"
  year <- 2018
  month <- 5
  day <- 1
  
  # Make map ---------------------------------------------------------------------
  #### Make raster
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A3",
                 date = "2021-10-01",
                 bearer = bearer)
  
  #### Prep data
  r <- r %>% mask(roi_sf) 
  
  r_df <- rasterToPoints(r, spatial = TRUE) %>% as.data.frame()
  names(r_df) <- c("value", "x", "y")
  
  ## Transform NTL
  r_df$value_adj <- log(r_df$value+1)
  
  ##### Map 
  p <- ggplot() +
    geom_raster(data = r_df, 
                aes(x = x, y = y, 
                    fill = value_adj)) +
    scale_fill_gradient2(low = "black",
                         mid = "yellow",
                         high = "red",
                         midpoint = 3.1) +
    labs(title = "Nighttime Lights: October 2021") +
    coord_quickmap() + 
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
  
  ggsave(p,
         filename = file.path("~/Documents/Github/blackmarbler/man/figures/ntl_gha.png"),
         height = 5, width = 6)
  
  # Extract timeseries -----------------------------------------------------------
  ntl_df <- bm_extract(roi_sf = roi_sf,
                       product_id = "VNP46A4",
                       date = 2012:2022,
                       bearer = bearer,
                       aggregation_fun = c("mean"))
  
  p <- ntl_df %>%
    ggplot() +
    geom_col(aes(x = date,
                 y = ntl_mean),
             fill = "darkorange") +
    facet_wrap(~NAME_1) +
    labs(x = NULL,
         y = "NTL Luminosity",
         title = "Ghana Admin Level 1: Annual Average Nighttime Lights") +
    scale_x_continuous(labels = seq(2012, 2022, 4),
                       breaks = seq(2012, 2022, 4)) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold")) 
  
  ggsave(p,
         filename = file.path("~/Documents/Github/blackmarbler/man/figures/ntl_trends_gha.png"),
         height = 5, width = 6)
  
  # Test daily extraction updating -----------------------------------------------------------
  setwd("~/Desktop")
  
  # Create directories to store data
  dir.create(file.path(getwd(), "bm_files"))
  dir.create(file.path(getwd(), "bm_files", "daily"))
  
  bm_extract(roi_sf = roi_sf,
             product_id = "VNP46A2",
             date = seq.Date(from = ymd("2023-01-01"), to = ymd("2023-01-02"), by = 1),
             bearer = bearer,
             output_location_type = "file",
             file_dir = file.path(getwd(), "bm_files", "daily"))
  
  # Append daily-level datasets into one file
  file.path(getwd(), "bm_files", "daily") %>%
    list.files(pattern = "*.Rds",
               full.names = T) %>%
    map_df(readRDS) %>%
    saveRDS(file.path(getwd(), "bm_files", "ntl_daily.Rds"))
  
}