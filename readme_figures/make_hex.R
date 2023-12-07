# Make Hexagon Logo

if(T){
  
  # Setup ----------------------------------------------------------------------
  library(blackmarbler)
  library(hexSticker)
  library(ggplot2)
  library(tidyverse)
  library(raster)
  library(rgeos)
  library(sp)
  library(sf)
  
  bearer <- read_csv("~/Desktop/bearer_bm.csv") %>%
    pull(token)
  
  # Make ROI -------------------------------------------------------------------
  loc_sf <- data.frame(id = 1,
                       lat = 39.10968435276649, 
                       lon = -84.51499766385844)
  
  coordinates(loc_sf) <- ~lon+lat
  crs(loc_sf) <- CRS("+init=epsg:4326")
  loc_sf <- gBuffer(loc_sf, byid = T, width = 400/111) %>% st_as_sf()
  
  r <- bm_raster(roi_sf = loc_sf,
                 product_id = "VNP46A4",
                 date = 2021,
                 bearer = bearer)
  
  #### Prep data
  #r <- r #%>% mask(roi_sf) 
  
  r_df <- rasterToPoints(r, spatial = TRUE) %>% as.data.frame()
  names(r_df) <- c("value", "x", "y")
  
  ## Remove very low values of NTL; can be considered noise 
  r_df$value[r_df$value <= 3] <- 0
  
  ## Distribution is skewed, so log
  r_df$value_adj <- log(r_df$value+1)
  
  ##### Map 
  p <- ggplot() +
    geom_raster(data = r_df, 
                aes(x = x, y = y, 
                    fill = value_adj)) +
    # scale_fill_gradient2(low = "black",
    #                      mid = "yellow",
    #                      high = "red",
    #                      midpoint = 5) +
    scale_fill_gradient2(low = "#000A33",
                         mid = "yellow",
                         high = "firebrick",
                         midpoint = 6) +
    coord_quickmap() + 
    theme_void() +
    theme(legend.position = "none")
  
  #p
  
  sticker(p, 
          package="blackmarbler", 
          spotlight = F,
          #l_alpha = 1, #0.15,
          p_size=23, #7 
          p_y = 1.40,
          p_family = "sans",
          p_fontface = "italic",
          s_x=1, 
          s_y=0.9, 
          s_width=2.5, 
          s_height=2.5,
          p_color = "white",
          h_fill = "black",
          h_color = "black",
          white_around_sticker = T,
          l_y = 1.4,
          l_x = 0.93,
          l_width = 3,
          l_height = 3,
          filename="~/Documents/Github/blackmarbler/man/figures/hex.png")
  
  # sticker(p, 
  #         package="blackmarblepy", 
  #         spotlight = F,
  #         #l_alpha = 1, #0.15,
  #         p_size=20, #7 
  #         p_y = 1.40,
  #         p_family = "sans",
  #         p_fontface = "italic",
  #         s_x=1, 
  #         s_y=0.9, 
  #         s_width=2.5, 
  #         s_height=2.5,
  #         p_color = "white",
  #         h_fill = "black",
  #         h_color = "black",
  #         white_around_sticker = T,
  #         l_y = 1.4,
  #         l_x = 0.93,
  #         l_width = 3,
  #         l_height = 3,
  #         filename="~/Documents/Github/blackmarblepy/docs/images/logo.png")
  
 
}