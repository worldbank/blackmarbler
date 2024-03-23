library(sf)

black_marble_tiles_sf <- read_sf("https://raw.githubusercontent.com/worldbank/blackmarbler/main/data/blackmarbletiles.geojson")

usethis::use_data(black_marble_tiles_sf,
                  overwrite = T,
                  internal = T)
