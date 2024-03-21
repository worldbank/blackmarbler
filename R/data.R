#' Black Marble Tiles Spatial Data
#'
#' Spatial dataset containing black marble tiles obtained from NASA's Black Marble project.
#'
#' @format A spatial dataframe (\code{sf}) with 648 rows and 2 columns:
#' \describe{
#'   \item{TileID}{Identification of the tile associated with the VIIRS satellites mapping.}
#'   \item{geometry}{Spatial geometry representing the tile boundaries.}
#' }
#'
#' @source
#' \url{https://blackmarble.gsfc.nasa.gov/}
#'
#' The dataset is obtained from the Black Marble project by NASA,
#' specifically the VIIRS (Visible Infrared Imaging Radiometer Suite) satellites.
#' The TileID column represents the identification of the tile associated with the
#' VIIRS satellites mapping. For more details on the dataset and its
#' characteristics, refer to the Black Marble User Guide:
#' \url{https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf}
#'
#' @references
#' Tan, Bin, et al. "Night-time light imagery analysis for urban population
#' mapping in China." International Journal of Remote Sensing 27.9 (2006): 1943-1950.
#'
#' Wolfe, Robert E., et al. "Achieving accuracy requirements for forest
#' biomass mapping: A data fusion method for combining the strengths of
#' spaceborne lidar and Landsat TM data." Remote Sensing of Environment
#' 80.2 (2002): 253-267.
#'
#' @seealso
#' \url{https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20220916.pdf}
#'
#' @examples
#' library(sf)
#' black_marble_tiles_sf <- read_sf("https://raw.githubusercontent.com/worldbank/blackmarbler/main/data/blackmarbletiles.geojson")
#' head(black_marble_tiles_sf)
"black_marble_tiles_sf"
