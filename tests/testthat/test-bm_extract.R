test_that("multiplication works", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

# Daily data: raster for October 3, 2021
ken_20210205_r <- bm_extract(roi_sf = roi_sf,
                            product_id = "VNP46A2",
                            date = "2021-10-03",
                            bearer = bearer)
})
