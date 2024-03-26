# Test for daily data
test_that("Test raster for VNP46A2 works", {
  # Define bearer token
  bearer <- Sys.getenv('BEARER_NASA_TOKEN')

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |>  sf::st_as_sf()

  ken_20210205_r <- bm_raster(roi_sf = roi_sf,
                              product_id = "VNP46A2",
                              date = "2021-10-03",
                              bearer = bearer)
  expect_true(class(ken_20210205_r)[1] == "SpatRaster",
              info = "ken_20210205_r is not a SpatRaster object")
})

# Test for monthly data
test_that("Test raster for VNP46A3 monthly data works", {
  # Define bearer token
  bearer <- Sys.getenv('BEARER_NASA_TOKEN')

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |>  sf::st_as_sf()

  ken_202103_r <- bm_raster(roi_sf = roi_sf,
                            product_id = "VNP46A3",
                            date = "2021-03-01",
                            bearer = bearer)

  expect_true(class(ken_202103_r) == "SpatRaster",
              info = "ken_202103_r is not a SpatRaster object")
})

# Test for annual data
test_that("Test raster for VNP46A4 anual data works", {
  # Define bearer token
  bearer <- Sys.getenv('BEARER_NASA_TOKEN')

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |>  sf::st_as_sf()

  ken_2021_r <- bm_raster(roi_sf = roi_sf,
                          product_id = "VNP46A4",
                          date = 2021,
                          bearer = bearer)
  expect_true(class(ken_2021_r) == "SpatRaster",
              info = "ken_2021_r is not a SpatRaster object")
})
