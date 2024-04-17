# Test make raster stack of nighttime lights across multiple time periods

test_that("Test raster stack for VNP46A2 daily data works", {
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data in March 2021
  ken_202103_r <- bm_raster(
    roi_sf = roi_sf,
    product_id = "VNP46A2",
    date = seq.Date(from = lubridate::ymd("2021-03-01"), to = lubridate::ymd("2021-03-03"), by = "day"),
    bearer = bearer,
    quiet = TRUE
  )

  expect_true(class(ken_202103_r) == "SpatRaster",
              info = "ken_202103_r is not a SpatRaster object"
  )
})

test_that("Test raster stack for VNP46A3 monthly data works", {
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Monthly aggregated data in 2021 and 2022
  ken_2021_2022_r <- bm_raster(
    roi_sf = roi_sf,
    product_id = "VNP46A3",
    date = seq.Date(from = lubridate::ymd("2021-01-01"), to = lubridate::ymd("2021-03-01"), by = "month"),
    bearer = bearer,
    quiet = TRUE
  )

  expect_true(class(ken_2021_2022_r) == "SpatRaster",
              info = "ken_2021_2022_r is not a SpatRaster object"
  )
})

test_that("Test raster stack for VNP46A4 anual data works", {
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Yearly aggregated data in 2012 and 2015
  ken_2012_2021_r <- bm_raster(
    roi_sf = roi_sf,
    product_id = "VNP46A4",
    date = 2012:2015,
    bearer = bearer,
    quiet = TRUE
  )

  expect_true(class(ken_2012_2021_r) == "SpatRaster",
              info = "ken_2012_2021_r is not a SpatRaster object"
  )
})

# pending test to run with quality flags
