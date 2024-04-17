
# Multiple vs single dates ------------------------------------------------

test_that("Test extract for VNP46A2 pasing multiple dates works", {
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data in March 2021
  ken_202103_r <- bm_extract(
    roi_sf = roi_sf,
    product_id = "VNP46A2",
    date = seq.Date(from = lubridate::ymd("2021-03-01"), to = lubridate::ymd("2021-03-03"), by = "day"),
    bearer = bearer,
    quiet = TRUE
  )

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )
})

test_that("Test extract for VNP46A3 pasing multiple dates works", {
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Monthly aggregated data in 2021 and 2022
  ken_202103_r <- bm_extract(
    roi_sf = roi_sf,
    product_id = "VNP46A3",
    date = seq.Date(from = lubridate::ymd("2021-01-01"), to = lubridate::ymd("2021-03-01"), by = "month"),
    bearer = bearer,
    quiet = TRUE
  )

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )
})



test_that("Test extract for VNP46A4 pasing multiple dates works", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A4",
                             date = 2012:2015,
                             bearer = bearer,
                             quiet = TRUE)

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})
