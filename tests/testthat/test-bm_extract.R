test_that("Test extract for VNP46A2 works", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

# Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                            product_id = "VNP46A2",
                            date = "2021-10-03",
                            bearer = bearer,
                            quiet = TRUE)

expect_true(class(ken_202103_r) == "data.frame",
            info = "ken_202103_r is not a data.frame object"
)

})

test_that("Test extract for VNP46A3 works.. fails", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202303_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A3",
                             date = "2023-03-01",
                             bearer = bearer,
                             quiet = TRUE,
                             check_all_tiles_exist = TRUE)

  expect_true(class(ken_202303_r) == "data.frame",
              info = "ken_202303_r is not a data.frame object"
  )

})

test_that("Test extract for VNP46A4 works", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A4",
                             date = 2021,
                             bearer = bearer,
                             quiet = TRUE)

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})

