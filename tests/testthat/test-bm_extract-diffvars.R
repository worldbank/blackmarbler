
# Changing the variable from black marble (using another NTL varia --------

# Changing the variable from black marble (using another NTL variable), or something like cloud cover [changing using the "variable" param]
test_that("Test extract for cloud cover variable works on VNP46A2", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A2",
                             variable = "QF_Cloud_Mask",
                             date = "2021-10-03",
                             bearer = bearer,
                             quiet = TRUE)

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})

test_that("Test extract for AllAngle_Composite_Snow_Free variable works on VNP46A3", {
  #https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.3_Sep_2022.pdf

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  ken_202303_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A3",
                             variable = "AllAngle_Composite_Snow_Free", # no var  named AllAngle_Composite_Snow_Free
                             date = "2023-03-01",
                             bearer = bearer,
                             quiet = TRUE,
                             check = TRUE)

  expect_true(class(ken_202303_r) == "data.frame",
              info = "ken_202303_r is not a data.frame object"
  )

})

test_that("Test extract for AllAngle_Composite_Snow_Free variable works on VNP46A4", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A4",
                             variable = "AllAngle_Composite_Snow_Free",
                             date = 2022,
                             bearer = bearer,
                             quiet = TRUE)

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})

test_that("Test extract for snow cover variable works on VNP46A4", {

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A4",
                             variable = "AllAngle_Composite_Snow_Covered",
                             date = 2022,
                             bearer = bearer,
                             quiet = TRUE)

  expect_true(class(ken_202103_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})



test_that("Test extract for AllAngle_Composite_Snow_Covered variable works on VNP46A3", {
  #https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/archives/Document%20Archive/Science%20Data%20Product%20Documentation/VIIRS_Black_Marble_UG_v1.3_Sep_2022.pdf

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202303_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A3",
                             variable = "AllAngle_Composite_Snow_Covered", # no var  named QF_Cloud_Mask
                             date = "2023-03-01",
                             bearer = bearer,
                             quiet = TRUE,
                             check = TRUE)

  expect_true(class(ken_202303_r) == "data.frame",
              info = "ken_202103_r is not a data.frame object"
  )

})
