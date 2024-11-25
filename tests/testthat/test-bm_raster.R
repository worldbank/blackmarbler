test_that("Query VNP46A1", {

  skip()

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Switzerland (covers 2 black marble tiles)
  roi_sf <- geodata::gadm(country = "CHE", level = 0, path = tempdir())

  # Daily data
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A1",
                 date = "2021-10-03",
                 bearer = bearer)

  expect_true(class(r) == "SpatRaster",
              info = "r is not a SpatRaster object"
  )

})

test_that("Query VNP46A2", {

  skip()

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Switzerland (covers 2 black marble tiles)
  roi_sf <- geodata::gadm(country = "CHE", level = 0, path = tempdir())

  # Daily data
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A2",
                 date = "2021-10-03",
                 bearer = bearer)

  expect_true(class(r) == "SpatRaster",
              info = "r is not a SpatRaster object"
  )

})

test_that("Query VNP46A3", {

  skip()

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Switzerland (covers 2 black marble tiles)
  roi_sf <- geodata::gadm(country = "CHE", level = 0, path = tempdir())

  # Daily data
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A3",
                 date = "2021-10",
                 bearer = bearer)

  expect_true(class(r) == "SpatRaster",
              info = "r is not a SpatRaster object"
  )

})

test_that("Query VNP46A4", {

  skip()

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Switzerland (covers 2 black marble tiles)
  roi_sf <- geodata::gadm(country = "CHE", level = 0, path = tempdir())

  # Daily data
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A4",
                 date = 2021,
                 bearer = bearer)

  expect_true(class(r) == "SpatRaster",
              info = "r is not a SpatRaster object"
  )

})

test_that("Get NASA token", {
  
  skip()

  expect_error(get_nasa_token(123, "black"), "username must be a character string")
  expect_error(get_nasa_token("black", 123), "password must be a character string")
  expect_error(get_nasa_token("black", "marble"), "Incorrect username or password")

  username <- Sys.getenv("NASA_USERNAME")
  password <- Sys.getenv("NASA_PASSWORD")

  token <- get_nasa_token(username, password)
  expect_type(token, "character")
  expect_gt(nchar(token), 100)

})


