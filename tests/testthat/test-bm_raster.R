test_that("Query VNP46A1", {
  
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
  
  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")
  
  # sf polygon of Switzerland (covers 2 black marble tiles)
  roi_sf <- geodata::gadm(country = "CHE", level = 0, path = tempdir())
  
  # Daily data
  r <- bm_raster(roi_sf = roi_sf,
                 product_id = "VNP46A3",
                 date = 2021,
                 bearer = bearer)
  
  expect_true(class(r) == "SpatRaster",
              info = "r is not a SpatRaster object"
  )
  
})

