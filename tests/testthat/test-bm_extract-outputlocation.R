# Output location (memory vs file), checking both with multiple da --------

# * Output location (memory vs file), checking both with multiple dates. (eg, for "file" and multiple dates, outputs separate files)
test_that("Test output location for 1 date works", {

  temp_location <- fs::dir_create("output_location_tests")

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A2",
                             date = "2021-10-03",
                             output_location = "file",
                             file_dir = temp_location,
                             file_prefix = NULL, # add warnign of requeired arguemtns
                             bearer = bearer,
                             quiet = TRUE)

  files_in_temp <- fs::dir_ls(temp_location)


  # test length is 2
  expect_equal(length(files_in_temp), 2,
               info = "Length of files in temp location is not 2"
  )

  # test .tif and .Rds files exist
  expect_true(any(fs::path_ext(files_in_temp) == "tif"),
              info = "No .tif files in temp location"
  )
  expect_true(any(fs::path_ext(files_in_temp) == "Rds"),
              info = "No .rds files in temp location"
  )

  fs::dir_delete(temp_location)

})

test_that("Test extract for multiple dates works", {

  temp_location <- fs::dir_create("output_location_tests")

  # Define bearer token
  bearer <- Sys.getenv("BEARER_NASA_TOKEN")

  # sf polygon of Ghana
  roi_sf <- geodata::gadm(country = "GHA", level = 1, path = tempdir()) |> sf::st_as_sf()

  dates_to_run <- seq.Date(from = lubridate::ymd("2021-03-01"), to = lubridate::ymd("2021-03-03"), by = "day")

  # Daily data: raster for October 3, 2021
  ken_202103_r <- bm_extract(roi_sf = roi_sf,
                             product_id = "VNP46A2",
                             date = dates_to_run,
                             output_location = "file",
                             file_dir = temp_location,
                             file_prefix = NULL, # add warnign of requeired arguemtns
                             bearer = bearer,
                             quiet = TRUE)

  # test length is 2
  expect_equal(length(files_in_temp), length(dates_to_run)*2,
               info = "Length of files in temp location is incorrect"
  )

  # test .tif and .Rds files exist
  expect_true(any(fs::path_ext(files_in_temp) == "tif"),
              info = "No .tif files in temp location"
  )
  expect_true(any(fs::path_ext(files_in_temp) == "Rds"),
              info = "No .rds files in temp location"
  )

  fs::dir_delete(temp_location)

})
