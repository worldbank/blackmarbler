
# Define the test cases
test_that("julian_to_month returns correct month for given Julian dates", {
  expect_equal(julian_to_month("001"), 1)
  expect_equal(julian_to_month("032"), 2)
  expect_equal(julian_to_month("060"), 3)
  expect_equal(julian_to_month("091"), 4)
  expect_equal(julian_to_month("121"), 5)
  expect_equal(julian_to_month("152"), 6)
  expect_equal(julian_to_month("182"), 7)
  expect_equal(julian_to_month("213"), 8)
  expect_equal(julian_to_month("244"), 9)
  expect_equal(julian_to_month("274"), 10)
  expect_equal(julian_to_month("305"), 11)
  expect_equal(julian_to_month("335"), 12)
})



# Define the unit test
test_that("remove_fill_value_from_satellite_data correctly removes artifact values", {
  # Create a sample dataset for testing
  sample_data <- c(255,
                   -999.9,
                   -32768,
                   65535)

  # Test for removal of artifact value 255
  cleaned_data_255 <- remove_fill_value_from_satellite_data(sample_data[1], "Granule")

  expect_equal(sum(is.na(cleaned_data_255)), 1,
               info =     "Expected 1 NA value after removing artifact value 255")

  # Test for removal of artifact value -999.9
  cleaned_data_999 <- remove_fill_value_from_satellite_data(sample_data[2], "UTC_Time")
  expect_equal(sum(is.na(cleaned_data_999)), 1,
               info =     "Expected 1 NA value after removing artifact value -999.9")

  # Test for removal of artifact value -32768
  cleaned_data_32768 <- remove_fill_value_from_satellite_data(sample_data[3], "Sensor_Azimuth")

  expect_equal(sum(is.na(cleaned_data_32768)), 1,
               info =  "Expected 1 NA value after removing artifact value -32768")

  # Test for removal of artifact value 65535
  cleaned_data_65535 <- remove_fill_value_from_satellite_data(sample_data[4],
                                                              "BrightnessTemperature_M12")

  expect_equal(sum(is.na(cleaned_data_65535)), 1,
               info =   "Expected 1 NA value after removing artifact value 65535")
})


# Define the unit test
test_that("apply_scaling_factor_to_viirs_data correctly applies scaling factor", {
  # Create a sample dataset for testing

  sample_data <- 2024

  # Test for scaling factor applied to DNB_At_Sensor_Radiance variable

  scaled_data <- apply_scaling_factor_to_viirs_data(sample_data,
                                                    "DNB_At_Sensor_Radiance",
                                                    quiet = FALSE)

  expect_equal(scaled_data, 202.4,
               info = "Expected scaled value for DNB_At_Sensor_Radiance to be 202.4")

  # Test for scaling factor applied to DNB_BRDF-Corrected_NTL variable
  scaled_data <- apply_scaling_factor_to_viirs_data(sample_data,
                                                    "DNB_BRDF-Corrected_NTL",
                                                    quiet = FALSE)

  expect_equal(scaled_data, 202.4,
               info = "Expected scaled value for DNB_BRDF-Corrected_NTL to be 202.4")

  # Test for scaling factor applied to AllAngle_Composite_Snow_Covered variable

  scaled_data <- apply_scaling_factor_to_viirs_data(sample_data, "AllAngle_Composite_Snow_Covered",
                                                    quiet = FALSE)


  expect_equal(scaled_data, 202.4,
               info = "Expected scaled value for AllAngle_Composite_Snow_Covered to be 202.4")
})

