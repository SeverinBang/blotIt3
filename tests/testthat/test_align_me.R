# align_me() - model parsing ----------------------------------------------

test_that("align_me() - model parsing", {

  sim_data_wide_file <- system.file(
    "extdata", "sim_data_wide.csv",
    package = "blotIt3"
  )

  input_data <- read_wide(sim_data_wide_file, description = seq_len(3))

  expect_error(
    align_me(data = input_data, fixed = "one_parameter"),
    "Left and right-hand side of formula 'fixed' is needed"
  )

  expect_error(
    align_me(data = input_data, scaling = "one_parameter"),
    "Left and right-hand side of formula 'scaling' is needed"
  )

  expect_error(
    align_me(data = input_data, error = "one_parameter"),
    "Left and right-hand side of formula 'error' is needed"
  )

})
