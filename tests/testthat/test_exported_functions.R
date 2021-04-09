
# split_for_scaling() -----------------------------------------------------

test_that("split_for_scaling()", {
  sim_data_wide_file <- system.file(
    "extdata", "sim_data_wide.csv",
    package = "blotIt3"
  )
  sim_data_long <- read_wide(sim_data_wide_file, description = seq_len(3))

  expect_equal(
    length(
      split_for_scaling(
        data = sim_data_long,
        distinguish_values = c("name", "time", "condition"),
        scaling_values = c("name", "replicate")
        )
      ),
    nrow(unique(sim_data_long[c("name", "replicate")]))
  )

})
