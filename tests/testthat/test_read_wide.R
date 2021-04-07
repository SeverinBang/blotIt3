
# read_wide() -------------------------------------------------------------


test_that("read_wide works properly", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )
    sim_data_wide <- read.csv(sim_data_wide_file)

    expect_equal(
        length(
            unique(read_wide(sim_data_wide_file, description = seq_len(3))$name)
        ),
        length(sim_data_wide) - 3
    )
    expect_error(
      read_wide(sim_data_wide_file),
      "Specify columns containing descriptions."
    )
})
