
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

    expect_warning(
        read_wide(
            sim_data_wide_file,
            description = c("time", "does_not_exist")
        ),
        paste0(
            "Not all columns proposed by argument 'description' are available",
            " in file.\nTaking the available ones."
        )
    )
})


