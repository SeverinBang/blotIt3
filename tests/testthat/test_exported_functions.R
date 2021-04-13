
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
                scaling_values = c("name", "ID"),
                input_scale = "linear",
                normalize_input = TRUE
            )
        ),
        14
    )

    expect_equal(
        length(
            split_for_scaling(
                data = sim_data_long,
                distinguish_values = c("name", "time", "condition"),
                scaling_values = "name",
                input_scale = "linear",
                normalize_input = TRUE
            )
        ),
        nrow(unique(sim_data_long["name"]))
    )

    expect_equal(
        length(
            split_for_scaling(
                data = sim_data_long,
                distinguish_values = c("name", "time", "condition"),
                scaling_values = "ID",
                input_scale = "linear",
                normalize_input = TRUE
            )
        ),
        1
    )

    expect_warning(
        split_for_scaling(
            data = sim_data_long,
            distinguish_values = c("name", "time", "condition"),
            scaling_values = c("name", "ID"),
            input_scale = "log2",
            normalize_input = TRUE
        ),
        paste0(
            "'normalize_input == TRUE' is only competable with ",
            "'input_scale == linear'. 'normalize_input' was ignored."
        )
    )
})
