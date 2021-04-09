# align_me() - model parsing ----------------------------------------------

test_that("align_me() - model parsing", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )

    input_data <- read_wide(sim_data_wide_file, description = seq_len(3))

    expect_error(
        align_me(data = input_data, distinguish = "one_parameter"),
        "Do not pass distinguish, scaling or error as string."
    )

    expect_error(
        align_me(data = input_data, distinguish = left ~ right, scaling = NULL),
        "Left and right-hand side of formula 'scaling' is needed"
    )

    expect_error(
        align_me(data = input_data, distinguish = NULL, error = left ~ right2),
        "Left and right-hand side of formula 'distinguish' is needed"
    )
})
