# align_me() - model parsing ----------------------------------------------

test_that("align_me() - model parsing", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )

    input_data <- read_wide(sim_data_wide_file, description = seq_len(3))

    expect_error(
        align_me(
            data = input_data,
            model = "yi / sj",
            error_model = "value * sigmaR",
            biological = "one_parameter",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1
        ),
        "Do not pass biological, scaling or error as string."
    )

    expect_error(
        align_me(
            data = input_data,
            model = "yi / sj",
            error_model <- "value * sigmaR",
            biological <- yi ~ name + time + condition,
            error = sigmaR ~ name + 1
        ),
        "All of model, error_model, biological, scaling, error must be set"
    )

    expect_error(
        align_me(
            data = input_data,
            biological = NULL,
            error = left ~ right2
        ),
        "All of model, error_model, biological, scaling, error must be set."
    )

    expect_error(
        align_me(
            data = input_data,
            model = "yi / sj",
            error_model = "value * sigmaR",
            biological = yi ~ condition,
            scaling = wrong ~ ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "linear"
        ),
        "Not all paramters are defined in either arguments
         'scaling', 'biological' or 'error'"
    )

    expect_error(
        align_me(
            data = input_data,
            model = "yi / sj",
            error_model = "value * sigmaR",
            biological = yi ~ condition,
            scaling = sj ~ ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "lineara"
        ),
        "'parameter_fit_scale' must be 'linear', 'log', 'log2' or 'log10'."
    )
})
