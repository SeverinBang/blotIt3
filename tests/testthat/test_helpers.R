
# get_symbols() -----------------------------------------------------------

test_that("read_symbols()", {
    expect_equal(
        get_symbols("left1 ~ right11 + right12"),
        c("left1", "right11", "right12")
    )
    expect_equal(
        get_symbols(
            "left2 ~ (right21 * right22 / (right23 * right24)) - right25"
        ),
        c("left2", "right21", "right22", "right23", "right24", "right25")
    )
})





# replace_symbols() -------------------------------------------------------

test_that("replace_symbols()", {
    expect_equal(
        replace_symbols(
            what = "before",
            by = "after",
            x = "this ~ before"
        ),
        "this~after"
    )

    expect_equal(
        replace_symbols(
            what = c("before1", "before2", "before3"),
            by = c("after1", "after2", "after3"),
            x = "this ~ before1*before2/before3"
        ),
        "this~after1*after2/after3"
    )
})



# analyze_blocks() --------------------------------------------------------

test_that("analyze_blocks()", {
    test_block_matrix <- matrix(
        c(
            1, 0, 0, 0,
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 1, 1, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
            0, 0, 0, 1
        ),
        nrow = 7,
        byrow = TRUE
    )
    expect_equal(
        analyze_blocks(
            test_block_matrix
        ),
        list(
            c(1, 2),
            c(3, 4, 5),
            c(6, 7)
        )
    )
    expect_equal(
        length(
            analyze_blocks(
                test_block_matrix
            )
        ),
        3
    )
})


# input_check() -----------------------------------------------------------

test_that("input_check()", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )
    sim_data_long <- read_wide(sim_data_wide_file, description = seq_len(3))

    expect_equal(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            distinguish = yi ~ name + time + condition,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            input_scale = "linear"
        ),
        "All input checks passed."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            distinguish = yi ~ name + time + condition + wrong,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            input_scale = "linear"
        ),
        paste0(
            "Not all column names set in 'distinguish', 'scaling' and 'error' ",
            "are present in 'data'."
        )
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            # distinguish = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            input_scale = "linear"
        ),
        "All of model, error_model, distinguish, scaling, error must be set."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            # model = "yi / sj",
            error_model = "value * sigmaR",
            # distinguish = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            # error = sigmaR ~ name + 1,
            input_scale = "linear"
        ),
        "All of model, error_model, distinguish, scaling, error must be set."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            distinguish = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            input_scale = "wrong"
        ),
        "'input_scale' must be 'linear', 'log', 'log2' or 'log10'."
    )
})


# generate_initial_pars() -------------------------------------------------
test_that("generate_initial_pars()", {
    parameters <- c(
        "distinguish" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    levels_list <- list(
        distinguish = c(
            "distinguish_1",
            "distinguish_2",
            "distinguish_3",
            "distinguish_4",
            "distinguish_5",
            "distinguish_6",
            "distinguish_7"
        ),
        scaling = c(
            "scaling_1",
            "scaling_2"
        ),
        error = "error1"
    )


    expect_equal(
        generate_initial_pars(
            parameters = parameters,
            input_scale = "linear",
            levels_list
        ),
        c(
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "sj" = 1,
            "sj" = 1,
            "sigmaR" = 1
        )
    )

    expect_equal(
        generate_initial_pars(
            parameters = parameters,
            input_scale = "log",
            levels_list
        ),
        c(
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "sj" = 0,
            "sj" = 0,
            "sigmaR" = 0
        )
    )

    expect_equal(
        generate_initial_pars(
            parameters = parameters,
            input_scale = "log2",
            levels_list
        ),
        c(
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "sj" = 0,
            "sj" = 0,
            "sigmaR" = 0
        )
    )

    expect_equal(
        generate_initial_pars(
            parameters = parameters,
            input_scale = "log2",
            levels_list = list(
                distinguish = c(
                "pAKT_0_0Uml Epo",
                "pAKT_5_0Uml Epo",
                "pAKT_10_0Uml Epo"
                ),
                scaling = c(
                    "scaling_1",
                    "scaling_2"
                ),
                error = "error1"
            )
        ),
        c(
            "yi" = 0,
            "yi" = 0,
            "yi" = 0,
            "sj" = 0,
            "sj" = 0,
            "sigmaR" = 0
        )
    )

    expect_equal(
        generate_initial_pars(
            parameters = parameters,
            input_scale = "linear",
            levels_list = list(
                distinguish = c(
                    "distinguish_a",
                    "distinguish_b",
                    "distinguish_c"
                ),
                scaling = c(
                    "scaling_a",
                    "scaling_b",
                    "scaling_c"
                ),
                error = c(
                    "error_a",
                    "error_b",
                    "error_c"
                )
            )
        ),
        c(
            "yi" = 1,
            "yi" = 1,
            "yi" = 1,
            "sj" = 1,
            "sj" = 1,
            "sj" = 1,
            "sigmaR" = 1,
            "sigmaR" = 1,
            "sigmaR" = 1
        )
    )
})
