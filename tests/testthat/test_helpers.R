
# get_symbols() -----------------------------------------------------------

test_that("get_symbols()", {
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


# split_for_scaling() -----------------------------------------------------

test_that("split_for_scaling()", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )
    sim_data_long <- read_wide(sim_data_wide_file, description = seq_len(3))

    effects_values <- list(
        biological_values = c("name", "time", "condition"),
        scaling_values = c("name", "ID"),
        error_values = "name"
    )

    effects_values_1 <- list(
        biological_values = c("name", "time", "condition"),
        scaling_values = c("name"),
        error_values = "name"
    )

    effects_values_2 <- list(
        biological_values = c("name", "time", "condition"),
        scaling_values = c("ID"),
        error_values = "name"
    )

    expect_equal(
        length(
            split_for_scaling(
                data = sim_data_long,
                effects_values,
                parameter_fit_scale = "linear",
                normalize_input = TRUE
            )
        ),
        14
    )

    expect_equal(
        length(
            split_for_scaling(
                data = sim_data_long,
                effects_values_1,
                parameter_fit_scale = "linear",
                normalize_input = TRUE
            )
        ),
        nrow(unique(sim_data_long["name"]))
    )

    expect_equal(
        length(
            split_for_scaling(
                data = sim_data_long,
                effects_values_2,
                parameter_fit_scale = "linear",
                normalize_input = TRUE
            )
        ),
        1
    )

    expect_warning(
        split_for_scaling(
            data = sim_data_long,
            effects_values_2,
            parameter_fit_scale = "log2",
            normalize_input = TRUE
        ),
        paste0(
            "'normalize_input == TRUE' is only competable with ",
            "'parameter_fit_scale == linear'. 'normalize_input' was ignored."
        )
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
            biological = yi ~ name + time + condition,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "linear"
        ),
        "All input checks passed."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            biological = yi ~ name + time + condition + wrong,
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "linear"
        ),
        paste0(
            "Not all column names set in 'biological', 'scaling' and 'error' ",
            "are present in 'data'."
        )
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            # biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "linear"
        ),
        "All of model, error_model, biological, scaling, error must be set."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            # model = "yi / sj",
            error_model = "value * sigmaR",
            # biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            # error = sigmaR ~ name + 1,
            parameter_fit_scale = "linear"
        ),
        "All of model, error_model, biological, scaling, error must be set."
    )

    expect_error(
        input_check(
            data = sim_data_long,
            model = "yi / sj",
            error_model = "value * sigmaR",
            biological = "yi ~ name + time + condition",
            scaling = sj ~ name + ID,
            error = sigmaR ~ name + 1,
            parameter_fit_scale = "wrong"
        ),
        "'parameter_fit_scale' must be 'linear', 'log', 'log2' or 'log10'."
    )
})


# generate_initial_pars() -------------------------------------------------
test_that("generate_initial_pars()", {
    parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    levels_list <- list(
        biological = c(
            "biological_1",
            "biological_2",
            "biological_3",
            "biological_4",
            "biological_5",
            "biological_6",
            "biological_7"
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
            parameter_fit_scale = "linear",
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
            parameter_fit_scale = "log",
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
            parameter_fit_scale = "log2",
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
            parameter_fit_scale = "log2",
            levels_list = list(
                biological = c(
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
            parameter_fit_scale = "linear",
            levels_list = list(
                biological = c(
                    "biological_a",
                    "biological_b",
                    "biological_c"
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


# objective_function() ----------------------------------------------------
test_that("objective_function()", {
    test_initial_parameters <- c(
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        sj = 1,
        sj = 1,
        sigmaR = 1
    )

    test_data_fit <- data.frame(
        name = c(
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT"
        ),
        time = c(
            0, 5, 10, 20, 30, 60, 240, 0
        ),
        value = c(
            1.0210929,
            1.2130989,
            0.8659901,
            0.9323670,
            1.0052727,
            0.9780511,
            1.1585450,
            0.8255822
        ),
        sigma = c(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        ),
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo",
            "pAKT_0_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_2"
        ),
        error = c(
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1"
        )
    )
    row.names(test_data_fit) <- c(
        "pAKT1",
        "pAKT2",
        "pAKT3",
        "pAKT4",
        "pAKT5",
        "pAKT6",
        "pAKT7",
        "pAKT8"
    )




    test_parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    test_levels_list <- list(
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test", "pAKT_2"
        ),
        error = c(
            "pAKT_1"
        )
    )

    test_effects_pars <- list(
        biological_pars = "yi",
        scaling_pars = "sj",
        error_pars = "sigmaR"
    )

    test_c_strenth <- 1000

    test_all_levels <- c(
        "pAKT_0_0Uml Epo",
        "pAKT_5_0Uml Epo",
        "pAKT_10_0Uml Epo",
        "pAKT_20_0Uml Epo",
        "pAKT_30_0Uml Epo",
        "pAKT_60_0Uml Epo",
        "pAKT_240_0Uml Epo",
        "pAKT_test",
        "pAKT_2",
        "pAKT_1"
    )

    test_mask <- generate_mask(
        test_initial_parameters,
        test_parameters,
        test_all_levels,
        test_data_fit
    )

    test_model_expr <- parse(text = "yi[biological]/sj[scaling]")
    test_error_model <- parse(
        text = "(yi[biological]/sj[scaling])* sigmaR[error]"
    )

    test_model_jacobian_expr <- list(
        biological = parse(text = "1/sj[scaling]"),
        scaling = parse(text = "-(yi[biological]/sj[scaling]^2)"),
        error = parse(text = "0")
    )

    test_error_model_jacobian_expr <- list(
        biological = parse(text = "1/sj[scaling]*sigmaR[error]"),
        scaling = parse(
            text = "-(yi[biological]/sj[scaling]^2*sigmaR[error])"
        ),
        error = parse(text = "(yi[biological]/sj[scaling])")
    )

    test_constraint_expr <- parse(text = "1e3 * (mean( yi ) - 1)")





    test_pass_parameter_list <- list(
        parameters = test_parameters,
        parameter_fit_scale = "linear",
        c_strength = test_c_strenth,
        model_expr = test_model_expr,
        error_model_expr = test_error_model,
        model_jacobian_expr = test_model_jacobian_expr,
        error_model_jacobian_expr = test_error_model_jacobian_expr,
        constraint_expr = test_constraint_expr
    )

    test_pass_parameter_list2 <- list(
        data_fit = test_data_fit,
        levels_list = test_levels_list,
        mask = test_mask,
        effects_pars = test_effects_pars
    )


    # * actual tests ----------------------------------------------------------


    expect_equal(
        objective_function(
            current_parameters = test_initial_parameters,
            pass_parameter_list = test_pass_parameter_list,
            pass_parameter_list2 = test_pass_parameter_list2,
            calculate_derivative = TRUE
        )$value,
        0.12445657
    )


    expect_equal(
        round(
            objective_function(
                current_parameters = test_initial_parameters,
                pass_parameter_list = test_pass_parameter_list,
                pass_parameter_list2 = test_pass_parameter_list2,
                calculate_derivative = TRUE
            )$gradient,
            digits = 8
        ),
        c(
            4.244916840,
            1.482979920,
            2.232102490,
            2.126117550,
            1.989399000,
            2.042934290,
            1.632636970,
            -13.463094600,
            -2.287992460,
            15.751086860
        )
    )
})


# split_for_scaling -------------------------------------------------------

test_that("split_for_scaling()", {
    sim_data_wide_file <- system.file(
        "extdata", "sim_data_wide.csv",
        package = "blotIt3"
    )
    sim_data_long <- read_wide(sim_data_wide_file, description = seq_len(3))

    test_effect_values <- list(
        biological_values = c("name", "time", "condition"),
        scaling_values = c("name", "ID"),
        error_values = c("name")
    )

    expected_values <- c(
        116.83827,
        138.80850,
        99.09068,
        106.68584,
        115.02805,
        111.91323,
        132.56618,
        94.46702
    )

    expect_warning(
        split_for_scaling(
            data = sim_data_long,
            effects_values = test_effect_values,
            normalize_input = TRUE,
            parameter_fit_scale = "log"
        ),
        paste0(
            "'normalize_input == TRUE' is only competable with ",
            "'parameter_fit_scale == linear'. 'normalize_input' was ignored."
        )
    )


    expect_equal(
        split_for_scaling(
            data = sim_data_long,
            effects_values = test_effect_values,
            normalize_input = FALSE,
            parameter_fit_scale = "linear"
        )[[1]]$value,
        expected_values
    )

    expect_equal(
        split_for_scaling(
            data = sim_data_long,
            effects_values = test_effect_values,
            normalize_input = TRUE,
            parameter_fit_scale = "linear"
        )[[1]]$value,
        expected_values / mean(expected_values)
    )

    expect_equal(
        split_for_scaling(
            data = sim_data_long,
            effects_values = test_effect_values,
            normalize_input = FALSE,
            parameter_fit_scale = "log"
        )[[1]]$value,
        expected_values
    )
})


# resolve_function() -----------------------------------------------------
test_that("resolve_function()", {
    test_initial_parameters <- c(
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        yi = 1,
        sj = 1,
        sj = 1,
        sigmaR = 1
    )

    test_data_fit <- data.frame(
        name = c(
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT",
            "pAKT"
        ),
        time = c(
            0, 5, 10, 20, 30, 60, 240, 0
        ),
        value = c(
            1.0210929,
            1.2130989,
            0.8659901,
            0.9323670,
            1.0052727,
            0.9780511,
            1.1585450,
            0.8255822
        ),
        sigma = c(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            NaN
        ),
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo",
            "pAKT_0_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_test",
            "pAKT_2"
        ),
        error = c(
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1",
            "pAKT_1"
        )
    )
    row.names(test_data_fit) <- c(
        "pAKT1",
        "pAKT2",
        "pAKT3",
        "pAKT4",
        "pAKT5",
        "pAKT6",
        "pAKT7",
        "pAKT8"
    )




    test_parameters <- c(
        "biological" = "yi",
        "scaling" = "sj",
        "error" = "sigmaR"
    )

    test_levels_list <- list(
        biological = c(
            "pAKT_0_0Uml Epo",
            "pAKT_5_0Uml Epo",
            "pAKT_10_0Uml Epo",
            "pAKT_20_0Uml Epo",
            "pAKT_30_0Uml Epo",
            "pAKT_60_0Uml Epo",
            "pAKT_240_0Uml Epo"
        ),
        scaling = c(
            "pAKT_test", "pAKT_2"
        ),
        error = c(
            "pAKT_1"
        )
    )

    test_effects_pars <- list(
        biological_pars = "yi",
        scaling_pars = "sj",
        error_pars = "sigmaR"
    )

    test_c_strenth <- 1000

    test_all_levels <- c(
        "pAKT_0_0Uml Epo",
        "pAKT_5_0Uml Epo",
        "pAKT_10_0Uml Epo",
        "pAKT_20_0Uml Epo",
        "pAKT_30_0Uml Epo",
        "pAKT_60_0Uml Epo",
        "pAKT_240_0Uml Epo",
        "pAKT_test",
        "pAKT_2",
        "pAKT_1"
    )

    test_mask <- generate_mask(
        test_initial_parameters,
        test_parameters,
        test_all_levels,
        test_data_fit
    )

    test_model_expr <- parse(text = "yi[biological]/sj[scaling]")
    test_error_model <- parse(
        text = "(yi[biological]/sj[scaling])* sigmaR[error]"
    )

    test_model_jacobian_expr <- list(
        biological = parse(text = "1/sj[scaling]"),
        scaling = parse(text = "-(yi[biological]/sj[scaling]^2)"),
        error = parse(text = "0")
    )

    test_error_model_jacobian_expr <- list(
        biological = parse(text = "1/sj[scaling]*sigmaR[error]"),
        scaling = parse(
            text = "-(yi[biological]/sj[scaling]^2*sigmaR[error])"
        ),
        error = parse(text = "(yi[biological]/sj[scaling])")
    )

    test_constraint_expr <- parse(text = "1e3 * (mean( yi ) - 1)")





    test_pass_parameter_list <- list(
        parameters = test_parameters,
        parameter_fit_scale = "linear",
        c_strength = test_c_strenth,
        model_expr = test_model_expr,
        error_model_expr = test_error_model,
        model_jacobian_expr = test_model_jacobian_expr,
        error_model_jacobian_expr = test_error_model_jacobian_expr,
        constraint_expr = test_constraint_expr
    )

    test_pass_parameter_list2 <- list(
        data_fit = test_data_fit,
        levels_list = test_levels_list,
        mask = test_mask,
        effects_pars = test_effects_pars
    )



    # * actual tests ----------------------------------------------------------


    expect_equal(
        resolve_function(
            current_parameters = test_initial_parameters,
            pass_parameter_list = test_pass_parameter_list,
            pass_parameter_list2 = test_pass_parameter_list2,
            calculate_derivative = FALSE
        )$residuals,
        c(
            -0.0210929,
            -0.2130989,
            0.1340099,
            0.0676330,
            -0.0052727,
            0.0219489,
            -0.1585450,
            0.1744178,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            1.0000000,
            0.0000000
        )
    )
})
