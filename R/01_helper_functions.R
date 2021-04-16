
# get_symbols() -----------------------------------------------------------

#' Function to retrieve symbols from formula-formatted input
#'
#' Formula like strings are parsed, and returned in form of a list
#'
#' @param char input
#' @param exclude parts of the formula which will be excluded from output
#'
#' @return data frame with columns "name", "time", "value" and other
#' columns describing the measurements.
#'
#' @noRd
get_symbols <- function(char, exclude = NULL) {
    # Parse input data
    char <- char[char != "0"]
    out <- parse(text = char, keep.source = TRUE)
    out <- utils::getParseData(out)

    # Split the input string at non-caracter symbols
    names <- unique(out$text[out$token == "SYMBOL"])

    # Remove strings passed as "exclude"
    if (!is.null(exclude)) {
        names <- names[!names %in% exclude]
    }
    return(names)
}




# identify_effects() ------------------------------------------------------

#' Method to identify the names and values of passed todistinguish, scaling and
#' error parameters
#'
#' The formula-formatted input will be parsed by \link{get_symbols} and the
#' respective output then passed back
#'
#' @param distinguish input for distinguish parameters
#' @param scaling input for scaling parameters
#' @param error input for error parameters
#'
#' @return two lists: values and parameter names
#'
#' @noRd

identify_effects <- function(distinguish = NULL, scaling = NULL, error = NULL) {

    # Get distinguish scale and error parameters
    distinguish_values <- get_symbols(as.character(distinguish)[3])
    scaling_values <- get_symbols(as.character(scaling)[3])
    error_values <- get_symbols(as.character(error)[3])

    # Add intercepts
    if (attr(terms(distinguish), "intercept") != 0 &
        length(distinguish_values) == 2) {
        distinguish_values <- c(distinguish_values, "1")
    }
    if (attr(terms(scaling), "intercept") != 0 &
             length(scaling_values) == 1) {
        scaling_values <- c(scaling_values, "1")
    }
    if (attr(terms(error), "intercept") != 0 &
        length(error_values) == 1) {
        error_values <- c(error_values, "1")
    }

    # Determine to which class parameters belong
    distinguish_pars <- get_symbols(as.character(distinguish)[2])
    scaling_pars <- get_symbols(as.character(scaling)[2])
    error_pars <- get_symbols(as.character(error)[2])

    effects_values <- list(
        distinguish_values = distinguish_values,
        scaling_values = scaling_values,
        error_values = error_values
    )
    effects_pars <- list(
        distinguish_pars = distinguish_pars,
        scaling_pars = scaling_pars,
        error_pars = error_pars
    )

    return(list(effects_values = effects_values, effects_pars = effects_pars))
}



# replace_symbol ----------------------------------------------------------

#' Method for replacing parts of a string. Used to paste expresisons for the
#' (error-) model.
#'
#' The formula-formatted input will be parsed by \link{get_symbols} and the
#' respective output then passed back
#'
#' @param what string or vector of strings to be replaced
#' @param by string or vector of strings that will take the replaced place
#' @param x string in which \code{what} will be replaced by \code{by}
#'
#' @return string, \code{x} with the replaced object.
#'
#' @noRd
#'
replace_symbols <- function(what, by, x) {
    x_orig <- x
    is_not_zero <- which(x != "0")
    x <- x[is_not_zero]
    my_names <- names(x)
    x_parsed <- parse(text = x, keep.source = TRUE)
    data <- utils::getParseData(x_parsed)
    by <- rep(by, length.out = length(what))
    names(by) <- what
    data$text[data$text %in% what] <- by[data$text[data$text %in% what]]
    data <- data[data$token != "expr", ]
    breaks <- c(0, which(diff(data$line1) == 1), length(data$line1))
    out <- lapply(
        seq_len((length(breaks) - 1)),
        function(i) {
            paste(
                data$text[seq((breaks[i] + 1), (breaks[i + 1]))],
                collapse = ""
            )
        }
    )
    names(out) <- my_names
    out <- unlist(out)
    x_orig[is_not_zero] <- out
    return(x_orig)
}



# paste_() ----------------------------------------------------------------

#' Just a shorthand
#' @param ... the thing that should be pasted.
#'
#' @noRd
paste_ <- function(...) paste(..., sep = "_")




# analyze_blocks() --------------------------------------------------------

#' Method to analyze which elements of the given matrix with <number of unique
#' set of scaling parameters> columns, and <unique set of distinguish
#' parameters> rows, are on the same scale (same block) and pass a list of lists
#' of the respective indices.
#' Each row describes a unique set of distinguishable conditions, e.g. a
#' specific target measured under a specific condition at one time point.
#' The columns describe the sets of scaling parameters as target name and
#' gel (in case of western blot).
#'
#'
#' @param block_matrix Input matrix, with entries equal to 1 wherever the set of
#' distinguishable effects (row) is measured und the respective set of scaling
#' effects (column), i.e. each row has a one for each scaling under which it was
#' measured and each column has a one indicating which distinguishable effects
#' where measured at the respective scaling.
#' All entries are either one or zero.
#'
#' @return List with one entry per scale. Each entry contains a list of row
#' indices of the sets of distinguishable effects measured under the respective
#' scale.
#'
#' @noRd

analyze_blocks <- function(block_matrix) {
    out <- which(apply(block_matrix, 1, sum) == 0)
    if (length(out) > 0) {
        block_matrix <- block_matrix[-out, ]
        cat("matrix contains zero rows which have been eliminated\n")
    }

    number_unique_distinguish <- dim(block_matrix)[1]
    r_components <- list()
    c_components <- list()

    counter <- 0
    while (length(unlist(r_components)) < number_unique_distinguish) {
        counter <- counter + 1

        if (length(unlist(r_components)) == 0) {
            w <- 1
        } else {
            my_sample <- (1:number_unique_distinguish)[-unlist(r_components)]
            w <- min(my_sample)
        }

        repeat {
            v <- unique(
                rapply(
                    as.list(w),
                    function(i) which(block_matrix[i, ] == 1)
                )
            )
            w_new <- unique(
                rapply(
                    as.list(v),
                    function(j) which(block_matrix[, j] == 1)
                )
            )
            if (length(w_new) == length(w)) break
            w <- w_new
        }
        r_components[[counter]] <- w
        c_components[[counter]] <- v
    }

    return(r_components)
}



# input_check() -----------------------------------------------------------

#' Check input parameters for structural errors
#'
#'
#' @noRd

input_check <- function(data = NULL,
                        model = NULL,
                        error_model = NULL,
                        distinguish = NULL,
                        scaling = NULL,
                        error = NULL,
                        input_scale = NULL) {

    # Check if parameters are present
    if (is.null(model) | is.null(error_model) | is.null(distinguish) |
        is.null(scaling) | is.null(error)) {
        stop(
            "All of model, error_model, distinguish, scaling, error ",
            "must be set."
        )
    }

    # check data
    column_names <- names(data)
    expected_names <- union(
        get_symbols(as.character(distinguish)[3]),
        union(
            get_symbols(as.character(scaling)[3]),
            get_symbols(as.character(error)[3])
        )
    )

    if (!all(expected_names %in% column_names)) {
        stop(
            "Not all column names set in 'distinguish', 'scaling' and ",
            "'error' are present in 'data'."
        )
    }

    # Check scaling
    if (!(as.character(input_scale) %in% c("linear", "log", "log2", "log10"))) {
        stop(
            "'input_scale' must be 'linear', 'log', 'log2' or 'log10'."
        )
    }



    # Stop if formulas have the wrong specification
    if (length(as.character(distinguish)) == 1 |
        length(as.character(scaling)) == 1 | length(as.character(error)) == 1) {
        stop("Do not pass distinguish, scaling or error as string.")
    }
    if (length(as.character(distinguish)) < 3) {
        stop("Left and right-hand side of formula 'distinguish' is needed")
    }
    if (length(as.character(scaling)) < 3) {
        stop("Left and right-hand side of formula 'scaling' is needed")
    }
    if (length(as.character(error)) < 3) {
        stop("Left and right-hand side of formula 'error' is needed")
    }

    return_msg <- "All input checks passed."
    return(return_msg)
}



# scale_target() ----------------------------------------------------------

#' Scaling function applied to the current dataset
#'
#' Heart of the package. The passed set is assumed to belong to one scale. The
#' respective scaling factors are calculated by optimizing
#' \link{objective_function}.
#'
#' @param current_data current data set, belongs to one scale.
#' @param effects_values list of names of the columns defined as distinguish,
#' scaling and error effects, respectively
#' @param average_techn_rep logical, indicates if the technical replicates
#' @param input_scale character, defining in which scale the input is passed.
#' must be one of c('linear', 'log', 'log2', 'log10')
#' should be averaged
#'
#' @return A new data set including the original, scaled and
#' predicted data as well as the estimated parameters
#'
#' @noRd

scale_target <- function(
    current_data,
    pass_parameter_list
) {
    # developement helper only
    if (FALSE) {
        current_data <- to_be_scaled[[1]]
    }

    # Retrieve parameters from list
    effects_values <- pass_parameter_list$effects_values
    parameter_data <- pass_parameter_list$parameter_data
    average_techn_rep <- pass_parameter_list$average_techn_rep
    verbose <- pass_parameter_list$verbose
    covariates <- pass_parameter_list$covariates
    parameters <- pass_parameter_list$parameters
    input_scale <- pass_parameter_list$input_scale
    effects_pars <- pass_parameter_list$effects_pars
    model_expr <- pass_parameter_list$model_expr
    error_model_expr <- pass_parameter_list$error_model_expr
    constraint_expr <- pass_parameter_list$constraint_expr
    model_jacobian_expr <- pass_parameter_list$model_jacobian_expr
    error_model_jacobian_expr <- pass_parameter_list$error_model_jacobian_expr
    c_strength <- pass_parameter_list$c_strength
    normalize <- pass_parameter_list$normalize


    current_name <- unique(current_data$name)

    # Add column for error
    current_data$sigma <- NaN

    # Add a dummy column filled with 1
    current_data[["1"]] <- "1"

    # Initialize the parameter for the current target
    current_parameter <- NULL

    if (!is.null(parameter_data)) {
        current_parameter <- parameter_data[
            parameter_data$name %in% current_data$name,
            ]
    }


    # Build a list of string of distinguish/scaling values present
    data_fit_distinguish <- do.call(
        paste_, current_data[, effects_values[[1]], drop = FALSE]
    )

    data_fit_scaling <- do.call(
        paste_, current_data[, effects_values[[2]], drop = FALSE]
    )

    # Reducing datapoints
    if (average_techn_rep) {
        if (verbose) {
            cat("Analyzing technical replicates ... ")
        }
        groups <- interaction(data_fit_distinguish, data_fit_scaling)
        if (any(duplicated(groups))) {
            current_data <- do.call(rbind, lapply(
                unique(groups),
                function(g) {
                    subdata <- current_data[groups == g, ]
                    outdata <- subdata[1, ]
                    outdata$value <- mean(subdata$value)
                    return(outdata)
                }
            ))
            cat(
                "data points that could not be distinguished by either ",
                "distinguish\nor scaling variables have been averaged.\n"
            )
        } else {
            if (verbose) {
                cat("none found.\n")
            }
        }
    }

    # compile dataset fit for fitting
    data_fit <- data.frame(
        current_data[
            ,
            union(c("name", "time", "value", "sigma"), covariates)
        ],
        distinguish = do.call(
            paste_,
            current_data[, effects_values[[1]], drop = FALSE]
        ),
        scaling = do.call(
            paste_,
            current_data[, effects_values[[2]], drop = FALSE]
        ),
        error = do.call(
            paste_,
            current_data[, effects_values[[3]], drop = FALSE]
        ),
        stringsAsFactors = FALSE
    )

    # Retrieve (unique) list of distinguish, scaling and error values
    levels_list <- list(
        distinguish = unique(as.character(data_fit$distinguish)),
        scaling = unique(as.character(data_fit$scaling)),
        error = unique(as.character(data_fit$error))
    )

    # Combine above lists all_levels will therefore have the length of the
    # sum of all different distinguish, scaling and error values (in that order)
    all_levels <- unlist(
        lapply(
            seq_along(parameters),
            function(k) {
                switch(names(parameters)[k],
                    distinguish = levels_list[[1]],
                    scaling = levels_list[[2]],
                    error = levels_list[[3]]
                )
            }
        )
    )
#
    initial_parameters <- generate_initial_pars(
        parameters,
        input_scale,
        levels_list
    )
#
    mask <- generate_mask(
        initial_parameters,
        parameters,
        all_levels,
        data_fit
    )

    fit_pars_distinguish <- NULL

    if (verbose) {
        cat("Starting fit\n")
    }

    pass_parameter_list2  <-  list(
        data_fit = data_fit,
        levels_list = levels_list,
        fit_pars_distinguish = fit_pars_distinguish,
        parameters = parameters,
        effects_pars = effects_pars,
        c_strength = c_strength,
        mask = mask,
        initial_parameters = initial_parameters
    )

# * call of trust() function ----------------------------------------------


    fit_result <- trust::trust(
        objfun =  objective_function,
        parinit = initial_parameters,
        rinit = 1,
        rmax = 10,
        blather = verbose,
        # # folowing: additional parameters for objective_function()
        # fit_pars_distinguish = fit_pars_distinguish,
        # calculate_derivative = TRUE,
        # data_fit = data_fit,
        # parameters = parameters,
        # levels_list = levels_list,
        # effects_pars = effects_pars,
        # mask = mask,
        # c_strength = c_strength,
        # model_expr = model_expr,
        # error_model_expr = error_model_expr,
        # constraint_expr = constraint_expr,
        # model_jacobian_expr = model_jacobian_expr,
        # error_model_jacobian_expr = error_model_jacobian_expr,
        # input_scale = input_scale
        # # parameterlists
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2
    )

    if (!fit_result$converged) {
        warning(paste("Non-converged fit for target", current_name))
    } else {
        cat("Fit converged.\n")
    }

    residuals_fit <- residual_function(
        current_parameters =  fit_result$argument,
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2,
        calculate_derivative = FALSE
        )

    bessel <- sqrt(
        nrow(data_fit) / (nrow(data_fit) - length(initial_parameters) + normalize)
        )

    # Get singular values (roots of non negative eigenvalues of M^* \cdot M)
    # here it is synonymous with eigenvalue.
    single_values <- svd(fit_result[["hessian"]])[["d"]]

    # Define a tollerance threshold, the root of the machine precission is
    # a usual value for this threshold.
    tol <- sqrt(.Machine$double.eps)

    # Define nonidentifiable as being "to small to handle", judged by the
    # above defined threshold
    non_identifiable <- which(single_values < tol * single_values[1])
    if (length(non_identifiable) > 0) {
        warning("Eigenvalue(s) of Hessian below tolerance. Parameter
                  uncertainties might be underestimated.")
    }

    # Generate parameter table
    parameter_table <- data.frame(
        name = current_name,
        level = c(
            rep(levels_list[[1]], length(effects_pars[[1]])),
            rep(levels_list[[2]], length(effects_pars[[2]])),
            rep(levels_list[[3]], length(effects_pars[[3]]))
        ),
        parameter = names(fit_result$argument),
        value = fit_result$argument,

        # Calculating the error from the inverse of the Fisher information
        # matrix which is in this case the Hessian, to which the above
        # calculated Bessel correction is applied.
        sigma = as.numeric(
            sqrt(diag(2 * MASS::ginv(fit_result$hessian)))
        ) * bessel,
        nll = fit_result$value,
        no_pars = length(fit_result$argument) - normalize,
        no_data = nrow(data_fit)
    )

    if (input_scale == "log") {
        parameter_table$value <- exp(parameter_table$value)
        parameter_table$sigma <- parameter_table$value * parameter_table$sigma
    } else if (input_scale == "log2") {
        parameter_table$value <- 2^(parameter_table$value)
        parameter_table$sigma <- parameter_table$value * parameter_table$sigma
    } else if (input_scale == "log10") {
        parameter_table$value <- 10^(parameter_table$value)
        parameter_table$sigma <- parameter_table$value * parameter_table$sigma
    }

    if (verbose) {
        cat("Estimated parameters on non-log scale:\n")
        print(parameter_table)
        cat(
            "converged:", fit_result$converged, ", iterations:",
            fit_result$iterations, "\n"
        )
        cat("-2*LL: ", fit_result$value, "on", nrow(data_fit) +
                normalize - length(fit_result$argument), "degrees of freedom\n")
    }

    attr(parameter_table, "value") <- fit_result$value
    attr(parameter_table, "df") <- nrow(data_fit) + normalize -
        length(data_fit$argument)

    # Predicted data
    out_predicted <- data_fit
    out_predicted$sigma <- fit_result$sigma * bessel
    out_predicted$value <- fit_result$prediction

    # Initialize list for scaled values
    initial_values_scaled <- rep(0, nrow(data_fit))
    if (verbose) {
        cat("Inverting model ... ")
    }


    # rootSolve::multiroot(
    #     f = evaluate_model,
    #     start = initial_values_scaled,
    #     jacfunc = evaluate_model_jacobian,
    #     par_list = residuals_fit[-1],
    #     verbose = TRUE,
    #     pass_parameter_list = pass_parameter_list,
    #     pass_parameter_list2 = pass_parameter_list2
    # )
    # print("DEEEBUG")

    values_scaled <- try(
        # my_multiroot
        rootSolve::multiroot(
            f = evaluate_model,
            start = initial_values_scaled,
            jacfunc = evaluate_model_jacobian,
            par_list = residuals_fit[-1],
            verbose = FALSE,
            pass_parameter_list = pass_parameter_list,
            pass_parameter_list2 = pass_parameter_list2
        )$root,
        silent = TRUE
    )
    if (verbose) {
        cat("done\n")
    }

    if (inherits(values_scaled, "try-error")) {
        out_scaled <- NULL
        warning("Rescaling to common scale not possible.
                  Equations not invertible.")
    } else {
        my_deriv <- abs(
            evaluate_model_jacobian(
                values_scaled,
                residuals_fit[-1],
                pass_parameter_list = pass_parameter_list,
                pass_parameter_list2 = pass_parameter_list2
            )
        )
        sigmas_scaled <- fit_result$sigma * bessel / my_deriv
        if (input_scale == "log") {
            values_scaled <- exp(values_scaled)
            sigmas_scaled <- values_scaled * sigmas_scaled
        } else if (input_scale == "log2") {
            values_scaled <- 2^(values_scaled)
            sigmas_scaled <- values_scaled * sigmas_scaled
        } else if (input_scale == "log10") {
            values_scaled <- 10^(values_scaled)
            sigmas_scaled <- values_scaled * sigmas_scaled
        }
        out_scaled <- current_data
        out_scaled$value <- values_scaled
        out_scaled$sigma <- sigmas_scaled
    }

    # Aligned
    no_initial <- length(levels_list[[1]])

    # Use one datapoint per unique set of fixed parameters
    out_aligned <- current_data[
        !duplicated(data_fit$distinguish),
        intersect(effects_values[[1]], colnames(current_data))
    ]

    # The values are the fitted parameters for the respective fixed
    # parameter ensembles.
    out_aligned$value <- residuals_fit[[effects_pars[[1]][1]]]
    if (input_scale == "log") {
        out_aligned$value <- exp(out_aligned$value)
    } else if (input_scale == "log2") {
        out_aligned$value <- 2^(out_aligned$value)
    } else if (input_scale == "log10") {
        out_aligned$value <- 10^(out_aligned$value)
    }

    # The sigmas are calculated by evaluating the fisher information matrix
    out_aligned$sigma <- as.numeric(
        sqrt(diag(2 * MASS::ginv(fit_result$hessian)))
    )[seq_len(no_initial)] * bessel
    if (input_scale != "linear") {
        out_aligned$sigma <- out_aligned$value * out_aligned$sigma
    }


    # Get the original data
    out_original_with_parameters <- current_data
    for (k in seq_along(parameters)) {
        effect <- names(parameters)[k]
        my_levels <- as.character(data_fit[[effect]])
        index0 <- which(as.character(parameter_table$parameter) == parameters[k])
        index1 <- match(my_levels, as.character(parameter_table$level[index0]))
        index <- index0[index1]
        out_original_with_parameters[[parameters[k]]] <- parameter_table$value[index]
    }
    out <- list(
        out_predicted,
        out_scaled,
        out_aligned,
        current_data,
        out_original_with_parameters,
        parameter_table
    )

    return(out)
}



# generate_initial_pars() -------------------------------------------------

#' Method to generate a set of initial parameters for \link{scale_target}
#'
#' @noRd

generate_initial_pars <- function(parameters,
                                  input_scale,
                                  levels_list) {
    # Create a named vector with initial values for all parameters. The name
    # is the parameter, and the initial value is 0 for log and 1 for linear.
    # each parameter name-value entry is repeated as many times as there are
    # measurements of the i'th target with this parameter.
    # Example: distinguish is passed as "ys ~ condition", so the entry with name
    #   "ys" will be repeated for as many times as there are measurements
    #   with the same "condition" entry.
    initial_parameters <- do.call(
        "c",
        lapply(
            seq_along(parameters),
            # Go through all parameters with current parameter n
            function(n) {
                if (input_scale != "linear") {
                    v <- 0
                } else {
                    v <- 1
                }
                # Let l be the number of measurements of same type (distinguish,
                # scaling or error) for the current n
                l <- switch(
                    # Get "distinguish", "scaling" or "error" from the current n
                    names(parameters)[n],
                    distinguish = length(levels_list[[1]]),
                    scaling = length(levels_list[[2]]),
                    error = length(levels_list[[3]])
                )

                # Get the n'th parameter
                p <- as.character(parameters[n])

                # The variable "parameter_values" will have l times the entry
                # "v" (0 for log, 1 for linear) with the corresponding parameter
                # as the name.
                parameter_values <-
                    structure(
                        rep(v, l),
                        names = rep(p, l)
                    )
                return(parameter_values)
            }
        )
    )
    return(initial_parameters)
}


# generate_mask() ---------------------------------------------------------

#' Creates a identifier list
#'
#' A list with one element per entry in initial_parameters is generated. For
#' each of those elements, it is checked which elements of the data set
#' \code{data_fit} coincides with the current element of \code{all_levels}.
#' The entry is a list with 1 indicating where the current element of
#' \code{all_levels} is found in the respective effect-column of
#' \code{data_fit}.
#' Keep in mind, that the \code{initial_parameters} and \code{all_levels} are
#' structured analogue.
#'
#' @param initial_parameters named vector, set of parameters as output of
#' \link{generate_initial_pars}
#' @param parameters named vector, with the 'values' of the three effects
#' @param all_levels character vector, strings describing all levels in the
#' current data set
#' @param data_fit data frame containing the data that will be fitted
#'
#' @noRd

generate_mask <- function(initial_parameters,
                          parameters,
                          all_levels,
                          data_fit) {
    mask <- lapply(
        seq_along(initial_parameters),
        function(k) {
            effect <- names(
                parameters[
                    match(
                        names(initial_parameters[k]),
                        parameters
                    )
                ]
            )
            mask_vector <- as.numeric(
                data_fit[[effect]] == all_levels[k]
            )
            return(mask_vector)
        }
    )
}



# residual_function() -----------------------------------------------------

#' Calculate residuals for optimization
#'
#' Residuals of model evaluations for a set of \code{current_parameters} and
#' \code{fit_pars_distinguish} are calculated by evaluation of \link{rss_model}.
#'
#' @param current_parameters named vector of vectors to be tested currently
#' @param  fit_pars_distinguish named vector that will contain the fitted values
#' of the distinguish parameters.
#'
#' @return large list
#'
#' @noRd
residual_function <- function(current_parameters,
                              pass_parameter_list,
                              pass_parameter_list2,
                              calculate_derivative) {
    if (FALSE) {
        current_parameters <- initial_parameters
        fit_pars_distinguish <- NULL
    }

    fit_pars_distinguish <- pass_parameter_list2$fit_pars_distinguish
    parameters <- pass_parameter_list$parameters
    levels_list <- pass_parameter_list2$levels_list
    effects_pars <- pass_parameter_list2$effects_pars
    data_fit <- pass_parameter_list2$data_fit
    model_expr <- pass_parameter_list$model_expr
    error_model_expr <- pass_parameter_list$error_model_expr
    constraint_expr <- pass_parameter_list$constraint_expr
    model_jacobian_expr <- pass_parameter_list$model_jacobian_expr
    error_model_jacobian_expr <- pass_parameter_list$error_model_jacobian_expr







    pars_all <- c(current_parameters, fit_pars_distinguish)
    par_list <- lapply(
        parameters,
        function(n) {
            subpar <- pars_all[names(pars_all) == n]
            if (n %in% effects_pars[1]) {
                # Rename entries of the first list (fixed)
                names(subpar) <- levels_list[[1]]
            }
            if (n %in% effects_pars[2]) {
                # Rename entries of the second list (latent)
                names(subpar) <- levels_list[[2]]
            }
            if (n %in% effects_pars[3]) {
                # Rename entries of the third list (error)
                names(subpar) <- levels_list[[3]]
            }
            return(subpar)
        }
    )

    # Rename the thre lists
    names(par_list) <- parameters

    # Create list of residuals and attach the corresponding current
    # parameters
    c(
        list(
            residuals = rss_model(
                par_list,
                data_fit,
                model_expr,
                error_model_expr,
                constraint_expr,
                model_jacobian_expr,
                error_model_jacobian_expr,
                calculate_derivative = calculate_derivative
            )
        ),
        par_list
    )
}



# rss_model() -------------------------------------------------------------

#' Calculate model residual sum of squares
#'
#' @return list of residuals
#'
#' @noRd
rss_model <- function(
    par_list,
    data_fit,
    model_expr,
    error_model_expr,
    constraint_expr,
    model_jacobian_expr,
    error_model_jacobian_expr,
    calculate_derivative = TRUE) {
    # First argument: data, second argument: expressions which are evalued with
    # the data of arg 1

    if (FALSE) {
        name <- as.list(data_fit)$name
        time <- as.list(data_fit)$time
        value <- as.list(data_fit)$value
        sigma <- as.list(data_fit)$sigma
        fixed <- as.list(data_fit)$fixed
        latent <- as.list(data_fit)$latent
        error <- as.list(data_fit)$error
        yi <- par_list$ys
        sj <- par_list$sj
        sigmaR <- par_list$sigmaR
        calculate_derivative = TRUE
    }
    with(
        # Gives list with entries: name, time, value, sigma, fixed, latent, error,
        # ys, sj, sigmaR, the last three are from "parlist".
        c(as.list(data_fit), par_list),
        {
            # Generates list with entry <NoOfMeasurements> x 1 as entries, save
            # this list as "var" and initialize the "prediction" list with it
            prediction <- var <- rep(1, length(value))

            # Get the prediction by evaluating the model
            prediction[seq_along(prediction)] <- eval(model_expr)

            # Create residuals: differences between prediction and measurements
            res <- prediction - value

            # Generate variances by squaring the evaluated error model
            var[seq_along(var)] <- eval(error_model_expr)^2

            # Evaluate constraints
            constr <- eval(constraint_expr)

            # Create list of residuals
            residuals <- c(res / sqrt(var), var, constr)

            # Initialize variables for derivatives
            residual_deriv <- variance_deriv <- NULL

            # Calculate derivatives if wished
            if (calculate_derivative) {
                jac <- lapply(
                    seq_len(length(par_list)),
                    function(k) {
                        # Initialize lists
                        v_mod <- v_err <- rep(0, length(value))

                        # Get derivatives for model and errormodel
                        v_mod[seq_len(length(value))] <- eval(model_jacobian_expr[[k]])
                        v_err[seq_len(length(value))] <- eval(error_model_jacobian_expr[[k]])

                        residual_deriv <- v_mod / sqrt(var) - v_err * res / var
                        variance_deriv <- v_err * 2 * sqrt(var)
                        list(residual_deriv, variance_deriv)
                    }
                )

                residual_deriv <- lapply(jac, function(j) j[[1]])
                variance_deriv <- lapply(jac, function(j) j[[2]])
            }
            attr(residuals, "residual_deriv") <- residual_deriv
            attr(residuals, "variance_deriv") <- variance_deriv

            return(residuals)
        }
    )
}


# objective_function() ----------------------------------------------------

#' Objective function for trust region optimizer
#'
objective_function <- function(current_parameters,
                               # fit_pars_distinguish,
                               # data_fit,
                               # parameters,
                               # levels_list,
                               # effects_pars,
                               # c_strength,
                               # mask,
                               # par_list = par_list,
                               # model_expr,
                               # error_model_expr,
                               # constraint_expr,
                               # model_jacobian_expr,
                               # error_model_jacobian_expr,
                               # input_scale,
                               pass_parameter_list,
                               pass_parameter_list2,
                               calculate_derivative = TRUE
                               ) {

    if (FALSE) {
        current_parameters <- initial_parameters
        calculate_derivative = TRUE
    }

    # Retrieve parameters from list
    effects_values <- pass_parameter_list$effects_values
    parameter_data <- pass_parameter_list$parameter_data
    average_techn_rep <- pass_parameter_list$average_techn_rep
    verbose <- pass_parameter_list$verbose
    covariates <- pass_parameter_list$covariates
    parameters <- pass_parameter_list$parameters
    input_scale <- pass_parameter_list$input_scale
    effects_pars <- pass_parameter_list$effects_pars
    model_expr <- pass_parameter_list$model_expr
    error_model_expr <- pass_parameter_list$error_model_expr
    constraint_expr <- pass_parameter_list$constraint_expr
    model_jacobian_expr <- pass_parameter_list$model_jacobian_expr
    error_model_jacobian_expr <- pass_parameter_list$error_model_jacobian_expr
    c_strength <- pass_parameter_list$c_strength
    normalize <- pass_parameter_list$normalize




    data_fit <- pass_parameter_list2$data_fit
    levels_list <- pass_parameter_list2$levels_list
    fit_pars_distinguish <- pass_parameter_list2$fit_pars_distinguish
    parameters <- pass_parameter_list2$parameters
    effects_pars <- pass_parameter_list2$effects_pars
    c_strength <- pass_parameter_list2$c_strength
    mask <- pass_parameter_list2$mask

    no_data <- nrow(data_fit)

    # Recover residuals from output of res_fn()
    calculated_residuals <- residual_function(
        current_parameters = current_parameters,
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2,
        calculate_derivative = calculate_derivative
        )$residuals

    # Retrieve residuals of model, errormodel and constraint as well as
    # the derivatives.
    residuals <- calculated_residuals[1:no_data]
    variances <- calculated_residuals[(no_data + 1):(2 * no_data)]
    constraint <- calculated_residuals[2 * no_data + 1]
    residual_deriv <- attr(calculated_residuals, "residual_deriv")
    variance_deriv <- attr(calculated_residuals, "variance_deriv")

    # Set bessel correction factor to 1
    bessel <- 1

    # Calculate the current value
    value <- sum(residuals^2) + bessel * sum(log(variances)) + constraint^2
    gradient <- NULL
    hessian <- NULL

    # Get derivatives for the measurements by applying the predefined mask
    if (calculate_derivative) {
        calculated_residuals_jacobian <- do.call(cbind, lapply(
            seq_along(current_parameters),
            function(k) {

                # Get the "class" (fixed, latent, or error) of the current k'th
                # parameter
                which_par <- match(names(current_parameters)[k], parameters)

                # Apply the mask to the residuals
                residual_jacobian <- residual_deriv[[which_par]] * mask[[k]]
                variance_jacobian <- variance_deriv[[which_par]] * mask[[k]]
                constrain_jacobian <- as.numeric(names(current_parameters[k]) ==
                                             parameters["distinguish"]) * c_strength / length(levels_list[[1]])

                # Convert to log if wanted
                if (input_scale == "log") {
                    constrain_jacobian <- constrain_jacobian * exp(current_parameters[k])
                } else if (input_scale == "log2") {
                    constrain_jacobian <- constrain_jacobian * 2^(current_parameters[k])
                } else if (input_scale == "log10") {
                    constrain_jacobian <- constrain_jacobian * 10^(current_parameters[k])
                }

                # Stitch the results together
                c(residual_jacobian, variance_jacobian, constrain_jacobian)
            }
        ))

        # Split the above list into corresponding parts
        residual_jacobian <- calculated_residuals_jacobian[1:no_data, , drop = FALSE]
        jac_vars <- calculated_residuals_jacobian[
            (no_data + 1):(2 * no_data), ,
            drop = FALSE
        ]
        constrain_jacobian <- calculated_residuals_jacobian[2 * no_data + 1, , drop = FALSE]

        # Compose to gradient vector and hessian matrix
        gradient <- as.vector(2 * residuals %*% residual_jacobian +
                                  (bessel / variances) %*% jac_vars + 2 * constraint * constrain_jacobian)
        hessian <- 2 * t(rbind(residual_jacobian, constrain_jacobian)) %*%
            (rbind(residual_jacobian, constrain_jacobian))
    }

    # Takes residual (ultimative res <- prediction - value and res / sqrt(var)
    # from rss_model()) to retrive the presdiction
    prediction <- residuals * sqrt(variances) + data_fit$value

    # Calculate errors
    sigma <- sqrt(variances)

    # Compose all of the above to one list
    list(
        value = value, gradient = gradient, hessian = hessian,
        prediction = prediction, sigma = sigma
    )
}


# evaluate_model() --------------------------------------------------------

#' Model evaluation method
#'
#' @noRd
#'

evaluate_model <- function(initial_parameters,
                           par_list,
                           pass_parameter_list,
                           pass_parameter_list2
                           # = calculated_residuals[-1],
                           # effects_pars = effects_pars,
                           # model_expr = model_expr,
                           # data_fit = data_fit
                           ) {
    # Create a list with with values from 1 to the number of parameters

    model_expr <- pass_parameter_list$model_expr
    data_fit <- pass_parameter_list2$data_fit
    # initial_parameters <- pass_parameter_list2$initial_parameters
    effects_pars <- pass_parameter_list$effects_pars



    distinguish <- seq_along(initial_parameters)

    # Generate a list with the entries:
    #   parameters, the sequence saved generated as "fixed", name, time,
    #   value, sigma, fixed, latent, error, ys, sj, sigmaR
    my_list <- c(
        list(initial_parameters, distinguish = distinguish),
        as.list(data_fit),
        par_list
    )
    names(my_list)[1] <- effects_pars[[1]][1]
    # to evlaluate manually
    # yi <- my_list$yi
    # fixed <- my_list$fixed
    # name <- my_list$name
    # time <- my_list$time
    # value <- my_list$value
    # sigma <- my_list$sigma
    # distinguish <- my_list$distinguish
    # scaling <- my_list$scaling
    # error <- my_list$error
    # yi <- my_list$yi
    # sj <- my_list$sj
    # sigmaR <- my_list$sigmaR
    values <- with(my_list, eval(model_expr) - data_fit$value)

    return(values)
}



# evaluate_model_jacobian() -----------------------------------------------

#' Model jacobian evaluation method
#'
#' @noRd
#'

evaluate_model_jacobian <- function(initial_parameters,
                                    par_list,
                                    pass_parameter_list,
                                    pass_parameter_list2
                                    # ,
                                    # effects_pars = effects_pars,
                                    # model_expr = model_expr,
                                    # data_fit = data_fit
                                    ) {

    data_fit <- pass_parameter_list2$data_fit
    # initial_parameters <- pass_parameter_list2$initial_parameters
    model_derivertive_expr <- pass_parameter_list$model_derivertive_expr
    effects_pars <- pass_parameter_list$effects_pars


    distinguish <- seq_along(initial_parameters)
    my_list <- c(
        list(initial_parameters, distinguish = distinguish), as.list(data_fit),
        par_list
    )
    names(my_list)[1] <- effects_pars[[1]][1]
    derivative_values <- with(my_list, eval(model_derivertive_expr))

    return(derivative_values)
}
