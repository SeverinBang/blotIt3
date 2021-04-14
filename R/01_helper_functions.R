
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

scale_target <- function(current_data,
                         effects_values,
                         average_techn_rep,
                         input_scale) {
    # developement helper only
    if (FALSE) {
        current_data <- to_be_scaled[[1]]
    }

    # Add column for error
    current_data$sigma <- NaN

    # Add a dummy column filled with 1
    current_data[["1"]] <- "1"

    # Initialize the parameter for the current target
    current_parameter <- NULL

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

    initial_parameters <- generate_initial_pars(
        parameters,
        input_scale,
        levels_list
    )

    mask <- generate_mask(
        initial_parameters,
        parameters,
        all_levels,
        data_fit
    )

    fit_pars_distinguish <- NULL
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
#' Residuals of model evaluations for a set of \code{test_parameters} and
#' \code{fit_pars_distinguish} are calculated by evaluation of \link{rss_model}.
#'
#' @param test_parameters named vector of vectors to be tested currently
#' @param  fit_pars_distinguish named vector that will contain the fitted values
#' of the distinguish parameters.
#'
#' @return large list
#'
#' @noRd
residual_function <- function(test_parameters,
                              fit_pars_distinguish,
                              parameters,
                              levels_list,
                              effects_pars,
                              deriv = TRUE) {
    if (FALSE) {
        test_parameters <- initial_parameters
        fit_pars_distinguish <- NULL
    }
    pars_all <- c(test_parameters, fit_pars_distinguish)
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
            res = rss_model(
                parlist, # generated in "res_fn"
                data_fit,
                model_expr,
                errmodel_expr,
                constraint_expr,
                jac_model_expr,
                jac_errmodel_expr,
                deriv = deriv
            )
        ),
        parlist
    )
}
