
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
#' @importFrom utils getParseData
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
#' @importFrom utils getParseData
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



# split_for_scaling()  ----------------------------------------------------


#' split_data
#'
#' Split data in independent blocks according to distinguish and scaling
#' variables as being defined for \link{align_me}. Each block will be given an
#' individual scaling factor.
#'
#' @param data data frame with columns "name", "time", "value" and others
#' @param effects_values two-sided formula, see \link{align_me}
#' @param scaling two-sided formula, see \link{align_me}
#' @param normalize_input logical, if set to TRUE, the input data will be
#' normalized by dividing all entries belonging to one scaling factor by their
#'  respective mean. This prevents convergence failure on some hardware when the
#'  data for different scaling effects differ by to many orders of magnitude.
#' @return list of data frames
#'
#' @noRd
split_for_scaling <- function(data,
                              effects_values,
                              normalize_input,
                              parameter_fit_scale) {
    if (!"1" %in% colnames(data)) {
        data["1"] <- 1
        intercept <- FALSE
    } else {
        intercept <- TRUE
    }

    # Construnct strings containing the values of the respective effects
    scaling_strings <- Reduce(paste_, data[effects_values[[2]]])
    distinguish_strings <- Reduce(paste_, data[effects_values[[1]]])

    scaling_strings_unique <- unique(scaling_strings)
    distinguish_strings_unique <- unique(distinguish_strings)

    # Initialize matrix
    block_matrix <- matrix(
        0,
        ncol = length(scaling_strings_unique),
        nrow = length(distinguish_strings_unique)
    )

    # For every datapoint set the entry in the matrix to 1, corresponding to the
    # index in the respective unique lists of fixed (row) and specific (col)
    for (i in seq_len(nrow(data))) {
        myrow <- which(distinguish_strings_unique == distinguish_strings[i])
        mycol <- which(scaling_strings_unique == scaling_strings[i])
        block_matrix[myrow, mycol] <- 1
    }

    # compile list of scalings
    list_of_scalings <- analyze_blocks(block_matrix)

    # Remove the "1"-column added above
    if (!intercept) {
        data <- data[, -which(colnames(data) == "1")]
    }

    # Compile the list
    list_out <- lapply(
        list_of_scalings,
        function(l) {
            data[distinguish_strings %in% distinguish_strings_unique[l], ]
        }
    )

    # normalize the data
    if (normalize_input) {
        if (parameter_fit_scale == "linear") {
            for (i in seq_len(length(list_out))) {
                list_out[[i]]$value <- list_out[[i]]$value /
                    mean(list_out[[i]]$value)
            }
        } else {
            warning(
                "'normalize_input == TRUE' is only competable with ",
                "'parameter_fit_scale == linear'. 'normalize_input' was ignored."
            )
        }
    }

    return(list_out)
}



# input_check() -----------------------------------------------------------

#' Check input parameters for structural errors
#'
#' the input of \link{align_me} is checked for user mistakes.
#'
#' @param data Input data, usually output of \link{read_wide}
#' @param model Model definition as a string
#' @param error_model Error model definition as a string
#' @param distinguish Definition of the distinguish effects
#' @param scaling Definition of the scaling effects
#' @param error Definition of the error effects
#' @param parameter_fit_scale String, either of c('linear', 'log', 'log2', 'log10').
#'
#' @noRd

input_check <- function(data = NULL,
                        model = NULL,
                        error_model = NULL,
                        distinguish = NULL,
                        scaling = NULL,
                        error = NULL,
                        parameter_fit_scale = NULL,
                        output_scale = NULL) {

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
    if (!(as.character(parameter_fit_scale) %in%
          c("linear", "log", "log2", "log10")) |
        !(as.character(output_scale) %in%
          c("linear", "log", "log2", "log10"))) {
        stop(
            "'parameter_fit_scale' and 'output_scale' must be 'linear', ",
            "'log' 'log2' or 'log10'."
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
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effects_values}}{
#'      Named list of vectors. The names correspond to the effects and the
#'      vectors contain the 'values'; the right hand side of the respective
#'      effects parameters passed to \link{align_me} the names of the columns
#'      associated with the respective effect.
#'  }
#'  \item{\code{parameter_data}}{
#'      If the \code{data} parameter passed to \link{align_me} is already an
#'      output of \link{align_me}, it contains the fitted parameters, otherwise
#'      \code{NULL}.
#'  }
#'  \item{\code{average_techn_rep}}{
#'      Logical, if \code{TRUE} technical replicates, that can not be separated
#'      by the respective \code{distinguish} and \code{scaling} values will be
#'      averaged.
#'  }
#'  \item{\code{verbose}}{
#'      Logical, if \code{TRUE} additional output will be generated.
#'  }
#'  \item{\code{covariates}}{
#'      String, the covariates as defined by the (error-) model expressions as
#'      defined in the respective parameters of \link{align_me}
#'  }
#'  \item{\code{parameters}}{
#'      Named vector, contains the variables for the three effects, and the
#'      respective effects as names.
#'  }
#'  \item{\code{parameter_fit_scale}}{
#'      String, defining the scale of the data see \code{parameter_fit_scale} in
#'      \link{align_me}.
#'  }
#'  \item{\code{effects_pars}}{
#'      Similar to \code{parameters}, but as a named list.
#'  }
#'  \item{\code{model_expr}}{
#'      What is passed to \code{model} in \link{align_me} but with the entries
#'      of \code{parameters} replaced by the respective name (e.g. 'yi' is
#'      replaced by 'distinguish'). The result is parsed as an expression.
#'  }
#'  \item{\code{error_model_expr}}{
#'      Same as \code{model_expr} but with the value of \code{error_model} from
#'      \link{align_me}.
#'  }
#'  \item{\code{constraint_expr}}{
#'      If \code{normalize} is set to \code{TRUE} in \link{align_me},
#'      \code{constraint_expr} contains expression of the normalization
#'      constraint, otherwise it is empty.
#'  }
#'  \item{\code{model_jacobian_expr}}{
#'      Named list with the derivatives of \code{model_expr} with respect to the
#'      respective effects as entries.
#'  }
#'  \item{\code{error_model_jacobian_expr}}{
#'      Analogue to \code{model_jacobian_expr} but with \code{error_model_expr}.
#'  }
#'  \item{\code{c_strength}}{
#'      Numerical 1000, if \code{normalize} is set to \code{TRUE} in
#'      \link{align_me}, otherwise unused.
#'  }
#'  \item{\code{normalize}}{
#'      Logical, direct taken from the \code{normalize} parameter of
#'      \link{align_me}.
#'  }
#' }
#'
#' @return A list of \code{data.frame}S with the entries:
#' \describe{
#'  \item{\code{out_prediction}}{
#'      \code{data.frame} with the columns \code{name}, \code{time},
#'      \code{value}, \code{sigma}, \code{distinguish}, \code{scaling} and
#'      \code{error}.
#'
#'      \code{values} contains the predictions by evaluation of the (error-)
#'      model with the fitted parameters on the original scale (if
#'      \code{normalize_input} is set to \code{FALSE} in \link{align_me}). The
#'      last three columns contain strings pasted from the respective values of
#'      the current entry.
#'  }
#'  \item{\code{out_scaled}}{
#'      \code{data.frame} with the original columns and additionally
#'      \code{sigma} and \code{1}.
#'
#'      The \code{values} are the original ones scaled to common scale by
#'      applying the inverse of the model with the fitted parameters. The
#'      \code{sigma} entries are the results of the evaluated error model scaled
#'      to common scale adhering to Gaussian error propagation.
#'  }
#'  \item{\code{out_aligned}}{
#'      \code{data.frame} with the original column names, the column
#'      \code{value} contains the estimated distinguish parameters.
#'
#'      The \code{values} are estimated distinguish parameters, while the errors
#'      \code{sigma} are from the Fisher information.
#'  }
#'  \item{\code{current_data}}{
#'      \code{data.frame} with the original data with added columns \code{sigma}
#'      and \code{1}.
#'  }
#'  \item{\code{out_orig_w_parameters}}{
#'      Same \code{data.frame} as \code{current_data} but with added columns for
#'      the \code{parameters} ad the respective fitted values.
#'  }
#'  \item{\code{parameter_table}}{
#'      \code{data.frame} with the columns:
#'      \itemize{
#'          \item{\code{level}:}{
#'          String pasting all the unique effect entries.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{parameter}:}{
#'          Parameter associated with the current effect (compare with
#'          \code{parameters}).
#'          }
#'      }
#'      \itemize{
#'          \item{\code{value}:}{
#'          Value of the estimated parameters. The ones corresponding to the
#'          distinguished parameters coincide with the \code{values} column of
#'          \code{out_aligned}.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{sigma}:}{
#'          Errors of the estemated parameter values by utilizing the Fisher
#'          information
#'          }
#'      }
#'      \itemize{
#'          \item{\code{nll}:}{
#'          Negative of twice the log likelihood of the fit in which the
#'          paramters are estimated.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{no_pars}:}{
#'          Number of fitted parameters (distinguish, scaling and error), minus
#'          one if \code{normalize = TRUE} in the call of \link{align_me}.
#'          }
#'      }
#'      \itemize{
#'          \item{\code{no_data}:}{
#'          Length of the current data file, i.e. number of measurements in the
#'          set that is currently scaled.
#'          }
#'      }
#'  }
#' }
#'
#' @importFrom rootSolve multiroot
#' @importFrom trust trust
#' @importFrom MASS ginv
#'
#' @noRd

scale_target <- function(current_data,
                         pass_parameter_list) {
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
    parameter_fit_scale <- pass_parameter_list$parameter_fit_scale
    effects_pars <- pass_parameter_list$effects_pars
    normalize <- pass_parameter_list$normalize
    output_scale <-  pass_parameter_list$output_scale
    iterlim <- pass_parameter_list$iterlim


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
        parameter_fit_scale,
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

    pass_parameter_list2 <- list(
        data_fit = data_fit,
        levels_list = levels_list,
        fit_pars_distinguish = fit_pars_distinguish,
        effects_pars = effects_pars,
        mask = mask,
        initial_parameters = initial_parameters
    )

    # * call of trust() function ----------------------------------------------


    fit_result <- trust::trust(
        objfun = objective_function,
        parinit = initial_parameters,
        rinit = 1,
        rmax = 10,
        iterlim = iterlim,
        blather = verbose,
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2
    )

    if (!fit_result$converged) {
        warning(paste("Non-converged fit for target", current_name))
    } else {
        cat("Fit converged.\n")
    }

    residuals_fit <- resolve_function(
        current_parameters = fit_result$argument,
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2,
        calculate_derivative = FALSE
    )

    bessel <- sqrt(
        nrow(data_fit) / (nrow(data_fit) - length(initial_parameters) +
            normalize)
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

    # * parameter table -------------------------------------------------------
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

    parameter_table$upper <- parameter_table$value + parameter_table$sigma
    parameter_table$lower <- parameter_table$value - parameter_table$sigma

    scaling_list_parameter_table <- scale_values(
        parameter_fit_scale = parameter_fit_scale,
        output_scale = output_scale,
        value = parameter_table$value,
        upper = parameter_table$upper,
        lower = parameter_table$lower,
        sigma = parameter_table$sigma

    )

    parameter_table$value <- scaling_list_parameter_table$value
    parameter_table$upper <- scaling_list_parameter_table$upper
    parameter_table$lower <- scaling_list_parameter_table$lower
    parameter_table$sigma <- scaling_list_parameter_table$sigma
#
#     if (parameter_fit_scale == "log") {
#         parameter_table$value <- exp(parameter_table$value)
#         parameter_table$upper <- exp(parameter_table$upper)
#         parameter_table$lower <- exp(parameter_table$lower)
#         parameter_table$sigma <- parameter_table$value * parameter_table$sigma
#     } else if (parameter_fit_scale == "log2") {
#         parameter_table$value <- 2^(parameter_table$value)
#         parameter_table$upper <- 2^(parameter_table$upper)
#         parameter_table$lower <- 2^(parameter_table$lower)
#         parameter_table$sigma <- parameter_table$value * parameter_table$sigma
#     } else if (parameter_fit_scale == "log10") {
#         parameter_table$value <- 10^(parameter_table$value)
#         parameter_table$upper <- 10^(parameter_table$upper)
#         parameter_table$lower <- 10^(parameter_table$lower)
#         parameter_table$sigma <- parameter_table$value * parameter_table$sigma
#     }

    if (verbose) {
        cat("Estimated parameters on non-log scale:\n")
        print(parameter_table)
        cat(
            "converged:", fit_result$converged, ", iterations:",
            fit_result$iterations, "\n"
        )
        cat("-2*LL: ", fit_result$value, "on", nrow(current_data) +
            normalize - length(fit_result$argument), "degrees of freedom\n")
    }

    attr(parameter_table, "value") <- fit_result$value
    attr(parameter_table, "df") <- nrow(current_data) + normalize -
        length(fit_result$argument)


    # * out_prediction --------------------------------------------------------
    # Predicted data
    out_prediction <- current_data
    out_prediction$sigma <- fit_result$sigma * bessel
    out_prediction$value <- fit_result$prediction


    # * out_scaled ------------------------------------------------------------
        # Initialize list for scaled values
    initial_values_scaled <- rep(0, nrow(current_data))
    if (verbose) {
        cat("Inverting model ... ")
    }


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

        upper_scaled <- values_scaled + sigmas_scaled
        lower_scaled <- values_scaled - sigmas_scaled

        scaling_list_scaled <- scale_values(
            parameter_fit_scale = parameter_fit_scale,
            output_scale = output_scale,
            value = values_scaled,
            upper = upper_scaled,
            lower = lower_scaled,
            sigma = sigmas_scaled
        )

        out_scaled <- current_data
        out_scaled$value <- scaling_list_scaled$value
        out_scaled$sigma <- scaling_list_scaled$sigma
        out_scaled$upper <- scaling_list_scaled$upper
        out_scaled$lower <- scaling_list_scaled$lower

        # values_scaled <- scaling_list_scaled$value
        #
        # if (parameter_fit_scale == "log") {
        #     values_scaled <- exp(values_scaled)
        #     upper_scaled <- exp(upper_scaled)
        #     lower_scaled <- exp(lower_scaled)
        #     sigmas_scaled <- values_scaled * sigmas_scaled
        #     upper_scaled <- values_scaled * upper_scaled
        #     lower_scaled <- values_scaled * lower_scaled
        # } else if (parameter_fit_scale == "log2") {
        #     values_scaled <- 2^(values_scaled)
        #     sigmas_scaled <- values_scaled * sigmas_scaled
        #     upper_scaled <- values_scaled * upper_scaled
        #     lower_scaled <- values_scaled * lower_scaled
        # } else if (parameter_fit_scale == "log10") {
        #     values_scaled <- 10^(values_scaled)
        #     sigmas_scaled <- values_scaled * sigmas_scaled
        #     upper_scaled <- values_scaled * upper_scaled
        #     lower_scaled <- values_scaled * lower_scaled
        # }
#
#         out_scaled <- current_data
#         out_scaled$value <- values_scaled
#         out_scaled$sigma <- sigmas_scaled
#         out_scaled$upper <- upper_scaled
#         out_scaled$lower <- lower_scaled

    }


    # * out_aligned -----------------------------------------------------------
    no_initial <- length(levels_list[[1]])

    # Use one datapoint per unique set of fixed parameters
    out_aligned <- current_data[
        !duplicated(data_fit$distinguish),
        intersect(effects_values[[1]], colnames(current_data))
    ]

    # The values are the fitted parameters for the respective fixed
    # parameter ensembles.
    out_aligned$value <- residuals_fit[[effects_pars[[1]][1]]]

    # The sigmas are calculated by evaluating the fisher information matrix
    out_aligned$sigma <- as.numeric(
        sqrt(diag(2 * MASS::ginv(fit_result$hessian)))
    )[seq_len(no_initial)] * bessel

    out_aligned$upper <- out_aligned$value + out_aligned$sigma
    out_aligned$lower <- out_aligned$value - out_aligned$sigma

    scaling_list_aligned <- scale_values(
        parameter_fit_scale = parameter_fit_scale,
        output_scale = output_scale,
        value = out_aligned$value,
        upper = out_aligned$upper,
        lower = out_aligned$lower,
        sigma = out_aligned$sigma
    )

    out_aligned$value <- scaling_list_aligned$value
    out_aligned$upper <- scaling_list_aligned$upper
    out_aligned$lower <- scaling_list_aligned$lower
    out_aligned$sigma <- scaling_list_aligned$sigma
    #
    #
    #
    #
    # if (parameter_fit_scale == "log") {
    #     out_aligned$value <- exp(out_aligned$value)
    #     out_aligned$upper <- exp(out_aligned$upper)
    #     out_aligned$lower <- exp(out_aligned$lower)
    #     out_aligned$sigma <- out_aligned$value * out_aligned$sigma
    # } else if (parameter_fit_scale == "log2") {
    #     out_aligned$value <- 2^(out_aligned$value)
    #     out_aligned$upper <- 2^(out_aligned$upper)
    #     out_aligned$lower <- 2^(out_aligned$lower)
    #     out_aligned$sigma <- out_aligned$value * out_aligned$sigma
    # } else if (parameter_fit_scale == "log10") {
    #     out_aligned$value <- 10^(out_aligned$value)
    #     out_aligned$upper <- 10^(out_aligned$upper)
    #     out_aligned$lower <- 10^(out_aligned$lower)
    #     out_aligned$sigma <- out_aligned$value * out_aligned$sigma
    # }
#
#
#     if (parameter_fit_scale != "linear") {
#         out_aligned$sigma <- out_aligned$value * out_aligned$sigma
#         out_aligned$upper <- out_aligned$value * out_aligned$upper
#         out_aligned$lower <- out_aligned$value * out_aligned$lower
#     }




    # * out_orig_w_parameters -------------------------------------------------
    out_orig_w_parameters <- current_data
    for (k in seq_along(parameters)) {
        effect <- names(parameters)[k]
        my_levels <- as.character(data_fit[[effect]])
        index0 <- which(
            as.character(parameter_table$parameter) == parameters[k]
        )
        index1 <- match(my_levels, as.character(parameter_table$level[index0]))
        index <- index0[index1]
        out_orig_w_parameters[[parameters[k]]] <- parameter_table$value[index]
    }
    out <- list(
        out_prediction = out_prediction,
        out_scaled = out_scaled,
        out_aligned = out_aligned,
        original = current_data,
        out_orig_w_parameters = out_orig_w_parameters,
        parameter_table = parameter_table
    )




    return(out)
}



# generate_initial_pars() -------------------------------------------------

#' Method to generate a set of initial parameters for \link{scale_target}
#'
#' @param parameters Named vector, contains the variables for the three effects,
#' and the respective effects as names.
#' @param parameter_fit_scale String: 'linear', 'log', 'log2' or 'log10'. Describing
#' what the initial value will be
#' @param levels_list Named list with one entry per effect containing a vector
#' of strings composed with the respective pasted entries of the effects columns
#'
#' @return Named vector with one entry per parameter, all are 1 if
#' \code{parameter_fit_scale = 'linear'} and 0 if not. The names are the entries of
#' \code{parameters} of the corresponding effect. Each entry is repeated as
#' often as there are parameters of the respective effect.
#'
#' @noRd

generate_initial_pars <- function(parameters,
                                  parameter_fit_scale,
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
                if (parameter_fit_scale != "linear") {
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
#' @return list with one entry per entry of \code{all_levels}. Each of this
#' entries is a numerical list of the length of \code{data_fit}. The entries are
#' either 1, if the corresponding entry of \code{data_fit} is resembled by the
#' current entry of \code{all_levels}.
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

    return(mask)
}



# resolve_function() -----------------------------------------------------

#' Function to sort the current set of parameters into the three effect
#' categories.
#'
#' Additionally, residuals of model evaluations for the set of
#' \code{current_parameters} \code{fit_pars_distinguish} are calculated by
#' evaluation of \link{rss_model}.
#'
#' @param current_parameters named vector of vectors to be tested currently
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{parameters}}{Named list of the parameters of the three effects,
#'  the names are the respective effects}
#' }
#' @param pass_parameter_list2 from the the \code{pass_parameter_list2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{fit_pars_distinguish}}{TODO: Check if it is always \code{NULL}?}
#'  \item{\code{levels_list}}{Named list of vectors. One entry per effect with
#'  the respective name. The entry contains a list of the unique strings
#'  composed from the entries of \code{effect_values} of the respective effect,
#'  i.e. the entries of the columns containing e.g. the distinguish-effects
#'  (name, time, condition etc.).}
#'  \item{\code{effects_pars}}{Named list of vectors. One entry per effect with
#'  the respective name. The entry then contains the string of the effect
#'  parameter, i.e. the variable name (e.g. "sj" for scaling).}
#' }
#'
#' @return named list of the residuals calculated in \link{rss_model} and the
#' \code{current_parameters} sorted in the effects with the respective entries
#' from \code{levels_list} as names.
#'
#' @noRd
#'
resolve_function <- function(current_parameters,
                             pass_parameter_list,
                             pass_parameter_list2,
                             calculate_derivative) {
    if (FALSE) {
        current_parameters <- initial_parameters
    }


    parameters <- pass_parameter_list$parameters
    fit_pars_distinguish <- pass_parameter_list2$fit_pars_distinguish
    levels_list <- pass_parameter_list2$levels_list
    effects_pars <- pass_parameter_list2$effects_pars

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
                pass_parameter_list,
                pass_parameter_list2,
                calculate_derivative = calculate_derivative
            )
        ),
        par_list
    )
}



# rss_model() -------------------------------------------------------------

#' Calculate model residuals
#'
#' Residuals are calculated by the difference (prediction - value), where
#' \code{prediction} is the evaluation of the model with the current set of
#' parameters and \code{value} the respective measurements. Optionally, the
#' derivatives are also calculated
#'
#'
#' @param par_list Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'distinguish')
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{model_expr}}{The argument of the \code{model} parameter of
#'  \link{align_me} parsed as an executable expression.}
#'  \item{\code{error_model_expr}}{The argument of the \code{error_model}
#'  parameter of \link{align_me} parsed as en executable expression.}
#'  \item{\code{constraint_expr}}{If the logical argument \code{normalize} of
#'  \link{align_me} is set to \code{TRUE}, an executable constraint expression
#'  is passed. If \code{normalize = FALSE} it is empty.}
#'  \item{\code{model_jacobian_expr}}{The Jacobian of the model as a named list
#'  with the respective derivatives as entries passed as expression.}
#'  \item{\code{error_model_jacobian_expr}}{The Jacobian of the error model as a
#'  named list with the respective derivatives as entries passed as expression.}
#' }
#'
#' @param pass_parameter_list2 from the the \code{pass_parameter_list2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{data_fit}}{The current dataset}
#' }
#'
#' @param calculate_derivative logical, if \code{TRUE}, derivatives will also be
#' calculated.
#'
#' @return list of residuals with the derivatives (optional) as attributes.
#'
#' @noRd
rss_model <- function(par_list,
                      pass_parameter_list,
                      pass_parameter_list2,
                      calculate_derivative = TRUE) {
    model_expr <- pass_parameter_list$model_expr
    error_model_expr <- pass_parameter_list$error_model_expr
    constraint_expr <- pass_parameter_list$constraint_expr
    model_jacobian_expr <- pass_parameter_list$model_jacobian_expr
    error_model_jacobian_expr <- pass_parameter_list$error_model_jacobian_expr

    data_fit <- pass_parameter_list2$data_fit

    with(
        # paste list with entries: name, time, value, sigma, the effects and
        # their paramters
        c(as.list(data_fit), par_list), {
            # Generate lists var and prediction with <NoOfMeasurements> entries
            # and initialize the it with 1
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
                        v_mod[seq_len(length(value))] <- eval(
                            model_jacobian_expr[[k]]
                        )
                        v_err[seq_len(length(value))] <- eval(
                            error_model_jacobian_expr[[k]]
                        )

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
#' @param current_parameters Current set of parameters to be tested. Will be
#' iteratively changed by the objective function.
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{parameters}}{
#'      Named vector, contains the variables for the three effects, and the
#'      respective effects as names.
#'  }
#'  \item{\code{parameter_fit_scale}}{
#'      String, defining the scale of the data see \code{parameter_fit_scale} in
#'      \link{align_me}.
#'  }
#'  \item{\code{c_strength}}{
#'      Numerical 1000, if \code{normalize} is set to \code{TRUE} in
#'      \link{align_me}, otherwise unused.
#'  }
#' }
#'
#' @param pass_parameter_list2 from the the \code{pass_parameter_list2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{data_fit}}{
#'      The current data set
#'  }
#'  \item{\code{levels_list}}{Named list of vectors. One entry per effect with
#'      the respective name. The entry contains a list of the unique strings
#'      composed from the entries of \code{effect_values} of the respective
#'      effect, i.e. the entries of the columns containing e.g. the distinguish
#'      effects (name, time, condition etc.).
#'  }
#'  \item{\code{mask}}{
#'      Output of \link{generate_mask}
#'  }
#' }
#'
#' @param calculate_derivative Logical, indicates if derivatives should be also
#' calculated.
#'
#' @return named list with the entries \code{value}, \code{gradient} and
#' \code{hessian}. The \code{value}S are the residual sum of squares with the
#' residuals as calculated by \link{rss_model} (via \link{resolve_function}).
#'
#' @noRd
objective_function <- function(current_parameters,
                               pass_parameter_list,
                               pass_parameter_list2,
                               calculate_derivative = TRUE) {
    if (FALSE) {
        current_parameters <- initial_parameters
        calculate_derivative <- TRUE
    }

    parameters <- pass_parameter_list$parameters
    parameter_fit_scale <- pass_parameter_list$parameter_fit_scale
    c_strength <- pass_parameter_list$c_strength

    data_fit <- pass_parameter_list2$data_fit
    levels_list <- pass_parameter_list2$levels_list
    mask <- pass_parameter_list2$mask

    no_data <- nrow(data_fit)

    # Recover residuals from output of res_fn()
    calculated_residuals <- resolve_function(
        current_parameters = current_parameters,
        pass_parameter_list = pass_parameter_list,
        pass_parameter_list2 = pass_parameter_list2,
        calculate_derivative = calculate_derivative
    )$residuals

    # Retrieve residuals of model, errormodel and constraint as well as
    # the derivatives.
    residuals <- calculated_residuals[1:no_data]
    variances <- calculated_residuals[seq((no_data + 1), (2 * no_data))]
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
                constrain_jacobian <- as.numeric(
                    names(current_parameters[k]) == parameters["distinguish"]
                ) * c_strength / length(levels_list[[1]])

                # Convert to log if wanted
                if (parameter_fit_scale == "log") {
                    constrain_jacobian <- constrain_jacobian *
                        exp(current_parameters[k])
                } else if (parameter_fit_scale == "log2") {
                    constrain_jacobian <- constrain_jacobian *
                        2^(current_parameters[k])
                } else if (parameter_fit_scale == "log10") {
                    constrain_jacobian <- constrain_jacobian *
                        10^(current_parameters[k])
                }

                # Stitch the results together
                c(residual_jacobian, variance_jacobian, constrain_jacobian)
            }
        ))

        # Split the above list into corresponding parts
        residual_jacobian <- calculated_residuals_jacobian[
            1:no_data, ,
            drop = FALSE
        ]
        jac_vars <- calculated_residuals_jacobian[
            (no_data + 1):(2 * no_data), ,
            drop = FALSE
        ]
        constrain_jacobian <- calculated_residuals_jacobian[
            2 * no_data + 1, ,
            drop = FALSE
        ]

        # Compose to gradient vector and hessian matrix
        gradient <- as.vector(
            2 * residuals %*% residual_jacobian + (bessel / variances) %*%
                jac_vars + 2 * constraint * constrain_jacobian
        )
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
#' Method to evaluate the model with a given set of parameters.
#'
#' @param initial_parameters Parameters used for model evaluation
#'
#' @param par_list Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'distinguish')
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effects_pars}}{
#'      Named list with one entry per effect. The values are the respective
#'      effect parameters.
#'  }
#'  \item{\code{model_expr}}{
#'      The argument of the \code{model} parameter of \link{align_me} parsed as
#'      an executable expression.
#'  }
#' }
#'
#' @param pass_parameter_list2 from the the \code{pass_parameter_list2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{data_fit}}{The current dataset.}
#' }
#'
#'
#' @noRd

evaluate_model <- function(initial_parameters,
                           par_list,
                           pass_parameter_list,
                           pass_parameter_list2) {
    # Create a list with with values from 1 to the number of parameters

    effects_pars <- pass_parameter_list$effects_pars
    model_expr <- pass_parameter_list$model_expr
    data_fit <- pass_parameter_list2$data_fit




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

    values <- with(my_list, eval(model_expr) - data_fit$value)

    return(values)
}



# evaluate_model_jacobian() -----------------------------------------------

#' Model jacobian evaluation method
#'
#' Method to evaluate the model jacobian with a given set of parameters.
#'
#' @param initial_parameters Parameters used for model jacobian evaluation
#'
#' @param par_list Named list with one entry per effect (with the respective
#' name). The entry contains a named vector with the unique stings of the
#' respective effect values i.e. the entries of the respective columns (e.g.
#' name, time, condition etc. for 'distinguish')
#' @param pass_parameter_list from the the \code{pass_parameter_list} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{effects_pars}}{
#'      Named list with one entry per effect. The values are the respective
#'      effect parameters.
#'  }
#'  \item{\code{model_derivertive_expr}}{
#'      The jacobian of the argument of the \code{model} parameter of
#'      \link{align_me} parsed as an executable expression.
#'  }
#' }
#'
#' @param pass_parameter_list2 from the the \code{pass_parameter_list2} argument
#' the following parameters are used:
#' \describe{
#'  \item{\code{data_fit}}{The current dataset.}
#' }
#'
#' @noRd
#'

evaluate_model_jacobian <- function(initial_parameters,
                                    par_list,
                                    pass_parameter_list,
                                    pass_parameter_list2) {
    effects_pars <- pass_parameter_list$effects_pars
    model_derivertive_expr <- pass_parameter_list$model_derivertive_expr

    data_fit <- pass_parameter_list2$data_fit


    distinguish <- seq_along(initial_parameters)
    my_list <- c(
        list(initial_parameters, distinguish = distinguish), as.list(data_fit),
        par_list
    )
    names(my_list)[1] <- effects_pars[[1]][1]
    derivative_values <- with(my_list, eval(model_derivertive_expr))

    return(derivative_values)
}



# ymaximal() --------------------------------------------------------------

#' Determine the maximal y values for plots
#'
#' A method producing a list of y values that will result in a eye pleasing
#' set of y tics when used in ggplot. Part of \link{plot_align_me}
#'
#' @param x numerical vector containing the y data which will be used to
#' determine the set of y max
#'
#' @return numerical vector with the maximal y values.
#'
#' @noRd
ymaximal <- function(x) {
    rx <- round(x, 1)
    lower <- floor(x)
    upper <- ceiling(x)
    lowdiff <- rx - lower
    uppdiff <- upper - rx

    if (x > 5) {
        if (upper %% 2 == 0) {
            out <- upper
        } else {
            out <- lower + 1
        }
    } else {
        if (x < 2) {
            if (rx > x) {
                if (rx %% 0.2 == 0) {
                    out <- rx
                } else {
                    out <- rx + 0.1
                }
            } else {
                if (rx %% 0.2 != 0) {
                    out <- rx + 0.1
                } else {
                    out <- rx + 0.2
                }
            }
        } else {
            if (lowdiff < uppdiff) {
                out <- lower + 0.5
            } else {
                out <- upper
            }
        }
    }

    return(out)
}


# scale_values() ----------------------------------------------------------

#' scale the values according to the current parameter scale
#'
#' @param parameter_fit_scale the current scale
#' @param value value to be scaled
#' @param upper upper error bound, i.e. value + sigma
#' @param lower lower error bound, i.e. value - sigma
#' @param sigma sigma to be scaled
#'
#' @noRd
scale_values <- function(parameter_fit_scale,
                         output_scale,
                         value,
                         upper,
                         lower,
                         sigma) {
    if (parameter_fit_scale == "log") {
        value <- exp(value)
        upper <- exp(upper)
        lower <- exp(lower)
        sigma <- value * sigma
    } else if (parameter_fit_scale == "log2") {
        value <- 2^(value)
        upper <- 2^(upper)
        lower <- 2^(lower)
        sigma <- value * sigma
    } else if (parameter_fit_scale == "log10") {
        value <- 10^(value)
        upper <- 10^(upper)
        lower <- 10^(lower)
        sigma <- value * sigma
    } else if (parameter_fit_scale == "linear") {
        value <- value
        upper <- upper
        lower <- lower
        sigma <- sigma
    }

    upper <- value + sigma
    lover <- value - sigma





    return(list(value = value, upper = upper, lower = lower, sigma = sigma))
}


#' # plot_time_course() ------------------------------------------------------
#'
#' #' Plotting subroutine for time course data
#' #'
#' #' The parameters are inherited from \link{plot_align_me}.
#' #'
#' #'
#' #' @return A ggplot update
#' #'
#' #' @noRd
#'
#' plot_time_course <- function(
#'     plot_list_points,
#'     plot_list_line,
#'     plot_points,
#'     plot_line,
#'     spline,
#'     scales,
#'     align_zeros,
#'     ncol,
#'     my_colors,
#'     xlab,
#'     ylab
#' ) {
#'     if (is.null(xlab)){
#'         x_label <- "Time"
#'     } else {
#'         x_label <- xlab
#'     }
#'
#'     if (is.null(ylab)){
#'         y_label <- "Signal"
#'     } else {
#'         y_label <- ylab
#'     }
#'     ## plot
#'     if (plot_points == "aligned" & plot_line == "aligned") {
#'         g <- ggplot(
#'             data = plot_list_points,
#'             aes(
#'                 x = time,
#'                 y = value,
#'                 group = distinguish,
#'                 color = distinguish,
#'                 fill = distinguish
#'             )
#'         )
#'         g <- g + facet_wrap(~name, scales = scales, ncol = ncol)
#'         if (is.null(my_colors)) {
#'             # my_colors <- scale_color_brewer()
#'             # # c(
#'             # # "#000000",
#'             # # "#C5000B",
#'             # # "#0084D1",
#'             # # "#579D1C",
#'             # # "#FF950E",
#'             # # "#4B1F6F",
#'             # # "#CC79A7",
#'             # # "#006400",
#'             # # "#F0E442",
#'             # # "#8B4513",
#'             # # rep("gray", 100)
#'             # # )
#'             # g <- g + scale_color_manual("Condition", values = my_colors) +
#'             #     scale_fill_manual("Condition", values = my_colors)
#'         } else {
#'             my_colors <- c(my_colors, rep("gray", 100))
#'             g <- g + scale_color_manual("Distinguished", values = my_colors) +
#'                 scale_fill_manual("Distinguished", values = my_colors)
#'         }
#'     } else {
#'         g <- ggplot(
#'             data = plot_list_points,
#'             aes(
#'                 x = time,
#'                 y = value,
#'                 group = scaling,
#'                 color = scaling,
#'                 fill = scaling
#'             )
#'         )
#'         g <- g + facet_wrap(
#'             ~ name * distinguish,
#'             scales = scales,
#'             ncol = ncol
#'         )
#'         if (is.null(my_colors)) {
#'             # my_colors <- scale_color_brewer()
#'             # # c(
#'             # # "#000000",
#'             # # "#C5000B",
#'             # # "#0084D1",
#'             # # "#579D1C",
#'             # # "#FF950E",
#'             # # "#4B1F6F",
#'             # # "#CC79A7",
#'             # # "#006400",
#'             # # "#F0E442",
#'             # # "#8B4513",
#'             # # rep("gray", 100)
#'             # # )
#'             # g <- g + scale_color_manual("Scaling", values = my_colors) +
#'             #     scale_fill_manual("Scaling", values = my_colors)
#'         } else {
#'             my_colors <- c(my_colors, rep("gray", 100))
#'             g <- g + scale_color_manual("Scaled", values = my_colors) +
#'                 scale_fill_manual("Scaled", values = my_colors)
#'         }
#'     }
#'
#'
#'     if (spline) {
#'         g <- g + geom_errorbar(
#'             data = plot_list_line,
#'             aes(
#'                 ymin = lower, # value - sigma,
#'                 ymax = upper # value + sigma
#'             ),
#'             width = 0
#'         )
#'         g <- g + geom_smooth(
#'             data = plot_list_line,
#'             se = FALSE,
#'             method = "lm",
#'             formula = y ~ poly(x, 3)
#'         )
#'     } else {
#'         if (plot_points == plot_line | plot_line == "prediction") {
#'             g <- g + geom_line(data = plot_list_line, size = 1)
#'             if (plot_line == "prediction") {
#'                 g <- g + geom_ribbon(
#'                     data = plot_list_line,
#'                     aes(
#'                         # probably this should be changed back
#'                         ymin = value - sigma, # lower
#'                         ymax = value + sigma # upper
#'                     ),
#'                     alpha = 0.1,
#'                     lty = 0
#'                 )
#'             }
#'         } else {
#'             g <- g + geom_line(data = plot_list_line, size = 1)#, color = "grey")
#'             g <- g + geom_ribbon(
#'                 data = plot_list_line,
#'                 aes(
#'                     ymin = lower, # value - sigma,
#'                     ymax = upper, # value + sigma#,
#'                     # fill = "grey",
#'                     # color = "grey"
#'                 ),
#'                 alpha = 0.3,
#'                 lty = 0
#'             )
#'         }
#'     }
#'
#'     g <- g + geom_point(data = plot_list_points, size = 2.5)
#'
#'     if(plot_line != "prediction"){
#'         errwidth <- max(plot_list_points$time)/50
#'         g <- g + geom_errorbar(
#'             data = plot_list_points,
#'             aes(
#'                 ymin = lower, # value - sigma,
#'                 ymax = upper # value + sigma
#'             ),
#'             size = 0.5,
#'             width = errwidth,
#'             alpha = 0.5
#'         )
#'     } else {
#'         errwidth <- max(plot_list_points$time)/50
#'         g <- g + geom_errorbar(
#'             data = plot_list_points,
#'             aes(
#'                 ymin = value - sigma, # ,
#'                 ymax = value + sigma #
#'             ),
#'             size = 0.5,
#'             width = errwidth,
#'             alpha = 0.5
#'         )
#'     }
#'
#'
#'     g <- g + theme_bw(base_size = 20) +
#'         theme(
#'             legend.position = "top",
#'             legend.key = element_blank(),
#'             strip.background = element_rect(color = NA, fill = NA),
#'             axis.line.x = element_line(size = 0.3, colour = "black"),
#'             axis.line.y = element_line(size = 0.3, colour = "black"),
#'             panel.grid.major.x = element_blank(),
#'             panel.grid.major.y = element_blank(),
#'             panel.grid.minor = element_blank(),
#'             panel.border = element_blank(),
#'             panel.background = element_blank(),
#'             plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
#'         )
#'     g <- g + xlab("\nTime") + ylab("Signal\n")
#'
#'     if (align_zeros) {
#'         if (plot_points != "original") {
#'             # scale y-axes (let them start at same minimum determined by
#'             # smallest value-sigma and end at individual ymax)
#'             plot_list_points <- as.data.table(plot_list_points)
#'             blank_data <- plot_list_points[
#'                 ,
#'                 list(ymax = max(upper), ymin = min(lower)),
#'                 by = c("name", "distinguish", "scaling")
#'             ]
#'             blank_data[, ":="(ymin = min(ymin))] # same minimum for all proteins
#'             blank_data[
#'                 ,
#'                 ":="(ymax = ymaximal(ymax)),
#'                 by = c("name", "distinguish", "scaling")
#'             ] # protein specific maximum
#'             blank_data <- melt(
#'                 blank_data,
#'                 id.vars = c("name", "distinguish", "scaling"),
#'                 measure.vars = c("ymax", "ymin"),
#'                 value.name = "value"
#'             )
#'             blank_data[, ":="(time = 0, variable = NULL)]
#'             g <- g + geom_blank(
#'                 data = as.data.frame(blank_data),
#'                 aes(x = time, y = value)
#'             )
#'         }
#'     }
#'
#'     return(g)
#' }
#'
#'
#' # plot_dose_response() ----------------------------------------------------
#'
#' #' Plotting subroutine for dose response data
#' #'
#' #' The parameters are inherited from \link{plot_align_me}.
#' #'
#' #'
#' #' @return A ggplot update
#' #'
#' #' @noRd
#'
#' plot_dose_response <- function(
#'     plot_list_points,
#'     plot_list_line,
#'     plot_points,
#'     plot_line,
#'     spline,
#'     scales,
#'     align_zeros,
#'     ncol,
#'     my_colors,
#'     xlab,
#'     ylab
#' ) {
#'     if (is.null(xlab)){
#'         x_label <- "Dose"
#'     } else {
#'         x_label <- xlab
#'     }
#'
#'     if (is.null(ylab)){
#'         y_label <- "Signal"
#'     } else {
#'         y_label <- ylab
#'     }
#'
#'
#'
#'
#'     if (!("dose" %in% names(plot_list_points))) {
#'         stop("'dose' must be set as a 'distinguish' parameter in align_me()\n")
#'     }
#'
#'     ## plot
#'     if (plot_points == "aligned" & plot_line == "aligned") {
#'         g <- ggplot(
#'             data = plot_list_points,
#'             aes(
#'                 x = dose,
#'                 y = value,
#'                 group = distinguish,
#'                 color = distinguish,
#'                 fill = distinguish
#'             )
#'         )
#'         g <- g + facet_wrap(~name, scales = scales, ncol = ncol)
#'         if (is.null(my_colors)) {
#'             # my_colors <- scale_color_brewer()
#'             # # c(
#'             # # "#000000",
#'             # # "#C5000B",
#'             # # "#0084D1",
#'             # # "#579D1C",
#'             # # "#FF950E",
#'             # # "#4B1F6F",
#'             # # "#CC79A7",
#'             # # "#006400",
#'             # # "#F0E442",
#'             # # "#8B4513",
#'             # # rep("gray", 100)
#'             # # )
#'             # g <- g + scale_color_manual("Condition", values = my_colors) +
#'             #     scale_fill_manual("Condition", values = my_colors)
#'         } else {
#'             my_colors <- c(my_colors, rep("gray", 100))
#'             g <- g + scale_color_manual("Distinguished", values = my_colors) +
#'                 scale_fill_manual("Distinguished", values = my_colors)
#'         }
#'     } else {
#'         g <- ggplot(
#'             data = plot_list_points,
#'             aes(
#'                 x = dose,
#'                 y = value,
#'                 group = scaling,
#'                 color = scaling,
#'                 fill = scaling
#'             )
#'         )
#'         g <- g + facet_wrap(
#'             ~ name * distinguish,
#'             scales = scales,
#'             ncol = ncol
#'         )
#'         if (is.null(my_colors)) {
#'             # my_colors <- scale_color_brewer()
#'             # # c(
#'             # # "#000000",
#'             # # "#C5000B",
#'             # # "#0084D1",
#'             # # "#579D1C",
#'             # # "#FF950E",
#'             # # "#4B1F6F",
#'             # # "#CC79A7",
#'             # # "#006400",
#'             # # "#F0E442",
#'             # # "#8B4513",
#'             # # rep("gray", 100)
#'             # # )
#'             # g <- g + scale_color_manual("Scaling", values = my_colors) +
#'             #     scale_fill_manual("Scaling", values = my_colors)
#'         } else {
#'             my_colors <- c(my_colors, rep("gray", 100))
#'             g <- g + scale_color_manual("Scaled", values = my_colors) +
#'                 scale_fill_manual("Scaled", values = my_colors)
#'         }
#'     }
#'
#'
#'     if (spline) {
#'         g <- g + geom_errorbar(
#'             data = plot_list_line,
#'             aes(
#'                 ymin = lower, # value - sigma,
#'                 ymax = upper # value + sigma
#'             ),
#'             width = 0
#'         )
#'         g <- g + geom_smooth(
#'             data = plot_list_line,
#'             se = FALSE,
#'             method = "lm",
#'             formula = y ~ poly(x, 3)
#'         )
#'     } else {
#'         if (plot_points == plot_line | plot_line == "prediction") {
#'             g <- g + geom_line(data = plot_list_line, size = 1)
#'             if (plot_line == "prediction") {
#'                 g <- g + geom_ribbon(
#'                     data = plot_list_line,
#'                     aes(
#'                         ymin = lower, # value - sigma,
#'                         ymax = upper # value + sigma
#'                     ),
#'                     alpha = 0.1,
#'                     lty = 0
#'                 )
#'             }
#'         } else {
#'             g <- g + geom_line(data = plot_list_line, size = 1)#, color = "grey")
#'             g <- g + geom_ribbon(
#'                 data = plot_list_line,
#'                 aes(
#'                     ymin = lower, # value - sigma,
#'                     ymax = upper#, # value + sigma,
#'                     # fill = "grey",
#'                     # color = "grey"
#'                 ),
#'                 alpha = 0.3,
#'                 lty = 0
#'             )
#'         }
#'     }
#'
#'     g <- g + geom_point(data = plot_list_points, size = 2.5)
#'     errwidth <- max(plot_list_points$dose)/50
#'     g <- g + geom_errorbar(
#'         data = plot_list_points,
#'         aes(
#'             ymin = lower, # value - sigma,
#'             ymax = upper # value + sigma
#'         ),
#'         size = 0.5,
#'         width = errwidth,
#'         alpha = 0.5
#'     )
#'
#'     g <- g + theme_bw(base_size = 20) +
#'         theme(
#'             legend.position = "top",
#'             legend.key = element_blank(),
#'             strip.background = element_rect(color = NA, fill = NA),
#'             axis.line.x = element_line(size = 0.3, colour = "black"),
#'             axis.line.y = element_line(size = 0.3, colour = "black"),
#'             panel.grid.major.x = element_blank(),
#'             panel.grid.major.y = element_blank(),
#'             panel.grid.minor = element_blank(),
#'             panel.border = element_blank(),
#'             panel.background = element_blank(),
#'             plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
#'         )
#'     g <- g + xlab(paste0("\n",x_label)) + ylab(paste0(y_label,"\n"))
#'
#'     if (align_zeros) {
#'         if (plot_points != "original") {
#'             # scale y-axes (let them start at same minimum determined by
#'             # smallest value-sigma and end at individual ymax)
#'             plot_list_points <- as.data.table(plot_list_points)
#'             blank_data <- plot_list_points[
#'                 ,
#'                 list(ymax = max(upper), ymin = min(lower)),
#'                 by = c("name", "distinguish", "scaling")
#'             ]
#'             blank_data[, ":="(ymin = min(ymin))] # same minimum for all proteins
#'             blank_data[
#'                 ,
#'                 ":="(ymax = ymaximal(ymax)),
#'                 by = c("name", "distinguish", "scaling")
#'             ] # protein specific maximum
#'             blank_data <- melt(
#'                 blank_data,
#'                 id.vars = c("name", "distinguish", "scaling"),
#'                 measure.vars = c("ymax", "ymin"),
#'                 value.name = "value"
#'             )
#'             blank_data[, ":="(dose = 0, variable = NULL)]
#'             g <- g + geom_blank(
#'                 data = as.data.frame(blank_data),
#'                 aes(x = dose, y = value)
#'             )
#'         }
#'     }
#'
#'     return(g)
#' }
