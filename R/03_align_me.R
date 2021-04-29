
#' Align time-course data based on an Mixed-Effects alignment model
#'
#' The function deals primarily with time-course data of different
#' targets which have been measured under different experimental
#' conditions and whose measured values might be on a different
#' scale, e.g. because of different amplification. The algorithm
#' determines the different scaling and estimates the time-course on
#' a common scale.
#'
#'
#' @param data data.frame containing the the data to be scaled. Usualy the
#' output of \link{read_wide} is used. Obligatory are the columns \code{name},
#' and \code{value}. In the case of time dependet data, a \code{time} column is
#' also necessary. Additionally, \code{data} should have further columns,
#' e.g. characterizing experimental conditions (the fixed effects) and
#' sources of variance between data of the same experimental condition
#' (the scaled variables).
#'
#' @param model character defining the model by which the values in
#' \code{data} can be described, e.g. "yi/sj"
#'
#' @param error_model character defining a model for the standard
#' deviation of a value, e.g. "sigma0 + value * sigmaR". This model
#' can contain parameters, e.g. "sigma0" and "sigmaR", or numeric
#' variables from \code{data}, e.g. "value" or "time".
#'
#' @param distinguish two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "yi", and "name1, ..."
#' refers to variables in \code{data}, e.g. "condition".
#' The parameters "par1, ..." are determined specific to the levels of
#' "name1, ...".
#'
#' @param scaling two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "sj", and "name1, ..."
#' refers to variables in \code{data}, e.g. "Experiment".
#' @param error two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameter contained in \code{error}, e.g. "sigma0" and "sigmaR", and
#' "name1, ..." refers to variables in \code{data}. If the same values
#' of "par1, ..." should be assumed for all data, "name1" can be "1".
#' @param input_scale character, describes how the input is scaled. Must be one
#' of \code{c("linear", "log", "log2", "log10")}.
#' @param normalize logical indicating whether the distinguishing effect
#' parameter should be normalized to unit mean.
#'
#' @param average_techn_rep logical, indicates if the technical replicates
#' should be averaged
#'
#' @param names_as_factors logical, indicates if the \code{name} columns of the
#' output should be parsed as factors instead of characters.
#'
#' @param verbose logical, print out information about each fit
#' @param normalize_input logical, if TRUE the input will be normalized before
#' scaling. see \code{split_data}.
#' @details Alignment of time-course data is achieved by an alignment
#' model which explains the observed data by a function mixing
#' fixed effects, usually parameters reflecting the "underlying"
#' time-course, and latent variables, e.g. scaling parameters taking
#' account for effects like different amplification or loading, etc.
#' Depending on the measurement technique, the data has constant
#' relative error, or constant absolute error or even a combination
#' of those. This error is described by an error function. The error
#' parameters are usually global, i.e. the same parameter values are
#' assumed for all data points.
#'
#' @return Object of class \code{aligned}, i.e. a data frame of the
#' alignment result containing an attribute "outputs":
#'  a list of data frames
#' \describe{
#' \item{aligned}{data.frame with the original column names plus the column
#'     \code{sigma}. Each set of unique distinguish effects i.e. biological
#'     different condition (e.g. time point, target and treatment) has one set
#'     of \code{value} and \code{sigma}. The values are the estimated true
#'     values i.e. the determined distinguish parameters. The errors in the
#'     \code{sigma} column are estimated by employing the fisher information
#'     to quantify the uncertainty of the respective fit. Both, the value and
#'     its error are on the common scale.}
#' \item{scaled}{The original measurements scaled to common scale by applying
#'     inverse model. The errors are the result of the evaluation of the error
#'     model and then also scaled to common scale by use of Gaussian error
#'     propagation.}
#' \item{prediction}{Original data with \code{value} replaced by the prediction
#'     (evaluation of the model with the fitted parameters), and \code{sigma}
#'     from the evaluation of the error model. Both are on the original scale.
#'     }
#'
#' \item{original}{The original data as passed as  \code{data}.}
#' \item{original_with_parameters}{The original data but with added columns
#'     containing the estimated parameters}
#'
#' \item{parameter}{Parameter table with the columns: \code{name}: the name of
#'     the current target, \code{level}: the pasted unique set of effects
#'     (distinguish, scaling or error), \code{parameter}: the parameter
#'     identifier as defined in the (error) model, \code{value} and \code{sigma}
#'     containing the determined values and corresponding errors, \code{nll}:
#'     twice the negative log-likelihood of the fit, \code{no_pars} and
#'     \code{no_data} containing the number of parameters and data points for
#'     the respected fit. This list entry also has two attributes: \code{value}
#'     containing the final value (residual sum of squares) passed to
#'     \link{trust::trust} by the objective function and \code{df} the degrees
#'     of freedom of the fitting process.}
#'
#' \item{distinguish}{Names of the columns containing the distinguish effects}
#' \item{scaling}{Names of the columns containing the scaling effects}
#' \item{error}{Names of the columns containing the error effects}
#' }
#'
#' The estimated parameters are returned by the attribute "parameters".
#' @example inst/examples/example_align_me.R
#' @seealso \link{read_wide} to read data in a wide column format and
#' get it in the right format for \code{align_me}.
#'
#' @importFrom stats D
#'
#' @export
align_me <- function(data,
                     model = NULL,
                     error_model = NULL,
                     distinguish = NULL,
                     scaling = NULL,
                     error = NULL,
                     input_scale = "linear",
                     normalize = TRUE,
                     average_techn_rep = FALSE,
                     verbose = FALSE,
                     names_as_factors = TRUE,
                     normalize_input = TRUE) {
    if (FALSE) {
        if (FALSE) {
            sim_data_wide_file <- system.file(
                "extdata", "sim_data_wide.csv",
                package = "blotIt3"
            )
            data <- read_wide(sim_data_wide_file, description = seq_len(3))
        } else {
            data <- read_wide(
                file = paste0(
                    "/home/severin/Documents/PhD/Projects/R/WesternblotDataSim",
                    "/data/2021-03-04_AS-E2_cells_Time_course_data.csv"
                ),
                description = 1:7,
                sep = ",",
                dec = "."
            )
            data <- data[c(1, 2, 7, 8, 9)]
            names(data) <- c("time", "condition", "ID", "name", "value")

            data <- subset(data, name == "pEPOR_au")
        }

        model <- "yi / sj"
        error_model <- "value * sigmaR"
        distinguish <- yi ~ name + time + condition
        scaling <- sj ~ name + ID
        error <- sigmaR ~ name + 1
        input_scale <- "linear"
        normalize <- TRUE
        average_techn_rep <- FALSE
        names_as_factors <- TRUE
        verbose <- TRUE
        normalize_input <- TRUE
    }



    # check input for mistakes
    input_check_report <- input_check(
        data = data,
        model = model,
        error_model = error_model,
        distinguish = distinguish,
        scaling = scaling,
        error = error,
        input_scale = input_scale
    )

    if (verbose) {
        cat(input_check_report, "\n")
    }

    # Check if data is already blotIt output
    parameter_data <- NULL

    if (inherits(data, "aligned")) {
        parameter_data <- data$parameters
        data <- data$original
    }

    ## read distinguishing, scaling and error effects from input
    effects <- identify_effects(
        distinguish = distinguish,
        scaling = scaling,
        error = error
    )

    effects_values <- effects$effects_values
    effects_pars <- effects$effects_pars

    ## prepare data
    data <- as.data.frame(data)

    to_be_scaled <- split_for_scaling(
        data,
        effects_values,
        normalize_input,
        input_scale
    )

    # Generate unique list of targets
    targets <- make.unique(
        vapply(to_be_scaled, function(d) as.character(d$name)[1],""),
        sep = "_"
    )

    # Rename names to unique if multiple scalings per target exist
    for (n in seq_along(targets)) {
        to_be_scaled[[n]]$name <- targets[n]
    }

    ## Get distinguish, scaling and error parameter from model and error model
    parameters <- get_symbols(c(model, error_model), exclude = colnames(data))

    # Include the normalization term as a constraint if saied so in the function
    # call
    if (normalize) {
        constraint <- paste("1e3 * (mean(", effects_pars[1][1], ") - 1)")
        c_strength <- 1000
    } else {
        constraint <- "0"
        c_strength <- 0
    }

    # Retrieve the covariates as the "remaining" model parameters, when fixed,
    # latent and errorparameters are excluded
    covariates <- union(
        get_symbols(
            model,
            exclude = c(effects_pars[1], effects_pars[2], effects_pars[3])
        ),
        get_symbols(
            error_model,
            exclude = c(effects_pars[1], effects_pars[2], effects_pars[3])
        )
    )
    cat("Covariates:", paste(covariates, sep = ", "), "\n")

    # Check if parameters from (error-) model and passed expressions coincide
    if (
        length(
            setdiff(
                c(
                    effects_pars[1],
                    effects_pars[2],
                    effects_pars[3]
                ),
                parameters
            )
        ) > 0
    ) {
        stop("Not all paramters are defined in either arguments
         'scaling', 'distinguish' or 'error'")
    }

    # Name the respective parameters fixed, latent and error
    names(parameters)[parameters %in% effects_pars[1]] <- "distinguish"
    names(parameters)[parameters %in% effects_pars[2]] <- "scaling"
    names(parameters)[parameters %in% effects_pars[3]] <- "error"

    # parse error model by replacing the "value" by the model
    error_model <- replace_symbols(
        "value",
        paste0("(", model, ")"),
        error_model
    )

    # Apply transformation according to the given 'input_scale' parameter
    if (input_scale == "log") {
        model <- replace_symbols(parameters, paste0(
            "exp(", parameters,
            ")"
        ), model)
        error_model <- replace_symbols(parameters, paste0(
            "exp(",
            parameters, ")"
        ), error_model)
        constraint <- replace_symbols(parameters, paste0(
            "exp(",
            parameters, ")"
        ), constraint)
        if (verbose) {
            cat("model, errormodel and constraint are scaled: x <- exp(x)\n")
        }
    } else if (input_scale == "log2") {
        model <- replace_symbols(parameters, paste0(
            "2^(", parameters,
            ")"
        ), model)
        error_model <- replace_symbols(parameters, paste0(
            "2^(",
            parameters, ")"
        ), error_model)
        constraint <- replace_symbols(parameters, paste0(
            "2^(",
            parameters, ")"
        ), constraint)
        if (verbose) {
            cat("model, errormodel and constraint are scaled: x <- 2^(x)\n")
        }
    } else if (input_scale == "log10") {
        model <- replace_symbols(parameters, paste0(
            "10^(", parameters,
            ")"
        ), model)
        error_model <- replace_symbols(parameters, paste0(
            "10^(",
            parameters, ")"
        ), error_model)
        constraint <- replace_symbols(parameters, paste0(
            "10^(",
            parameters, ")"
        ), constraint)
        if (verbose) {
            cat("model, errormodel and constraint are scaled: x <- 10^(x)\n")
        }
    } else if (input_scale == "linear") {
        if (verbose) {
            cat("model, errormodel and constraint remain linear scaled.\n")
        }
    }


    # Calculating the derivative
    model_derivertive <- deparse(
        D(parse(text = model), name = effects_pars[[1]][1])
    )

    # Calculating the (error) model jacobian
    model_jacobian <- lapply(
        parameters,
        function(p) {
            deparse(D(parse(text = model), name = p))
        }
    )
    error_model_jacobian <- lapply(
        parameters,
        function(p) {
            deparse(D(parse(text = error_model), name = p))
        }
    )

    # Construct math function expressions
    constraint_expr <- parse(text = constraint)

    # Replace the parameters used in the (error) model and the derivatives
    # by placeholders of the respective effect category
    for (n in seq_along(parameters)) {
        model <- replace_symbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), model)
        error_model <- replace_symbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), error_model)
        model_derivertive <- replace_symbols(parameters[n], paste0(
            parameters[n],
            "[", names(parameters)[n], "]"
        ), model_derivertive)
        model_jacobian <- lapply(model_jacobian, function(myjac) {
            replace_symbols(
                parameters[n],
                paste0(
                    parameters[n], "[", names(parameters)[n],
                    "]"
                ), myjac
            )
        })
        error_model_jacobian <- lapply(error_model_jacobian, function(myjac) {
            replace_symbols(
                parameters[n],
                paste0(
                    parameters[n], "[", names(parameters)[n],
                    "]"
                ), myjac
            )
        })
    }

    cat("Model:        ", model, "\n", sep = "")
    cat("Error Model:  ", error_model, "\n", sep = "")

    # Parse functions to an executable expression
    model_expr <- parse(text = model)
    model_derivertive_expr <- parse(text = model_derivertive)
    error_model_expr <- parse(text = error_model)

    model_jacobian_expr <- lapply(
        model_jacobian,
        function(myjac) parse(text = myjac)
    )
    error_model_jacobian_expr <- lapply(
        error_model_jacobian,
        function(myjac) parse(text = myjac)
    )

    pass_parameter_list <- list(
        effects_values = effects_values,
        parameter_data = parameter_data,
        average_techn_rep = average_techn_rep,
        verbose = verbose,
        covariates = covariates,
        parameters = parameters,
        input_scale = input_scale,
        effects_pars = effects_pars,
        model_expr = model_expr,
        error_model_expr = error_model_expr,
        constraint_expr = constraint_expr,
        model_jacobian_expr = model_jacobian_expr,
        error_model_jacobian_expr = error_model_jacobian_expr,
        c_strength = c_strength,
        normalize = normalize,
        model_derivertive_expr = model_derivertive_expr
    )
    out <- lapply(
        seq_along(to_be_scaled),
        function(i,
                 pass_parameter_list) {
            try(
                {
                    cat("Target ", i, "/", length(to_be_scaled), ":", targets[i], "\n", sep = "")
                    out <- scale_target(
                        current_data = to_be_scaled[[i]],
                        pass_parameter_list = pass_parameter_list
                    )

                    return(out)
                },
                silent = FALSE
            )
        },
        # additional parameters for scale_target
        pass_parameter_list = pass_parameter_list
    )

    # Sanitize results from failed fits
    # rbind parameter table from not-failed elements of out
    parameter_table <- do.call(
        rbind,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    o[[6]]
                } else {
                    NULL
                }
            }
        )
    )

    # add up the "values", i.e. the -2*LL
    attr(parameter_table, "value") <- do.call(
        sum,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    attr(o[[6]], "value")
                } else {
                    0
                }
            }
        )
    )

    # add up degrees of freedom
    attr(parameter_table, "df") <- do.call(
        sum,
        lapply(
            out,
            function(o) {
                if (!inherits(o, "try-error")) {
                    attr(o[[6]], "df")
                } else {
                    0
                }
            }
        )
    )

    # paste together the other results
    out_combined <- lapply(
        seq_len(5),
        function(i) {
            do.call(
                rbind,
                lapply(
                    out,
                    function(o) {
                        if (!inherits(o, "try-error")) {
                            o[[i]]
                        } else {
                            NULL
                        }
                    }
                )
            )
        }
    )

    names(out_combined) <- c(
        "prediction",
        "scaled",
        "aligned",
        "original",
        "original_with_parameters"
    )

    return_list <- list(
        aligned = out_combined$aligned,
        scaled = out_combined$scaled,
        prediction = out_combined$prediction,
        original = out_combined$original,
        original_with_parameters = out_combined$original_with_parameters,
        parameter = parameter_table,
        distinguish = effects_values[[1]],
        scaling = effects_values[[2]],
        error = effects_values[[3]]
    )

    if (names_as_factors == TRUE) {
        return_list$aligned$name <- factor(
            return_list$aligned$name,
            levels = unique(return_list$aligned$name)
        )
        return_list$scaled$name <- factor(
            return_list$scaled$name,
            levels = unique(return_list$scaled$name)
        )
        return_list$prediction$name <- factor(
            return_list$prediction$name,
            levels = unique(return_list$prediction$name)
        )
        return_list$original$name <- factor(
            return_list$original$name,
            levels = unique(return_list$original$name)
        )
        return_list$original_with_parameters$name <- factor(
            return_list$original_with_parameters$name,
            levels = unique(return_list$original_with_parameters$name)
        )
    }

    # class(return_list) <- c("aligned", "data.frame")

    return(return_list)
}
