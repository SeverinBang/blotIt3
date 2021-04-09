
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
#' @param input_is_log logical indicating whether input data is on log scale.
#' @param normalize logical indicating whether the distinguishing effect
#' parameter should be normalized to unit mean.
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
#' \item{prediction}{original data with value and sigma replaced by
#'                   the predicted values and sigmas}
#' \item{scaled}{original data with the values transformed according
#'               to the inverse model, i.e. \code{model} solved for
#'               the first parameter in \code{fixed}, e.g. "ys".
#'               Sigma values are computed by error propagation
#'               from the inverse model equation.}
#' \item{aligned}{the reduced data with the fixed effects and their
#'                uncertainty, only. The result of the alignment
#'                algorithm.}
#' \item{original}{the original data}
#' \item{parameter}{original data augmented by parameter columns.
#'                  Parameters in each row correspond to the levels of
#'                  fixed, latent or error as passed to \code{align_me()}.
#'                  Used for initialization or parameter values when
#'                  refitting with modified model.}
#' }
#'
#' The estimated parameters are returned by the attribute "parameters".
#' @example inst/examples/example_align_me.R
#' @seealso \link{read_wide} to read data in a wide column format and
#' get it in the right format for \code{align_me()}.
#' @export
align_me <- function(data,
                     model = NULL,
                     error_model = NULL,
                     distinguish = NULL,
                     scaling = NULL,
                     error = NULL,
                     input_is_log = FALSE,
                     normalize = TRUE,
                     verbose = TRUE,
                     normalize_input = TRUE) {
    if (FALSE) {
        sim_data_wide_file <- system.file(
            "extdata", "sim_data_wide.csv",
            package = "blotIt3"
        )
        data <- read_wide(sim_data_wide_file, description = seq_len(3))
        model <- "yi / sj"
        error_model <- "value * sigmaR"
        distinguish <- yi ~ name + time + condition
        scaling <- sj ~ name + ID
        error <- sigmaR ~ name + 1
        input_is_log <- FALSE
        normalize <- TRUE
        verbose <- TRUE
        normalize_input <- TRUE
    }

    if (is.null(model) | is.null(error_model) | is.null(distinguish) |
        is.null(scaling) | is.null(error)) {
        stop(
            "All of model, error_model, distinguish, scaling, error ",
            "must be set."
        )
    }

    ## read distinguishing, scaling and error effects from input
    effects <- identify_effects(
        distinguish = distinguish,
        scaling = scaling,
        error = error
    )

    scaling_values <- effects$effects_values$scaling_values
    distinguish_values <- effects$effects_values$distinguish_values
    error_values <- effects$effects_values$error_values

    scaling_pars <- effects$effects_pars$scaling_pars
    distinguish_pars <- effects$effects_pars$distinguish_pars
    error_pars <- effects$effects_pars$error_pars

    ## prepare data
    data <- as.data.frame(data)

    targets <- unique(data$name)

    to_be_scaled <- split_for_scaling(
        data,
        distinguish_values,
        scaling_values,
        normalize_input,
        input_is_log
    )

    ## normalize input per scaling group to avoid failure of scaling. Only for
    ## linear data.
    if (normalize_input & !input_is_log) {
        to_be_scaled <- lapply(
            to_be_scaled,
            function(i) {
                i$value <- i$value / mean(i$value)

                return(i)
            }
        )
    }

    ## Get distinguish, scaling and error parameter from model and error model
    parameters <- get_symbols(c(model, error_model), exclude = colnames(data))

    # Include the normalization term as a constraint if saied so in the function
    # call
    if (normalize) {
        constraint <- paste("1e3 * (mean(", distinguish_pars[1], ") - 1)")
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
            exclude = c(distinguish_pars, scaling_pars, error_pars)
        ),
        get_symbols(
            error_model,
            exclude = c(distinguish_pars, scaling_pars, error_pars)
        )
    )
    # cat("Covariates:", paste(covariates, sep = ", "), "\n")

    # Check if parameters from (error-) model and passed expressions coincide
    if (
        length(
            setdiff(c(distinguish_pars, scaling_pars, error_pars), parameters)
        ) > 0
    ) {
        stop("Not all paramters are defined in either arguments
         'scaling', 'distinguish' or 'error'")
    }

    # Name the respective parameters fixed, latent and error
    names(parameters)[parameters %in% distinguish_pars] <- "distinguish"
    names(parameters)[parameters %in% scaling_pars] <- "scaling"
    names(parameters)[parameters %in% error_pars] <- "error"

    # parse error model by replacing the "value" by the model
    error_model <- replace_symbols(
        "value",
        paste0("(", model, ")"),
        error_model
    )

    # Calculating the derivative
    model_derivertive <- deparse(
        D(parse(text = model), name = distinguish_pars[1])
    )
}
