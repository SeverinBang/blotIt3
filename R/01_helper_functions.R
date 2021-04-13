
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
