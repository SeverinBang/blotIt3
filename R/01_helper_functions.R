
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



# identify_parameters() ---------------------------------------------------

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

identify_parameters <- function(
    distinguish = NULL, scaling = NULL, error = NULL
    ) {

    if (length(as.character(distinguish)) == 1 |
        length(as.character(scaling)) == 1 | length(as.character(error)) == 1) {
        stop("Do not pass distinguish, scaling or error as string.")
    }

    # Stop if formulas have the wrong specification
    if (length(as.character(distinguish)) < 3) {
        stop("Left and right-hand side of formula 'distinguish' is needed")
    }
    if (length(as.character(scaling)) < 3) {
        stop("Left and right-hand side of formula 'scaling' is needed")
    }
    if (length(as.character(error)) < 3) {
        stop("Left and right-hand side of formula 'error' is needed")
    }

    # Get distinguish scale and error parameters
    distinguish <- get_symbols(as.character(distinguish)[3])
    scaling <- get_symbols(as.character(scaling)[3])
    error <- get_symbols(as.character(error)[3])

    # Determine to which class parameters belong
    distinguish_pars <- get_symbols(as.character(distinguish)[2])
    scaling_pars <- get_symbols(as.character(scaling)[2])
    error_pars <- get_symbols(as.character(error)[2])
}
