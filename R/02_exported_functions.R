# read_wide() -------------------------------------------------------------

#' Import measurement data
#'
#' Data from a given wide style csv-file is imported. While importing, the data
#' is converted into a long table
#'
#' @param file character, the name of the data file.
#' @param description numeric index vector of the columns containing
#' the description.
#' @param time numeric index of length 1: the time column.
#' @param header a logical value indicating whether the file contains
#' the names of the variables as its first line.
#' @param ... further arguments being passed to \link{read.csv}.
#'
#' @return data frame with columns "name", "time", "value" and other
#' columns describing the measurements.
#'
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' ## Import example data set
#' sim_data_wide_file <- system.file(
#'     "extdata", "sim_data_wide.csv",
#'     package = "blotIt3"
#' )
#' read_wide(sim_data_wide_file, description = seq_len(3))
read_wide <- function(file, description = NULL, time = 1, header = TRUE, ...) {
    my_data <- read.csv(file, header = header, ...)
    all_names <- colnames(my_data)

    if (is.null(description)) {
        stop("Specify columns containing descriptions.")
    }
    if (is.character(description)) {
        n_description <- which(colnames(my_data) %in% description)
    } else {
        n_description <- description
    }
    ## Check the index of the "time" column
    if (is.character(time)) {
        n_time <- which(colnames(my_data) == time)
    } else {
        n_time <- time
    }

    ## Check availability of description and time
    if (length(n_description) < length(description)) {
        warning(
            "Not all columns proposed by argument 'description' are available",
            " in file.\nTaking the available ones."
        )
    }

    if (length(n_time) == 0) {
        stop(
            "File did not contain a time column as proposed by 'time' argument."
        )
    }

    # Distinguish description data from measurement data
    description_entries <- my_data[, n_description]
    rest_long <- unlist(my_data[, -n_description])

    # Create output data frame
    new_data <- data.frame(
        description_entries,
        name = rep(all_names[-n_description], each = dim(my_data)[1]),
        value = rest_long
    )

    # Remove missing items
    new_data <- new_data[!is.nan(new_data$value), ]
    new_data <- new_data[!is.na(new_data$value), ]


    colnames(new_data)[n_time] <- "time"

    return(new_data)
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
                              input_scale) {
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
        if (input_scale == "linear") {
            for (i in seq_len(length(list_out))) {
                list_out[[i]]$value <- list_out[[i]]$value /
                    mean(list_out[[i]]$value)
            }
        } else {
            warning(
                "'normalize_input == TRUE' is only competable with ",
                "'input_scale == linear'. 'normalize_input' was ignored."
            )
        }
    }

    return(list_out)
}
