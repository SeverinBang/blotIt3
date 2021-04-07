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
read_wide <- function(file, description = 1:3, time = 1, header = TRUE,  ...) {
    mydata <- read.csv(file, header = header, ...)
    allnames <- colnames(mydata)

    # Translate characters to numeric if necessary
    if (is.character(description)) ndescription <- which(colnames(mydata) %in% description) else ndescription <- description
    if (is.character(time)) ntime <- which(colnames(mydata) == time) else ntime <- time

    # Check availability of description and time
    if (length(ndescription) < length(description)) {
        warning("Not all columns proposed by argument 'description' are available in file. \nTaking the available ones.")
    }

    if (length(ntime) == 0) {
        stop("File did not contain a time column as proposed by 'time' argument.")
    }

    # Distinguish description data from measurement data
    Description <- mydata[, ndescription]
    restLong <- unlist(mydata[, -ndescription])

    # Create output data frame
    newdata <- data.frame(
        Description,
        name = rep(allnames[-ndescription], each = dim(mydata)[1]),
        value = restLong
    )

    # Remove missing items
    newdata <- newdata[!is.nan(newdata$value), ]
    newdata <- newdata[!is.na(newdata$value), ]


    colnames(newdata)[ntime] <- "time"

    return(newdata)
}
