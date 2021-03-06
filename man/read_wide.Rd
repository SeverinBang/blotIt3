% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_exported_functions.R
\name{read_wide}
\alias{read_wide}
\title{Import measurement data}
\usage{
read_wide(file, description = NULL, time = 1, header = TRUE, ...)
}
\arguments{
\item{file}{character, the name of the data file.}

\item{description}{numeric index vector of the columns containing
the description.}

\item{time}{numeric index of length 1: the time column.}

\item{header}{a logical value indicating whether the file contains
the names of the variables as its first line.}

\item{...}{further arguments being passed to \link{read.csv}.}
}
\value{
data frame with columns "name", "time", "value" and other
columns describing the measurements.
}
\description{
Data from a given wide style csv-file is imported. While importing, the data
is converted into a long table
}
\examples{
## Import example data set
sim_data_wide_file <- system.file(
    "extdata", "sim_data_wide.csv",
    package = "blotIt3"
)
read_wide(sim_data_wide_file, description = seq_len(3))
}
