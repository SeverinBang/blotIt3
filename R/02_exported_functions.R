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



# plot_align_me() ---------------------------------------------------------

#' All-in-one plot function for blotIt3
#'
#' Takes the output of \link{align_me} and generates graphs. Which data will be
#' plotted can then be specified separately.
#'
#' @param out_list \code{out_list} file as produced by \link{align_me}, a list
#' of data.frames.
#' @param plot_points String to specify which data set should be plotted in form
#' of points with corresponding error bars. It must be one of
#' \code{c("original", "scaled", "prediction", "aligned")}.
#' @param plot_line Same as above but with a line and error band.
#' @param spline Logical, if set to \code{TRUE}, what is specified as
#' \code{plot_line} will be plotted as a smooth spline instead of straight lines
#' between points.
#' @param scales String passed as \code{scales} argument to \link{facet_wrap}.
#'
#' @param align_zeros Logical, if \code{TRUE}, the zero ticks are aligned
#' between the facets.
#' @param plot_caption Logical, if \code{TRUE}, a caption describing the plotted
#' data is added to the plot.
#' @param ncol Numerical passed as \code{ncol} argument to
#' \link{facet_wrap}.
#'
#' @param my_colors list of custom color values as taken by the \code{values}
#' argument in the \link{scale_color_manual} method for \code{ggplot} objects.
#'
#' @param duplicate_zero_points Logical, if set to \code{TRUE} all zero time
#' points are assumed to belong to the first condition. E.g. when the different
#' conditions consist of treatments added at time zero. Default is \code{FALSE}.
#'
#' @param my_order Optional list of target names in the custom order that will
#' be used for faceting.
#'
#' @param plot_scale_x character, defining the scale of the x axis
#'
#' @param plot_scale_y character, defining the scale of the y axis
#'
#' @param dose_response Logical, indicates if the plot should be dose response
#'
#' @param x_lab Optional, value passed to \code{xlab} parameter of \link{ggplot}
#' for the x-axis. Default is \code{NULL} leading to 'Time' or 'Dose',
#' respectively.
#'
#' @param y_lab Optional, value passed to \code{ylab} parameter of \link{ggplot}
#' for the y-axis. Default is \code{NULL} leading to 'Signal'.
#'
#' @param ... Logical expression used for subsetting the data frames, e.g.
#' \code{name == "pERK1" & time < 60}.
#'
#' \describe{
#' To reproduce the known function \code{plot1}, \code{plot2} and \code{plot3},
#' use:
#' \item{plot1}{
#' \code{plot_points} = 'original', \code{plot_line} = 'prediction'
#' }
#' \item{plot2}{
#' \code{plot_points} = 'scaled', \code{plot_line} = 'aligned'
#' }
#' \item{plot3}{
#' \code{plot_points} = 'aligned', \code{plot_line} = 'aligned'
#' }
#' }
#'
#' @import ggplot2 data.table
#'
#' @return ggplot object
#'
#' @export
#' @author Severin Bang and Svenja Kemmer

plot_align_me <- function(out_list,
                          ...,
                          plot_points = "aligned",
                          plot_line = "aligned",
                          spline = FALSE,
                          scales = "free",
                          align_zeros = TRUE,
                          plot_caption = TRUE,
                          ncol = NULL,
                          my_colors = NULL,
                          duplicate_zero_points = FALSE,
                          my_order = NULL,
                          plot_scale_y = NULL,
                          plot_scale_x = NULL,
                          dose_response = FALSE,
                          x_lab = NULL,
                          y_lab = NULL
                          ) {
    if (FALSE) {
        out_list <- out_list
        plot_points <- "aligned"
        plot_line <- "aligned"
        spline <- FALSE
        scales <- "free"
        align_zeros <- TRUE
        plot_caption <- TRUE
        ncol <- NULL
        my_colors <- NULL
        duplicate_zero_points <- FALSE
        my_order <- NULL
        plot_scale_y = NULL
        plot_scale_x = "log10"
        dose_response <- TRUE
        x_title <- "bingo"
        x_lab <- NULL
        y_lab <- NULL

    }
    if (!plot_points %in% c("original", "scaled", "prediction", "aligned") |
        !plot_line %in% c("original", "scaled", "prediction", "aligned")) {
        stop(
            "\n\t'plot_points' and 'plot_line' must each be one of
            c('original', 'scaled', 'prediction', 'aligned')\n"
        )
    }

    # if (class(out)[1] == "list") {
    #     cat("Data is in form of list, continuing\n")
    #     out <- out
    # } else {
    #     cat("Data is not yet in form of a list, passing it thorugh\n
    #     \t'blotIt_out_to_list(out,use_factors = T)'\n")
    #     out <- blotIt_out_to_list(out, use_factors = T)
    # }

    # change plotting order from default
    if (!is.null(my_order)) {
        if (length(setdiff(levels(out_list[[1]]$name), my_order)) != 0) {
            stop("my_order doesn't contain all protein names.")
        } else {
            out_list$aligned$name <- factor(
                out_list$aligned$name,
                levels = my_order
            )
            out_list$scaled$name <- factor(
                out_list$scaled$name,
                levels = my_order
            )
            out_list$prediction$name <- factor(
                out_list$prediction$name,
                levels = my_order
            )
            out_list$original$name <- factor(
                out_list$original$name,
                levels = my_order
            )
        }
    }


    distinguish <- out_list$distinguish
    scaling <- out_list$scaling

    plot_list <- out_list

    # duplicate 0 values for all doses
    if (duplicate_zero_points) {
        for (ndat in 1) {
            dat <- plot_list[[ndat]]
            subset_zeros <- copy(subset(dat, time == 0))
            mydoses <- setdiff(unique(dat$dose), 0)
            my_zeros_add <- NULL
            for (d in seq(1, length(mydoses))) {
                subset_zeros_d <- copy(subset_zeros)
                subset_zeros_d$dose <- mydoses[d]
                my_zeros_add <- rbind(my_zeros_add, subset_zeros_d)
            }
            dat <- rbind(dat, my_zeros_add)
            plot_list[[ndat]] <- dat
        }
    }

    # add columns containing the respective scaling and distinguish effects

    # aligned

    if (dose_response == TRUE) {
        x_value <- "dose"
    } else {
        x_value <- "time"
    }


    plot_list$aligned$distinguish <- do.call(
        paste0,
        plot_list$aligned[
            ,
            distinguish[!(distinguish %in% c("name", x_value))],
            drop = FALSE
        ]
    )
    plot_list$aligned$scaling <- NA

    # scaled
    plot_list$scaled$distinguish <- do.call(
        paste0,
        out_list$scaled[
            ,
            distinguish[!(distinguish %in% c("name", x_value))],
            drop = FALSE
        ]
    )
    plot_list$scaled$scaling <- do.call(
        paste0,
        out_list$scaled[, scaling[scaling != "name"], drop = FALSE]
    )

    # prediction
    plot_list$prediction$distinguish <- do.call(
        paste0,
        out_list$prediction[
            ,
            distinguish[!(distinguish %in% c("name", x_value))],
            drop = FALSE
        ]
    )
    plot_list$prediction$scaling <- do.call(
        paste0,
        out_list$prediction[, scaling[scaling != "name"], drop = FALSE]
    )

    # original
    plot_list$original$distinguish <- do.call(
        paste0,
        out_list$original[
            ,
            distinguish[!(distinguish %in% c("name", x_value))],
            drop = FALSE
        ]
    )
    plot_list$original$scaling <- do.call(
        paste0,
        out_list$original[, scaling[scaling != "name"], drop = FALSE]
    )

    plot_list_points <- plot_list[[plot_points]]
    plot_list_line <- plot_list[[plot_line]]

    plot_list_points <- subset(plot_list_points, ...)
    plot_list_line <- subset(plot_list_line, ...)

    # legend_name <- paste(scaling, collapse = ", ")

    # build Caption
    used_errors <- list(
        aligned = "Fisher Information",
        scaled = "Propergated error model to common scale",
        prediction = "Error model",
        original = "None"
    )

    used_data <- list(
        aligned = "Estimated true values",
        scaled = "Original data scaled to common scale",
        prediction = "Predictions from model evaluation on original scale",
        original = "Original data"
    )

    caption_text <- paste0(
        "Datapoints: ", used_data[[plot_points]], "\n",
        "Errorbars: ", used_errors[[plot_points]], "\n",
        "Line: ", used_data[[plot_line]], "\n",
        if (plot_points != plot_line) {
            paste0("Errorband: ", used_errors[[plot_line]], "\n")
        },
        "\n",
        "Date: ", Sys.Date()
    )

    # we want to keep the x ticks!
    if (scales == "distinguish") {
        scales <- "free_x"
    }



    # if (dose_response == FALSE ) {
    #     g <- plot_time_course(
    #         plot_list_points = plot_list_points,
    #         plot_list_line = plot_list_line,
    #         plot_points = plot_points,
    #         plot_line = plot_line,
    #         spline = spline,
    #         scales = scales,
    #         align_zeros = align_zeros,
    #         ncol = ncol,
    #         my_colors = my_colors,
    #         xlab = xlab,
    #         ylab = ylab
    #     )
    # } else {
    #     g <- plot_dose_response(
    #         plot_list_points = plot_list_points,
    #         plot_list_line = plot_list_line,
    #         plot_points = plot_points,
    #         plot_line = plot_line,
    #         spline = spline,
    #         scales = scales,
    #         align_zeros = align_zeros,
    #         ncol = ncol,
    #         my_colors = my_colors,
    #         xlab = xlab,
    #         ylab = ylab
    #     )
    # }



# * settings for dose response/time course --------------------------------

    if (dose_response == TRUE) {
        x_label <-  "Dose"
        x_variable <- "dose"
    } else {
        x_label <-  "Time"
        x_variable <- "time"
    }

    y_label <- "Signal"

    if (!is.null(x_lab)){
        x_label <- x_lab
    }

    if (!is.null(y_lab)){
        y_label <- y_lab
    }

    if (!is.null(plot_scale_x)) {
        errwidth <- 0
    } else {
        errwidth <- max(plot_list_points[x_variable])/50
    }


# * plotting --------------------------------------------------------------
    if (plot_points == "aligned" & plot_line == "aligned") {
        g <- ggplot(
            data = plot_list_points,
            aes_string(
                x = x_variable,
                y = "value",
                group = "distinguish",
                color = "distinguish",
                fill = "distinguish"
            )
        )
        g <- g + facet_wrap(~name, scales = scales, ncol = ncol)
        if (is.null(my_colors)) {
            # my_colors <- scale_color_brewer()
            # # c(
            # # "#000000",
            # # "#C5000B",
            # # "#0084D1",
            # # "#579D1C",
            # # "#FF950E",
            # # "#4B1F6F",
            # # "#CC79A7",
            # # "#006400",
            # # "#F0E442",
            # # "#8B4513",
            # # rep("gray", 100)
            # # )
            # g <- g + scale_color_manual("Condition", values = my_colors) +
            #     scale_fill_manual("Condition", values = my_colors)
        } else {
            my_colors <- c(my_colors, rep("gray", 100))
            g <- g + scale_color_manual("Distinguished", values = my_colors) +
                scale_fill_manual("Distinguished", values = my_colors)
        }
    } else {
        g <- ggplot(
            data = plot_list_points,
            aes_string(
                x = x_variable,
                y = "value",
                group = "scaling",
                color = "scaling",
                fill = "scaling"
            )
        )
        g <- g + facet_wrap(
            ~ name * distinguish,
            scales = scales,
            ncol = ncol
        )
        if (is.null(my_colors)) {
            # my_colors <- scale_color_brewer()
            # # c(
            # # "#000000",
            # # "#C5000B",
            # # "#0084D1",
            # # "#579D1C",
            # # "#FF950E",
            # # "#4B1F6F",
            # # "#CC79A7",
            # # "#006400",
            # # "#F0E442",
            # # "#8B4513",
            # # rep("gray", 100)
            # # )
            # g <- g + scale_color_manual("Scaling", values = my_colors) +
            #     scale_fill_manual("Scaling", values = my_colors)
        } else {
            my_colors <- c(my_colors, rep("gray", 100))
            g <- g + scale_color_manual("Scaled", values = my_colors) +
                scale_fill_manual("Scaled", values = my_colors)
        }
    }


    if (spline) {
        g <- g + geom_errorbar(
            data = plot_list_line,
            aes(
                ymin = lower, # value - sigma,
                ymax = upper # value + sigma
            ),
            width = 0
        )
        g <- g + geom_smooth(
            data = plot_list_line,
            se = FALSE,
            method = "lm",
            formula = y ~ poly(x, 3)
        )
    } else {
        if (plot_points == plot_line | plot_line == "prediction") {
            g <- g + geom_line(data = plot_list_line, size = 1)
            if (plot_line == "prediction") {
                g <- g + geom_ribbon(
                    data = plot_list_line,
                    aes(
                        # probably this should be changed back
                        ymin = lower,
                        ymax = upper
                    ),
                    alpha = 0.1,
                    lty = 0
                )
            }
        } else {
            g <- g + geom_line(data = plot_list_line, size = 1)#, color = "grey")
            g <- g + geom_ribbon(
                data = plot_list_line,
                aes(
                    ymin = lower, # value - sigma,
                    ymax = upper, # value + sigma#,
                    # fill = "grey",
                    # color = "grey"
                ),
                alpha = 0.3,
                lty = 0
            )
        }
    }

    g <- g + geom_point(data = plot_list_points, size = 2.5)

    if(plot_line != "prediction"){
        g <- g + geom_errorbar(
            data = plot_list_points,
            aes(
                ymin = lower, # value - sigma,
                ymax = upper # value + sigma
            ),
            size = 0.5,
            width = errwidth,
            alpha = 0.5
        )
    } else {
        g <- g + geom_errorbar(
            data = plot_list_points,
            aes(
                ymin = lower, # ,
                ymax = upper #
            ),
            size = 0.5,
            width = errwidth,
            alpha = 0.5
        )
    }


    g <- g + theme_bw(base_size = 20) +
        theme(
            legend.position = "top",
            legend.key = element_blank(),
            strip.background = element_rect(color = NA, fill = NA),
            axis.line.x = element_line(size = 0.3, colour = "black"),
            axis.line.y = element_line(size = 0.3, colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
        )
    g <- g + xlab(paste0("\n",x_label)) + ylab(paste0(y_label,"\n"))

    if (align_zeros) {
        if (plot_points != "original") {
            # scale y-axes (let them start at same minimum determined by
            # smallest value-sigma and end at individual ymax)
            plot_list_points <- as.data.table(plot_list_points)
            blank_data <- plot_list_points[
                ,
                list(ymax = max(upper), ymin = min(lower)),
                by = c("name", "distinguish", "scaling")
            ]
            blank_data[, ":="(ymin = min(ymin))] # same minimum for all proteins
            blank_data[
                ,
                ":="(ymax = ymaximal(ymax)),
                by = c("name", "distinguish", "scaling")
            ] # protein specific maximum
            blank_data <- melt(
                blank_data,
                id.vars = c("name", "distinguish", "scaling"),
                measure.vars = c("ymax", "ymin"),
                value.name = "value"
            )
            blank_data[, ":="(x_variable = 0, variable = NULL)]
            setnames(blank_data, "x_variable", x_variable)
            g <- g + geom_blank(
                data = as.data.frame(blank_data),
                aes_string(x = x_variable, y = "value")
            )
        }
    }


    if (plot_caption) {
        g <- g + labs(caption = caption_text)
    }

    if (is.null(plot_scale_y)) {
        if (out_list$output_scale != "linear") {
            g <- g + coord_trans(y = out_list$output_scale)
        }
    } else if (plot_scale_y %in% c("log", "log2", "log10")) {
        g <- g + scale_y_continuous(trans = plot_scale_y)
    }

    if (is.null(plot_scale_x)) {
        if (out_list$output_scale != "linear") {
            g <- g + coord_trans(x = out_list$output_scale)
        }
    } else if (plot_scale_x %in% c("log", "log2", "log10")) {
        g <- g + scale_x_continuous(trans = plot_scale_x)
    }


    return(g)
}


# llr_test() --------------------------------------------------------------

#' Method for hypothesis testing
#'
#' Two outputs of \link{align_me} can be tested as nested hypothesis. This can
#' be used to test if e.g. buffer material influence can be neglected or a
#' specific measurement point is an outlier.
#'
#' @param H0 output of \link{align_me} obeying the null hypothesis. A special
#' case of \code{H1}.
#' @param H1 output of \link{align_me}, the general case
#'
#' @return list with the log-likelihood ratio, statistical information and the
#' numerical p-value calculated by the evaluating the chi-squared distribution
#' at the present log-likelihood ratio with the current degrees of freedom.
#'
#' @examples
#' ## load provided example data file
#' lrr_data_path <- system.file(
#'     "extdata", "example_llr_test.csv",
#'     package = "blotIt3"
#' )
#'
#' ## import data
#' llr_data <- read_wide(
#'     file = lrr_data_path,
#'     description = seq(1, 4),
#'     sep = ",",
#'     dec = "."
#' )
#'
#' ## generate H0: the buffer column is not named as a distinguish effect e.g.
#' ## not considered as a biological different condition
#' H0 <- align_me(
#'     data = llr_data3,
#'     model = "yi / sj",
#'     error_model = "value * sigmaR",
#'     distinguish = yi ~ name + time + stimmulus,
#'     scaling = sj ~ name + ID,
#'     error = sigmaR ~ name + 1,
#'     parameter_fit_scale_log = FALSE,
#'     normalize = TRUE,
#'     average_techn_rep = FALSE,
#'     verbose = FALSE,
#'     normalize_input = TRUE
#' )
#'
#' ## generate H1: here the buffer column is named in the distinguish parameter
#' ## therefore different entries are considered as biologically different
#' H1 <- align_me(
#'     data = llr_data3,
#'     model = "yi / sj",
#'     error_model = "value * sigmaR",
#'     distinguish = yi ~ name + time + stimmulus + buffer,
#'     scaling = sj ~ name + ID,
#'     error = sigmaR ~ name + 1,
#'     parameter_fit_scale_log = FALSE,
#'     normalize = TRUE,
#'     average_techn_rep = FALSE,
#'     verbose = FALSE,
#'     normalize_input = TRUE
#' )
#'
#' ## perform test
#' llr_test(H0, H1)
#' @export
llr_test <- function(H0, H1, check = TRUE) {
    distinguish0 <- union(
        H0$distinguish[!(H0$distinguish %in% c("name", "time"))],
        "1"
    )
    distinguish1 <- union(
        H1$distinguish[!(H1$distinguish %in% c("name", "time"))],
        "1"
    )

    scaling0 <- union(H0$scaling[H0$scaling != "name"], "1")
    scaling1 <- union(H1$scaling[H1$scaling != "name"], "1")

    error0 <- union(H0$error[H0$error != "name"], "1")
    error1 <- union(H1$scaling[H1$scaling != "name"], "1")

    if (check) {
        if (
            !all(distinguish0 %in% distinguish1) |
                !all(scaling0 %in% scaling1) | !all(error0 %in% error1)
        ) {
            stop("H0 is not a special case of H1.")
        }
    }


    value0 <- attr(H0$parameter, "value")
    value1 <- attr(H1$parameter, "value")
    df0 <- attr(H0$parameter, "df")
    df1 <- attr(H1$parameter, "df")

    list(
        llr = value0 - value1,
        statistic = paste0(
            "chisquare with ", df0 - df1, " degrees of freedom."
        ),
        p.value = pchisq(value0 - value1, df = df0 - df1, lower.tail = FALSE)
    )
}
