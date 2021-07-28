myprofile <- function (obj, pars, whichPar, alpha = 0.05, limits = c(lower = -Inf,
                                                        upper = Inf), method = c("integrate", "optimize"), stepControl = NULL,
          algoControl = NULL, optControl = NULL, verbose = FALSE,
          cores = 1, ...)
{
  if (FALSE) {
    obj <- objective_function
    pars <- fit_result$argument
    whichPar <- 1
    dotArgs <- list(
      pass_parameter_list=pass_parameter_list,
      pass_parameter_list2=pass_parameter_list2
    )
    alpha <- 0.05
    limits <- c(
      lower = -Inf,
      upper = Inf)
    method <- "optimize"
    stepControl <- NULL
    algoControl <- NULL
    optControl <- NULL
    verbose <- FALSE
    cores <- 1
  }
  dotArgs <- list(...)
  #sanePars <- sanitizePars(pars, dotArgs$fixed)
  sanePars <- list(pars = pars, fixed = dotArgs$fixed)
  pars <- sanePars$pars
  fixed <- sanePars$fixed
  method <- match.arg(method)
  if (method == "integrate") {
    sControl <- list(stepsize = 1e-04, min = 1e-04, max = Inf,
                     atol = 0.01, rtol = 0.01, limit = 500, stop = "value")
    aControl <- list(gamma = 1, W = "hessian", reoptimize = FALSE,
                     correction = 1, reg = .Machine$double.eps)
    oControl <- list(rinit = 0.1, rmax = 10, iterlim = 10,
                     fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  if (method == "optimize") {
    sControl <- list(stepsize = 0.01, min = 1e-04, max = Inf,
                     atol = 0.1, rtol = 0.1, limit = 100, stop = "value")
    aControl <- list(gamma = 0, W = "identity", reoptimize = TRUE,
                     correction = 1, reg = 0)
    oControl <- list(rinit = 0.1, rmax = 10, iterlim = 100,
                     fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  cores <- min(length(whichPar), cores)
  cores <- sanitizeCores(cores)
  if (!is.null(stepControl))
    sControl[match(names(stepControl), names(sControl))] <- stepControl
  if (!is.null(algoControl))
    aControl[match(names(algoControl), names(aControl))] <- algoControl
  if (!is.null(optControl))
    oControl[match(names(optControl), names(oControl))] <- optControl
  if (cores > 1) {
    if (Sys.info()[["sysname"]] == "Windows") {
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl = cluster)
      parallel::clusterCall(cl = cluster, function(x) .libPaths(x),
                            .libPaths())
      varlist <- ls()
      varlist <- c("obj", "whichPar", "alpha", "limits",
                   "method", "verbose", "cores", "pars", "fixed",
                   "dotArgs", "sControl", "aControl", "oControl")
      parallel::clusterExport(cluster, envir = environment(),
                              varlist = varlist)
    }
    else {
      doParallel::registerDoParallel(cores = cores)
    }
    "%mydo%" <- foreach::"%dopar%"
  }
  else {
    "%mydo%" <- foreach::"%do%"
  }
  if (is.character(whichPar))
    whichPar <- which(names(pars) %in% whichPar)
  loaded_packages <- .packages()
  out <- foreach::foreach(whichIndex = whichPar, .packages = loaded_packages,
                          .inorder = TRUE, .options.multicore = list(preschedule = FALSE)) %mydo%
    {
      #loadDLL(obj)
      whichPar.name <- names(pars)[whichIndex]
      obj.opt <- obj
      obj.prof <- function(p, ...) {
        out <- obj(p, ...)
        Id <- diag(1/.Machine$double.eps, length(out$gradient))
        Id[whichIndex, whichIndex] <- 1
        colnames(Id) <- rownames(Id) <- names(out$gradient)
        W <- match.arg(aControl$W[1], c("hessian", "identity"))
        out$hessian <- switch(W, hessian = out$hessian,
                              identity = Id)
        return(out)
      }
      pseudoinverse <- function(m, tol) {
        msvd <- svd(m)
        index <- which(abs(msvd$d) > max(dim(m)) * max(msvd$d) *
                         tol)
        if (length(index) == 0) {
          out <- array(0, dim(m)[2:1])
        }
        else {
          out <- msvd$u[, index] %*% (1/msvd$d[index] *
                                        t(msvd$v)[index, ])
        }
        attr(out, "valid") <- 1:length(msvd$d) %in%
          index
        return(out)
      }
      constraint <- function(p) {
        value <- p[whichIndex] - pars[whichIndex]
        gradient <- rep(0, length(p))
        gradient[whichIndex] <- 1
        return(list(value = value, gradient = gradient))
      }
      lagrange <- function(y) {
        p <- y
        lambda <- 0
        out <- do.call(obj.prof, c(list(p = p), dotArgs))
        g.original <- constraint(p)
        g <- direction * g.original$value
        gdot <- direction * g.original$gradient
        ldot <- out$gradient
        lddot <- out$hessian
        M <- rbind(cbind(lddot, gdot), matrix(c(gdot,
                                                0), nrow = 1))
        v <- c(-rep(gamma, length(p)) * (ldot + lambda *
                                           gdot), 1)
        v0 <- c(-rep(0, length(p)) * (ldot + lambda *
                                        gdot), 1)
        W <- pseudoinverse(M, tol = aControl$reg)
        valid <- attr(W, "valid")
        if (any(!valid)) {
          dy <- try(as.vector(W %*% v)[1:length(p)],
                    silent = FALSE)
          dy0 <- try(as.vector(W %*% v0)[1:length(p)],
                     silent = FALSE)
          dy[!valid[1:length(p)]] <- dy0[!valid[1:length(p)]] <- 0
          dy[whichIndex] <- dy0[whichIndex] <- direction
          warning(paste0("Iteration ", i, ": Some singular values of the Hessian are below the threshold. Optimization will be performed."))
        }
        else {
          dy <- try(as.vector(W %*% v)[1:length(p)],
                    silent = FALSE)
          dy0 <- try(as.vector(W %*% v0)[1:length(p)],
                     silent = FALSE)
        }
        if (!inherits(dy, "try-error")) {
          names(dy) <- names(y)
          correction <- sqrt(sum((dy - dy0)^2))/sqrt(sum(dy^2))
        }
        else {
          dy <- NA
          correction <- 0
          warning(paste0("Iteration ", i, ": Impossible to invert Hessian. Trying to optimize instead."))
        }
        out.attributes <- attributes(out)[sapply(attributes(out),
                                                 is.numeric)]
        out.attributes.names <- names(out.attributes)
        return(c(list(dy = dy, value = out$value, gradient = out$gradient,
                      correction = correction, valid = valid, attributes = out.attributes.names),
                 out.attributes))
      }
      doIteration <- function() {
        optimize <- aControl$reoptimize
        if (is.na(dy[1])) {
          optimize <- TRUE
          y.try <- y
          y.try[whichIndex] <- y[whichIndex] + direction *
            stepsize
          rinit <- oControl$rinit
        }
        else {
          dy.norm <- sqrt(sum(dy^2))
          rinit <- min(c(oControl$rinit, 3 * dy.norm))
          y.try <- y + dy
          if (any(!lagrange.out$valid))
            optimize <- TRUE
        }
        if (optimize) {
          parinit.opt <- y.try[-whichIndex]
          fixed.opt <- c(fixed, y.try[whichIndex])
          arglist <- c(list(objfun = obj.opt, parinit = parinit.opt,
                            fixed = fixed.opt, rinit = rinit), oControl[names(oControl) !=
                                                                          "rinit"], dotArgs[names(dotArgs) != "fixed"])
          myfit <- try(do.call(trust::trust, arglist), silent = FALSE)
          if (!inherits(myfit, "try-error")) {
            y.try[names(myfit$argument)] <- as.vector(myfit$argument)
          }
          else {
            warning("Optimization not successful. Profile may be erroneous.")
          }
        }
        return(y.try)
      }
      doAdaption <- function() {
        lagrange.out.try <- lagrange(y.try)
        valid <- TRUE
        dobj.pred <- sum(lagrange.out$gradient * (y.try -
                                                    y))
        dobj.fact <- lagrange.out.try$value - lagrange.out$value
        correction <- lagrange.out.try$correction
        if (correction > aControl$correction)
          gamma <- gamma/2
        if (correction < 0.5 * aControl$correction)
          gamma <- min(c(aControl$gamma, gamma * 2))
        if (abs(dobj.fact - dobj.pred) > sControl$atol &
            stepsize > sControl$min) {
          stepsize <- max(c(stepsize/1.5, sControl$min))
          valid <- FALSE
        }
        if (abs(dobj.fact - dobj.pred) < 0.3 * sControl$atol |
            abs((dobj.fact - dobj.pred)/dobj.fact) < 0.3 *
            sControl$rtol) {
          stepsize <- min(c(stepsize * 2, sControl$max))
        }
        if (verbose) {
          diff.thres <- diff.steps <- diff.limit <- 0
          if (threshold < Inf)
            diff.thres <- 1 - max(c(0, min(c(1, (threshold -
                                                   lagrange.out.try$value)/delta))))
          if (sControl$limit < Inf)
            diff.steps <- i/sControl$limit
          diff.limit <- switch(as.character(sign(constraint.out$value)),
                               `1` = 1 - (limits[2] - constraint.out$value)/limits[2],
                               `-1` = diff.limit <- 1 - (limits[1] - constraint.out$value)/limits[1],
                               `0` = 0)
          percentage <- max(c(diff.thres, diff.steps,
                              diff.limit), na.rm = TRUE) * 100
          progressBar(percentage)
          myvalue <- format(substr(lagrange.out$value,
                                   0, 8), width = 8)
          myconst <- format(substr(constraint.out$value,
                                   0, 8), width = 8)
          mygamma <- format(substr(gamma, 0, 8), width = 8)
          myvalid <- all(lagrange.out$valid)
          cat("\tvalue:", myvalue, "constraint:", myconst,
              "gamma:", mygamma, "valid:", myvalid)
        }
        return(list(lagrange = lagrange.out.try, stepsize = stepsize,
                    gamma = gamma, valid = valid))
      }
      i <- 0
      direction <- 1
      gamma <- aControl$gamma
      stepsize <- sControl$stepsize
      ini <- pars
      lagrange.out <- lagrange(ini)
      constraint.out <- constraint(pars)
      delta <- qchisq(1 - alpha, 1)
      threshold <- lagrange.out[[sControl$stop]] + delta
      out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
      out <- c(value = lagrange.out$value, constraint = as.vector(constraint.out$value),
               stepsize = stepsize, gamma = gamma, whichPar = whichIndex,
               out.attributes, ini)
      if (verbose) {
        cat("Compute right profile\n")
      }
      direction <- 1
      gamma <- aControl$gamma
      stepsize <- sControl$stepsize
      y <- ini
      lagrange.out <- lagrange.out
      constraint.out <- constraint.out
      while (i < sControl$limit) {
        sufficient <- FALSE
        retry <- 0
        while (!sufficient & retry < 5) {
          dy <- stepsize * lagrange.out$dy
          y.try <- try(doIteration(), silent = TRUE)
          out.try <- try(doAdaption(), silent = TRUE)
          if (inherits(y.try, "try-error") | inherits(out.try,
                                                      "try-error")) {
            sufficient <- FALSE
            stepsize <- stepsize/1.5
            retry <- retry + 1
          }
          else {
            sufficient <- out.try$valid
            stepsize <- out.try$stepsize
          }
        }
        if (inherits(y.try, "try-error") | inherits(out.try,
                                                    "try-error"))
          break
        y <- y.try
        lagrange.out <- out.try$lagrange
        constraint.out <- constraint(y.try)
        stepsize <- out.try$stepsize
        gamma <- out.try$gamma
        out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
        out <- rbind(out, c(value = lagrange.out$value,
                            constraint = as.vector(constraint.out$value),
                            stepsize = stepsize, gamma = gamma, whichPar = whichIndex,
                            out.attributes, y))
        value <- lagrange.out[[sControl$stop]]
        if (value > threshold | constraint.out$value >
            limits[2])
          break
        i <- i + 1
      }
      if (verbose) {
        cat("\nCompute left profile\n")
      }
      i <- 0
      direction <- -1
      gamma <- aControl$gamma
      stepsize <- sControl$stepsize
      y <- ini
      lagrange.out <- lagrange(ini)
      constraint.out <- constraint(pars)
      while (i < sControl$limit) {
        sufficient <- FALSE
        retry <- 0
        while (!sufficient & retry < 5) {
          dy <- stepsize * lagrange.out$dy
          y.try <- try(doIteration(), silent = TRUE)
          out.try <- try(doAdaption(), silent = TRUE)
          if (inherits(y.try, "try-error") | inherits(out.try,
                                                      "try-error")) {
            sufficient <- FALSE
            stepsize <- stepsize/1.5
            retry <- retry + 1
          }
          else {
            sufficient <- out.try$valid
            stepsize <- out.try$stepsize
          }
        }
        if (inherits(y.try, "try-error") | inherits(out.try,
                                                    "try-error"))
          break
        y <- y.try
        lagrange.out <- out.try$lagrange
        constraint.out <- constraint(y.try)
        stepsize <- out.try$stepsize
        gamma <- out.try$gamma
        out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
        out <- rbind(c(value = lagrange.out$value, constraint = as.vector(constraint.out$value),
                       stepsize = stepsize, gamma = gamma, whichPar = whichIndex,
                       out.attributes, y), out)
        value <- lagrange.out[[sControl$stop]]
        if (value > threshold | constraint.out$value <
            limits[1])
          break
        i <- i + 1
      }
      out <- as.data.frame(out)
      out$whichPar <- paste0(whichPar.name) #, out$whichPar
      parframe(out, parameters = names(pars), metanames = c("value",
                                                            "constraint", "stepsize", "gamma", "whichPar"),
               obj.attributes = names(out.attributes))
    }
  if (Sys.info()[["sysname"]] == "Windows" & cores > 1) {
    parallel::stopCluster(cluster)
    doParallel::stopImplicitCluster()
  }
  do.call(rbind, out)
}

# adjusted profiles function ----------------------------------------------




# helpers -----------------------------------------------------------------


sanitizePars <- function(pars = NULL, fixed = NULL) {

  # Convert fixed to named numeric
  if (!is.null(fixed)) fixed <- structure(as.numeric(fixed), names = names(fixed))

  # Convert pars to named numeric
  if (!is.null(pars)) {
    pars <- structure(as.numeric(pars), names = names(pars))
    # remove fixed from pars
    pars <- pars[setdiff(names(pars), names(fixed))]
  }


  return(list(pars = pars, fixed = fixed))

}


sanitizeCores <- function(cores)  {

  max.cores <- parallel::detectCores()
  min(max.cores, cores)
  #
  # if (Sys.info()[['sysname']] == "Windows") cores <- 1
  # return(cores)

}

sanitizeConditions <- function(conditions) {

  new <- str_replace_all(conditions, "[[:punct:]]", "_")
  new <- str_replace_all(new, "\\s+", "_")
  return(new)

}


sanitizeData <- function(x, required = c("name", "time", "value"), imputed = c(sigma = NA, lloq = -Inf)) {

  all.names <- names(x)

  missing.required <- setdiff(required, all.names)
  missing.imputed <- setdiff(names(imputed), all.names)

  if (length(missing.required) > 0)
    stop("These mandatory columns are missing: ", paste(missing.required, collapse = ", "))

  if (length(missing.imputed) > 0) {

    for (n in missing.imputed) x[[n]] <- imputed[n]

  }

  list(data = x, columns = c(required, names(imputed)))

}

#' Generate a parameter frame
#'
#' @description A parameter frame is a data.frame where the rows correspond to different
#' parameter specifications. The columns are divided into three parts. (1) the meta-information
#' columns (e.g. index, value, constraint, etc.), (2) the attributes of an objective function
#' (e.g. data contribution and prior contribution) and (3) the parameters.
#' @seealso \link{profile}, \link{mstrust}
#' @param x data.frame.
#' @param parameters character vector, the names of the parameter columns.
#' @param metanames character vector, the names of the meta-information columns.
#' @param obj.attributes character vector, the names of the objective function attributes.
#' @return An object of class \code{parframe}, i.e. a data.frame with attributes for the
#' different names. Inherits from data.frame.
#' @details Parameter frames can be subsetted either by \code{[ , ]} or by \code{subset}. If
#' \code{[ , index]} is used, the names of the removed columns will also be removed from
#' the corresponding attributes, i.e. metanames, obj.attributes and parameters.
#' @example inst/examples/parlist.R
#' @export
parframe <- function(x = NULL, parameters = colnames(x), metanames = NULL, obj.attributes = NULL) {

  if (!is.null(x)) {
    rownames(x) <- NULL
    out <- as.data.frame(x)
  } else {
    out <- data.frame()
  }

  attr(out, "parameters") <- parameters
  attr(out, "metanames") <- metanames
  attr(out, "obj.attributes") <- obj.attributes
  class(out) <- c("parframe", "data.frame")

  return(out)

}



plotProfile.parframe <- function(profs, ..., maxvalue = 5, parlist = NULL) {

  if("parframe" %in% class(profs))
    arglist <- list(profs)
  else
    arglist <- as.list(profs)


  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }

  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")

    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }

    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])


    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }

    subdata <- do.call(rbind, lapply(names(proflist), function(n) {

      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue

      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)

      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }

      return(sub)
    }))
    return(subdata)
  }))

  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))

  data.zero <- subset(data, is.zero)

  threshold <- c(1, 2.7, 3.84)

  data <- droplevels.data.frame(subset(data, ...))

  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") +
    geom_hline(yintercept=threshold, lty=2, color="gray") +
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")

  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)

    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)

}


#' @export
#' @rdname plotProfile
plotProfile.list <- function(profs, ..., maxvalue = 5, parlist = NULL) {

  if("parframe" %in% class(profs))
    arglist <- list(profs)
  else
    arglist <- as.list(profs)


  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }

  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")

    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }

    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])


    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }

    subdata <- do.call(rbind, lapply(names(proflist), function(n) {

      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue

      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)

      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }

      return(sub)
    }))
    return(subdata)
  }))

  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))

  data.zero <- subset(data, is.zero)

  threshold <- c(1, 2.7, 3.84)

  data <- droplevels.data.frame(subset(data, ...))

  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") +
    geom_hline(yintercept=threshold, lty=2, color="gray") +
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")

  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)

    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)

}

#' Profile likelihood plot
#'
#' @param profs Lists of profiles as being returned by \link{profile}.
#' @param ... logical going to subset before plotting.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @return A plot object of class \code{ggplot}.
#' @details See \link{profile} for examples.
#' @export
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}



#' Profile uncertainty extraction
#'
#' @description extract parameter uncertainties from profiles
#' @param object object of class \code{parframe}, returned from \link{profile} function.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... not used right now.
#' @param val.column the value column used in the parframe, usually 'data'.
#' @export
confint.parframe <- function(object, parm = NULL, level = 0.95, ..., val.column = "data") {

  profile <- object
  obj.attributes <- attr(profile, "obj.attributes")

  if (is.null(parm))
    parm <- unique(profile[["whichPar"]])

  threshold <- qchisq(level, df = 1)

  # Reduce to profiles for parm
  profile <- profile[profile[["whichPar"]] %in% parm,]

  # Evaluate confidence intervals per parameter
  CIs <- lapply(split(profile, profile[["whichPar"]]), function(d) {

    # Get origin of profile
    origin <- which.min(abs(d[["constraint"]]))
    whichPar <- d[["whichPar"]][origin]

    # Define function to return constraint value where threshold is passed
    get_xThreshold <- function(branch) {

      y <- branch[[val.column]] - d[[val.column]][origin]
      x <- branch[["constraint"]]

      # If less than 3 points, return NA
      if (length(x) < 3)
        return(NA)

      # If threshold exceeded, take closest points below and above threshold
      # and interpolate
      if (any(y > threshold)) {
        i.above <- utils::head(which(y > threshold), 1)
        i.below <- utils::tail(which(y < threshold), 1)
        if (i.below > i.above) {
          return(NA)
        } else {
          slope <- (y[i.above] - y[i.below])/(x[i.above] - x[i.below])
          dy <- threshold - y[i.below]
          dx <- dy/slope
          x_threshold <- x[i.below] + dx
          return(x_threshold)
        }
      }

      # If threshold not exceeded,
      # take the last 20% of points (at least 3) an perform linear fit
      n_last20 <- max(3, length(which(x - x[1] > 0.8 * (max(x) - x[1]))))
      x <- tail(x, n_last20)
      y <- tail(y, n_last20)
      slope <- sum((x - mean(x))*(y - mean(y)))/sum((x - mean(x))^2)

      # If slope < 0, return Inf
      if (slope < 0)
        return(Inf)

      # Extrapolate until threshold is passed
      dy <- threshold - tail(y, 1)
      dx <- dy/slope
      x_threshold <- tail(x, 1) + dx

      # Test if extrapolation takes the point of passage very far
      # Set to Inf in that case
      if (x_threshold > 10*(max(x) - min(x)))
        x_threshold <- Inf

      return(x_threshold)


    }

    # Right profiles
    right <- d[d[["constraint"]] >= 0,]
    upper <- d[[whichPar]][origin] + get_xThreshold(right)

    # Left profile
    left <- d[d[["constraint"]] <= 0,]
    left[["constraint"]] <- - left[["constraint"]]
    left <- left[order(left[["constraint"]]), ]
    lower <- d[[whichPar]][origin] - get_xThreshold(left)


    data.frame(name = whichPar,
               value = d[[whichPar]][origin],
               lower = lower,
               upper = upper)


  })

  do.call(rbind, CIs)

}
