################################################################################
#  Survival data tranformation for fitting a Poisson model                     #
################################################################################
#                                                                              #
#  This function  tranform survival data into a format compatible with         #
#   the glm() function for fitting a Poisson model.                            #
#                                                                              #
#  Its parameters are                                                          #
#   - data          : a data frame with columns:                               #
#                     - id       : the patient identifyier                     #
#                     - time     : the event/censoring time                    #
#                     - status   : the event(1) or censoring(0) indicator      #
#                     - ...      : other factors such like the covariables     #
#                             needed in the regression model                   #
#   - all.breaks    : the breakpoints between time intervals                   #
#   - interval.width: the width of the time intervals on wich the risks        #
#                     will be assumed consant.                                 #
#                     To be used if intervals of the same length               #
#                     This parameter is ignored if all.breaks is specified     #
#   - nInts         : the number of intervals containing the same expected     #
#                     number of events                                         #
#                     This parameter is ignored if either all.breaks or        #
#                     interval.width is specified                              #
#   - factors       : a vector of characters, containing the names of the      #
#                     factors to be kept in the transformed data set           #
#   - compress      : a logical, indicating whether the record with the same   #
#                     factor profile should be summarized into one record,     #
#                     i.e. whether the data should be expressed in a short form#
#                                                                              #
################################################################################
#                                                                              #
# Author: Federico Rotolo <federico.rotolo@gustaveroussy.fr>                   #
# Original code by Stephanie Kovalchik                                         #
# http://r.789695.n4.nabble.com/exponential-proportional-hazard-model-td805536.html #
#                                                                              #
#   Date:                 June 22, 2015                                        #
#   Last modification on: June 15, 2016                                        #
#################################################################################
poissonize <- function(data,
                       all.breaks = NULL,
                       interval.width = NULL,
                       nInts = 8,
                       factors = NULL,
                       compress = TRUE) {
  if (is.null(all.breaks)) {
    if (!is.null(interval.width)) {
      all.breaks <- seq(0, max(data$time),
                        interval.width)
    } else {
      # library('survival')
      smod <- survfit(Surv(time, status) ~ 1, data = data)
      smin <- min(smod$surv)
      all.breaks <- c(0, sapply(1 - 1:nInts / nInts, function(p)
        smod$time[(smod$surv - smin) / (1 - smin) <= p][1]))
    }
  }
  tol <- max(data$time) / 1e4
  all.breaks <- sort(unique(setdiff(round(
    c(0, all.breaks, max(data$time)) / tol
  ),
  NA))) * tol
  
  # THE FOLLOWING FUNCTION COMPUTES ALL THE TIME BREAKS GIVEN THE EXIT TIME
  person.breaks <- Vectorize(function(stopT)
    unique(c(all.breaks[all.breaks < stopT], stopT) / tol) * tol, SIMPLIFY = FALSE)
  
  # NEXT WE GET A LIST OF EACH SUBJECT'S TIME PERIODS
  the.breaks <- person.breaks(data$time)
  
  # NOW WE OBTAIN THE EXPANDED PERIOD OF OBSERVATION
  startT <-
    lapply(the.breaks, function(x)
      x[-length(x)])  # LEFT TIME POINT
  stopT <-
    lapply(the.breaks, function(x)
      x[-1])    # RIGHT TIME POINTS
  
  # THE FOLLOWING ARE NEEDED TO COMPLETE THE LONG VERSION OF THE DATA SET
  count.per.id <- sapply(startT, length)
  index <- tapply(data$id, data$id, length)
  index <-
    cumsum(index) # INDEX OF LAST OBSERVATION FOR EACH PATIENT
  event <- rep(0, sum(count.per.id))
  event[cumsum(count.per.id)] <- data$status[index]
  
  # BRING ALL OF THIS TOGETHER TO CREATE THE EXPANDED DATASET
  expData <- cbind(data.frame(
    id = rep(data$id[index], count.per.id),
    startT = unlist(startT),
    stopT = unlist(stopT),
    event = event
  ),
  data[sapply(rep(data$id[index], count.per.id),
              function(x)
                which(data$id == x)), factors, drop = FALSE])
  
  
  
  # CREATE TIME VARIABLE WHICH INDICATES THE PERIOD OF OBSERVATION
  # THIS WILL BE THE OFFSET IN THE POISSON MODEL FIT
  expData$time <-
    expData$stopT - expData$startT # LENGTH OF OBSERVATION
  
  # NEXT WE CREATE A FACTOR FOR EACH INTERVAL THAT WILL ALLOW US
  # TO HAVE A DIFFERENT RATE FOR EACH PERIOD
  expData$interval <-
    factor(expData$start, levels = unique(expData$start))
  expData <- subset(expData, select = -c(startT, stopT))
  #expData$time <- expData$time
  if (compress) {
    m <- aggregate(x = expData$event,
                   by = expData[, c('interval', factors), drop = FALSE],
                   sum)
    colnames(m)[ncol(m)] <- 'm'
    Rt <- aggregate(x = expData$time,
                    by = expData[, c('interval', factors), drop = FALSE],
                    sum)
    colnames(Rt)[ncol(Rt)] <- 'Rt'
    N <- aggregate(x = expData$event,
                   by = expData[, c('interval', factors), drop = FALSE],
                   function(x)
                     sum(!is.na(x)))
    colnames(N)[ncol(N)] <- 'N'
    expData <- merge(m, merge(Rt, N, sort = FALSE), sort = FALSE)
  }
  
  attr(expData, 'interval.width') <- interval.width
  attr(expData, 'all.breaks') <- all.breaks
  return(expData)
}


################################################################################
#  Plots the hazard function and the survival curve                            #
#     estimated via a Poisson model                                            #
################################################################################
#                                                                              #
#  Parameters                                                                  #
#   - x          : the fitted Poisson model                                    #
#   - type       : the type of plot, either 'haz' for the hazard function      #
#                  or 'Surv', for the survival curve                           #
#                                                                              #
################################################################################
#                                                                              #
# Author: Federico Rotolo <federico.rotolo@gustaveroussy.fr>                   #
#                                                                              #
#   Date:                 June 22, 2015                                        #
#   Last modification on: November 14, 2016                                    #
################################################################################
plotsson <- function(x,
                     type = c('survival', 'hazard'),
                     add = FALSE,
                     xscale = 1,
                     by,
                     col,
                     ...) {
  type <- tolower(type)
  type <- match.arg(type)
  
  breaks <- attr(x$data, 'all.breaks')
  risks <- exp(coef(x)[grep('interval', names(coef(x)))]) *
    diff(breaks)
  
  if (type == 'survival') {
    risks <- exp(-cumsum(risks))
    Ylims <- 0:1
    Ylab <- 'Survival probability'
  } else {
    Ylims <- c(0, 1.1 * max(risks))
    Ylab <- 'Hazard'
  }
  
  startT <- breaks[-length(breaks)] / xscale
  stopT <- breaks[-1] / xscale
  
  if (!add)
    plot(
      0,
      xlim = c(0, 1.1 * max(stopT)),
      ylim = Ylims,
      xlab = 'time',
      ylab = Ylab,
      type = 'n',
      ...
    )
  
  byVals <- if (missing(by)) {
    0
  } else {
    unique(x$data[, by])
  }
  beta <- if (missing(by)) {
    0
  } else {
    coef(x)[by]
  }
  
  if (type == 'survival') {
    for (v in 1:length(byVals))
      segments(
        startT,
        c(1, risks[2:length(risks) - 1] ^ exp(beta * byVals[v])),
        stopT,
        risks ^ exp(beta * byVals[v]),
        col = ifelse(missing(col), v, col),
        ...
      )
  } else {
    for (v in 1:length(byVals))
      segments(
        startT,
        risks * exp(beta * byVals[v]),
        stopT,
        risks * exp(beta * byVals[v]),
        col = ifelse(missing(col), v, col),
        ...
      )
  }
  
}
