# estimateQ.WGD -----------------------------------------------------------

estimateQ.WGD = function (x, m, alleleSet, alleleFreq = NULL, history, init = NULL, 
                          maxiter = 100, tol = 1e-04, xGreaterZero = TRUE, useGradient = TRUE) 
{
  N <- length(x)
  possAlleles <- alleleSet
  lengthOfTheta <- length(possAlleles) - 1
  if (nrow(history) != length(possAlleles)) 
    stop("length of alleleSet must equal the number of rows of history.")
  Avec <- as.vector(history)
  
  easyType <- "no"
  if (identical(Avec, c(0, 1, 2, 0))) 
    easyType <- "CNLOH"
  if (identical(Avec, c(1, 1, 3, 0))) 
    easyType <- "SingleGain"
  if (identical(Avec, c(0,2,4,0)))
    easyType = "WGD"
  
  
  if (easyType == "CNLOH") 
    Ainv = matrix(c(0, 1, 0.5, 0), byrow = TRUE, ncol = 2)
  if (easyType == "SingleGain") 
    Ainv = matrix(c(0, 1, 1/3, -1/3), byrow = TRUE, ncol = 2)
  if (easyType=="WGD")
    Ainv = matrix(c(0, 0.5, 0.25, 0), ncol=2, nrow=2, byrow=TRUE)
  if (easyType == "no") 
    Ainv = matrix(c(0, 0, 1, 0, 1, 0, 0.25, -0.5, -0.25), 
                  byrow = TRUE, ncol = 3)
  
  if (is.null(init)) {
    init <- history %*% rep(1/length(possAlleles), length = length(possAlleles))
    init <- as.vector(init)/sum(init)
  }
  EStep <- function(q) {
    logPPerAllele <- mapply(possAlleles, q, FUN = function(a, 
                                                           qa) {
      dbinom(x = x, size = m, prob = a, log = TRUE) + log(qa)
    })
    rowMax <- apply(logPPerAllele, 1, max)
    logPPerAlleleMinusMax <- sweep(logPPerAllele, 1, rowMax, 
                                   "-")
    logdenom <- rowMax + log(apply(exp(logPPerAlleleMinusMax), 
                                   1, sum))
    logPPerAllele <- sweep(logPPerAllele, 1, logdenom, "-")
    return(exp(logPPerAllele))
  }
  if (easyType == "no" | easyType == "WGD") {
    uiL1 <- diag(rep(-1, lengthOfTheta), nrow = lengthOfTheta, 
                 ncol = lengthOfTheta)
    ciL1 <- rep(1, lengthOfTheta)
    uiGr0 <- diag(rep(1, lengthOfTheta), nrow = lengthOfTheta, 
                  ncol = lengthOfTheta)
    ciGr0 <- rep(0, lengthOfTheta)
    uiA <- Ainv %*% (rbind(diag(lengthOfTheta), rep(-1, lengthOfTheta)))
    ciA <- -Ainv %*% c(rep(0, lengthOfTheta), 1)
    ui <- rbind(uiL1, uiGr0, uiA)
    ci <- c(-ciL1, -ciGr0, ciA)
  }
  MStep <- function(PMat = NULL, qinit) {
    if (is.null(alleleFreq)) 
      Na <- colSums(PMat)
    else Na <- alleleFreq
    if (easyType %in% c("CNLOH", "SingleGain")) {
      if (easyType == "CNLOH") {
        a <- 0
        b <- 1
      }
      if (easyType == "SingleGain") {
        a <- 1/2
        b <- 1
      }
      qout <- Na[1]/sum(Na)
      if (qout > b) 
        qout <- b
      if (qout < a) 
        qout <- a
      qout <- c(qout, 1 - qout)
    }
    else {
      f <- function(q) {
        qS <- 1 - sum(q)
        ll <- sum(Na * log(c(q, qS)))
        if (xGreaterZero) {
          missZero <- 1 - rowSums(mapply(possAlleles, 
                                         c(q, qS), FUN = function(a, qa) {
                                           (1 - a)^m * qa
                                         }))
          ll <- ll - sum(missZero)
        }
        return(-ll)
      }
      aS <- tail(possAlleles, 1)
      num <- sapply(head(possAlleles, -1), FUN = function(a, 
                                                          qa) {
        (1 - aS)^m - (1 - a)^m
      })
      gr <- function(q) {
        qS <- 1 - sum(q)
        NS <- tail(Na, 1)
        g <- (head(Na, -1)/q - NS/qS)
        if (xGreaterZero) {
          missZero <- 1 - rowSums(mapply(possAlleles, 
                                         c(q, qS), FUN = function(a, qa) {
                                           (1 - a)^m * qa
                                         }))
          g <- g - colSums(sweep(num, 1, missZero, "/"))
        }
        return(-g)
      }
      theta <- head(qinit, -1)
      if (lengthOfTheta == 1) {
        upper <- min((ci/ui)[ui < 0])
        lower <- max((ci/ui)[ui > 0])
        out <- optimize(f = f, interval = c(lower, upper), 
                        tol = .Machine$double.eps)
        qout <- c(out$minimum, 1 - out$minimum)
      }
      else {
        out <- constrOptim(theta = theta, f = f, grad = if (useGradient) 
          gr
          else NULL, ui = ui, ci = ci)
        qout <- c(out$par, 1 - sum(out$par))
      }
    }
    names(qout) <- names(possAlleles)
    return(qout)
  }
  TotalStep <- function(q) {
    MStep(EStep(q), qinit = q)
  }
  if (is.null(alleleFreq)) {
    qOld <- init
    qNew <- TotalStep(qOld)
    allQ <- cbind(qOld, qNew)
    nIter <- 1
    while (sum(abs(qNew - qOld)) >= tol & nIter <= maxiter) {
      qOld <- qNew
      pMat <- EStep(qOld)
      qNew <- MStep(pMat, qinit = qOld)
      nIter <- nIter + 1
      allQ <- cbind(allQ, qNew)
    }
    pMat <- EStep(qNew)
  }
  else {
    qNew <- MStep(qinit = init)
    pMat <- NA
    nIter <- NA
    allQ <- NA
  }
  names(qNew) <- names(possAlleles)
  return(list(q = qNew, perLocationProb = pMat, alleleSet = possAlleles, 
              optimDetails = list(nIter = nIter, initial = init, pathOfQ = allQ)))
}

# eventTiming.WGD ---------------------------------------------------------

eventTiming.WGD <- function (x, m, history, totalCopy, method = c("fullMLE", "partialMLE", 
                                                                  "Bayes"), type = c("gain", "CNLOH"), seqError = 0, bootstrapCI = NULL, 
                             B = if (method == "Bayes") 10000 else 500, CILevel = 0.95, 
                             normCont = 0, verbose = TRUE, returnAssignments = FALSE, 
                             coverageCutoff = 1, minMutations = 10, init = NULL, maxiter = 100, 
                             tol = 1e-04, mutationId = 1:length(x), ...) 
{
  
  method <- match.arg(method)
  doCI <- !is.null(bootstrapCI)
  type <- match.arg(type)
  nMuts <- length(x)
  
  if (length(m) != nMuts) 
    stop("x and m must be of same length")
  if (length(mutationId) != length(unique(mutationId)) & verbose) 
    warning("mutationId not unique values; using index as id instead.")
  if (length(mutationId) != nMuts & verbose) 
    warning("mutationId an invalid length; using index as id instead.")
  if (returnAssignments) 
    returnData <- TRUE
  else returnData <- FALSE
  
  A <- history
  nCopies <- totalCopy
  if (is.null(dim(A))) 
    stop("'history' should be a matrix of size (nEvents +1) x (nEvents +1)")
  if (ncol(A) != nrow(A)) 
    stop("'history' should be a square matrix (hint: do not include the allele frequency 1 except for CNLOH events)")
  K <- ncol(A) - 1
  
  if (length(normCont) != 1 || normCont > 1 || normCont < 0) 
    stop("normCont must be given and should be single value between 0 and 1.")
  
  possAlleles <- allAF(nCopies, normCont = normCont, type = "mutation")[[1]]
  possAlleles = head(possAlleles, -2)
  
  if (nrow(A) != length(possAlleles)) 
    stop(length(possAlleles), "possible alleles for this scenario, but history only has", 
         nrow(A), "rows.")
  possAlleles <- errorAF(possAlleles, seqError = seqError)
  
  if (coverageCutoff < 1) 
    stop("coverageCutoff must be at least 1")
  nBelowCutoff <- length(which(m < coverageCutoff))
  nZero <- length(which(m == 0))
  if (any(m < coverageCutoff)) {
    if (verbose) 
      warning(nBelowCutoff, "mutations present below cutoff", 
              coverageCutoff, ", will be ignored")
    if (verbose & nZero > 0) 
      warning(nZero, "mutations have no reads at all")
    x <- x[which(m >= coverageCutoff)]
    mutationId <- mutationId[which(m >= coverageCutoff)]
    m <- m[which(m >= coverageCutoff)]
  }
  
  summaryTable <- cbind(nMuts, length(m), nZero, nBelowCutoff)
  colnames(summaryTable) <- c("Total Mutations", "Qualify for fitting", 
                              "Zero Read Coverage", "Coverage Below Cutoff")
  summaryTable <- t(summaryTable)
  colnames(summaryTable) <- "Number"
  call <- list(alleleSet = possAlleles, history = history, 
               totalCopy = totalCopy, type = type, method = method, 
               seqError = seqError, normCont = normCont, coverageCutoff = coverageCutoff, 
               minMutations = minMutations, B = B, init = init, maxiter = maxiter, 
               tol = tol)
  if (returnData) 
    inputData = data.frame(mutationId = mutationId, x = x, 
                           m = m)
  failReason <- NULL
  success <- TRUE
  if (length(m) < minMutations) {
    failReason <- paste(failReason, "not enough mutations above cutoff that satisfy coverage criteria", 
                        sep = ",")
    if (verbose) 
      warning("less than ", minMutations, " mutated locations meet coverage criteria; no estimate of pi will be calculated")
    success <- FALSE
  }
  dummyOutput <- list(pi = rep(NA, length = ncol(A)), alleleSet = possAlleles, 
                      optimDetails = list(nIter = NA, initial = NA, pathOfQ = NA), 
                      perLocationProb = NA, q = NA)
  names(dummyOutput$pi) <- paste("Stage", 0:K, sep = "")
  ci <- matrix(rep(NA, length = 2 * length(dummyOutput$pi)), 
               ncol = 2)
  row.names(ci) <- names(dummyOutput$pi)
  dummyOutput$piCI <- ci
  if (!is.null(failReason)) {
    output <- dummyOutput
  }
  else {
    rankA <- qr(A)$rank
    piId <- rankA == ncol(A)
    piEst <- rep(NA, ncol(A))
    if (piId) {
      if (method %in% c("fullMLE", "Bayes")) {
        output <- try(estimateQ.WGD(x, m, possAlleles, history = history, 
                                    init = init, maxiter = maxiter, tol = tol))
        if (output$optimDetails$nIter == maxiter) {
          failReason <- paste(failReason, "hit maximum number of iterations", 
                              sep = ",")
          success <- FALSE
        }
      }
      if (method == "partialMLE") {
        mleAlleles <- mleAF(x, m, totalCopy = nCopies, 
                            seqError = seqError, normCont = normCont, maxCopy = if (type == 
                                                                                    "gain") 
                              totalCopy - 1
                            else totalCopy)
        output <- try(.estimateQ(x, m, possAlleles, alleleFreq = mleAlleles$alleleSet$Freq, 
                                 history = history, init = init, maxiter = maxiter, 
                                 tol = tol))
      }
      if (!inherits(output, "try-error")) {
        piEst <- solve(A) %*% output$q
        piEst <- as.vector(piEst)
        piEst <- piEst/sum(piEst)
        if (method == "Bayes") {
          initBayes <- piEst
          bayesout <- .bayesEstimate(x, m, alleleSet = possAlleles, 
                                     history = A, init = initBayes, nsim = B, 
                                     CILevel = CILevel, ...)
          piEst <- bayesout$pi
          q <- A %*% piEst
          q <- q/sum(q)
          output <- list(piCI = bayesout$piCI, q = q, 
                         perLocationProb = NA, alleleSet = possAlleles, 
                         optimDetails = list(nIter = NA, initial = rbind(initBayes, 
                                                                         bayesout$mode), pathOfQ = NA))
        }
      }
      else {
        if (method == "Bayes") {
          initBayes <- seq(1, length(possAlleles))/length(possAlleles)
          bayesout <- .bayesEstimate(x, m, alleleSet = possAlleles, 
                                     history = A, init = initBayes, nsim = B, 
                                     ...)
          piEst <- bayesout$pi
          q <- A %*% piEst
          q <- q/sum(q)
          output <- list(piCI = bayesout$piCI, q = q, 
                         perLocationProb = NA, alleleSet = possAlleles, 
                         optimDetails = list(nIter = NA, initial = rbind(initBayes, 
                                                                         bayesout$mode), pathOfQ = NA))
        }
        else {
          output <- dummyOutput
          failReason <- paste(failReason, "error returned in estimating q", 
                              sep = ",")
          success <- FALSE
        }
      }
      if (method == "Bayes") 
        row.names(output$optimDetails$initial) <- c("mle", 
                                                    "postMode")
    }
    else {
      output <- dummyOutput
      failReason <- paste(failReason, paste("history matrix not invertible, rank=", 
                                            rankA, sep = ""), sep = ",")
      success <- FALSE
      if (verbose) 
        warning(paste("history matrix not invertible, rank=", 
                      rankA, sep = ""))
    }
    names(piEst) <- paste("Stage", 0:K, sep = "")
    output$pi <- piEst
    if (!returnAssignments) 
      output <- output[-grep("perLocationProb", names(output))]
    else {
      output[["perLocationProb"]] <- data.frame(inputData, 
                                                output[["perLocationProb"]], check.names = FALSE)
    }
    if (doCI & method != "Bayes") {
      if (!any(is.na(output$pi))) {
        bootsample <- bootstrapEventTiming(call = call, 
                                           B = B, x = x, m = m, type = bootstrapCI, pi = output$pi)
        ci <- t(apply(bootsample, 2, quantile, c((1 - 
                                                    CILevel)/2, 1 - (1 - CILevel)/2)))
      }
      else {
        ci <- matrix(rep(NA, length = 2 * length(output$pi)), 
                     ncol = 2)
        row.names(ci) <- names(output$pi)
      }
      output$piCI <- ci
      output <- output[c("pi", "piCI", "q", "perLocationProb", 
                         "optimDetails")]
    }
    else {
      if (method != "Bayes") 
        output <- output[c("pi", "q", "perLocationProb", 
                           "optimDetails")]
      else {
        if (!any(is.na(output$pi))) {
          row.names(output$piCI) <- names(output$pi)
        }
        output <- output[c("pi", "piCI", "q", "perLocationProb", 
                           "optimDetails")]
      }
    }
  }
  return(c(output, list(summaryTable = summaryTable, success = success, 
                        failReason = if (!is.null(failReason)) gsub("^,", "", 
                                                                    failReason) else failReason, call = call)))
}


# bootstrapEventTiming.WGD ------------------------------------------------


bootstrapEventTiming = function (eventOrdering, B, type = c("parametric", "nonparametric"), 
                                 pi, x, m, call) 
{
  type <- match.arg(type)
  if (missing(x) | missing(m)) {
    if (!"perLocationProb" %in% names(eventOrdering)) 
      stop("perLocationProb must be one of the named elements of 'eventOrdering' if 'x' and 'm' are not supplied.")
    x <- eventOrdering$perLocationProb$x
    m <- eventOrdering$perLocationProb$m
  }
  if (missing(pi)) {
    if (!"pi" %in% names(eventOrdering) | (any(is.na(eventOrdering$pi)) & 
                                           type == "parametric")) 
      stop("invalid values of 'pi' in 'eventOrdering' object")
    pi <- eventOrdering$pi
  }
  if (missing(call)) {
    if (!"call" %in% names(eventOrdering)) 
      stop("invalid values of 'call' in 'eventOrdering' object")
    call <- eventOrdering$call
  }
  N <- length(m)
  eventFunctionNames <- c("history", "totalCopy", "method", 
                          "seqError", "type", "normCont", "coverageCutoff", "minMutations", 
                          "init", "maxiter", "tol")
  requiredNames <- c("alleleSet", eventFunctionNames)
  if (!all(requiredNames %in% names(call))) 
    stop("Missing call information:", requiredNames[!requiredNames %in% 
                                                      names(call)])
  possAlleles <- call$alleleSet
  A <- call$history
  initial <- call$init
  if (type == "parametric") {
    if (ncol(A) != length(pi)) 
      stop("'history' from 'call' does not match length of 'pi'")
    q <- as.vector(A %*% pi)
    q <- q/sum(q)
    if (length(q) != length(possAlleles)) 
      stop("'alleleSet' and 'history' from 'call' does not match in size ")
    if (any(q < 0)) {
      if (any(q < -1e-10)) 
        stop("Programming error -- negative probabilities")
      else q[q < 0] <- 0
    }
    if (any(q > 1)) {
      if (any(q > 1 + 1e-10)) 
        stop("Programming error -- >1 probabilities")
      else q[q > 1] <- 1
    }
    bootData <- .parametricSimulation(B, m, q, alleleSet = possAlleles, 
                                      onlyReturnData = TRUE)
  }
  else {
    wh <- matrix(sample(1:N, N * B, replace = TRUE), nrow = N, 
                 ncol = B)
    bootData <- lapply(1:B, function(kk) {
      cbind(nMutAllele = x[wh[, kk]], nReads = m[wh[, kk]])
    })
  }
  piBoot <- lapply(bootData, function(z) {
    do.call(eventTiming.WGD, c(list(x = z[, "nMutAllele"], m = z[, 
                                                                 "nReads"]), call[eventFunctionNames]))$pi
  })
  piBoot <- do.call("rbind", piBoot)
  return(piBoot)
}


# eventTimingOverList.WGD -------------------------------------------------

# dfList = list of dataframes per sample, each with the following columns: segId, nMutAllele, nReads, type, mutationId
# "type" includes: CNLOH, SingleGain, DoubleGain, WGDs
# normCont = normal contamination
# eventArgs = arguments passed to eventTiming.WGD, e.g. minMutations

eventTimingOverList.WGD <- function (dfList, normCont, eventArgs = NULL) 
{
  # check all required column names
  if (any(sapply(dfList, function(x) {
    any(!c("segId", "type", "nMutAllele", "nReads", "mutationId") %in% 
        names(x))
  }))) 
    stop("dfList elements are missing required column names")
  
  # check that segment ID's do not contain full stops
  if (any(sapply(dfList, function(x) {
    length(grep("[.]", x$segId)) > 0
  }))) 
    stop("segId cannot contain a period")
  
  # If list names missing, add them as Sample 1, Sample 2, etc. 
  if (is.null(names(dfList))) 
    names(dfList) <- paste("Sample", 1:length(dfList), sep = "")
  
  # Need normal contamination per sample
  if (length(normCont) != length(dfList)) 
    stop("dfList and normCont must be of the same length (1 per sample)")
  
  # for each type of event, get segments
  whCall <- lapply(dfList, function(x) {
    tapply(x$segId, factor(x$type, levels = c("Other", "CNLOH", 
                                              "SingleGain", "Diploid", "DoubleGain", "WGD")), unique, 
           simplify = FALSE)
  })
  
  singleSampleFunction <- function(x, nc, dat) {
    
    # Time CNLOH
    if (length(x[["CNLOH"]]) > 0) {
      ACNLOH <- makeEventHistory(totalCopy = 2, type = "LOH")[[1]]
      eventCNLOH <- lapply(x[["CNLOH"]], function(segId) {
        subdat <- dat[which(dat$segId == segId), ]
        out <- do.call(eventTiming, c(list(x = subdat$nMutAllele, 
                                           m = subdat$nReads, history = ACNLOH, totalCopy = 2, 
                                           type = "CNLOH", mutationId = subdat$mutationId, 
                                           normCont = nc), eventArgs))
      })
      names(eventCNLOH) <- x[["CNLOH"]]
    }
    else eventCNLOH <- NULL
    
    # Time SG
    if (length(x[["SingleGain"]]) > 0) {
      AGain <- makeEventHistory(totalCopy = 3, type = "gain")[[1]]
      eventGain <- lapply(x[["SingleGain"]], function(segId) {
        subdat <- dat[which(dat$segId == segId), ]
        do.call(eventTiming, c(list(x = subdat$nMutAllele, 
                                    m = subdat$nReads, history = AGain, totalCopy = 3, 
                                    type = "gain", mutationId = subdat$mutationId, 
                                    normCont = nc), eventArgs))
      })
      names(eventGain) <- x[["SingleGain"]]
    }
    else eventGain <- NULL
    
    # Time Double Gains
    if (length(x[["DoubleGain"]]) > 0) {
      ADGain <- makeEventHistory(totalCopy = 4, type = "gain")[[1]]
      eventDGain <- lapply(x[["DoubleGain"]], function(segId) {
        subdat <- dat[which(dat$segId == segId), ]
        do.call(eventTiming, c(list(x = subdat$nMutAllele, 
                                    m = subdat$nReads, history = ADGain, totalCopy = 4, 
                                    type = "gain", mutationId = subdat$mutationId, 
                                    normCont = nc), eventArgs))
      })
      names(eventDGain) <- x[["DoubleGain"]]
    }
    else eventDGain <- NULL
    
    # Time WGD segments
    if (length(x[["WGD"]]) > 0) {
      WGDGain <- matrix(c(0,4,2,0), byrow=TRUE, ncol=2)
      eventWGDGain <- lapply(x[["WGD"]], function(segId) {
        subdat <- dat[which(dat$segId == segId), ]
        do.call(eventTiming.WGD, c(list(x = subdat$nMutAllele, 
                                        m = subdat$nReads, history = WGDGain, totalCopy = 4, 
                                        type = "gain", mutationId = subdat$mutationId, 
                                        normCont = nc), eventArgs))
      })
      names(eventWGDGain) <- x[["WGD"]]
    }
    else eventWGDGain <- NULL
    
    
    return(list(SingleGain = eventGain, CNLOH = eventCNLOH, 
                DoubleGain = eventDGain, WGD = eventWGDGain))
  }
  
  mapply(whCall, normCont, dfList, FUN = singleSampleFunction, 
         SIMPLIFY = FALSE)
}



.parametricSimulation = function (B, m, q, alleleSet, onlyReturnData = FALSE) 
{
  if (length(q) != length(alleleSet)) 
    stop("'alleleSet' and 'q' must be of the same length")
  if (is.null(dim(m))) {
    nMut <- length(m)
    m <- matrix(m, ncol = B, nrow = nMut, byrow = FALSE)
  }
  else {
    nMut <- nrow(m)
    if (B != ncol(m)) 
      stop("'B' must be same as number of columns of m")
  }
  if (nMut == 0) {
    data.frame(AF = rep(NA, nMut), nReads = rep(NA, nMut), 
               nMutAllele = rep(NA, nMut), obsAF = rep(NA, nMut))
  }
  nAlleles <- length(q)
  nPerCategory <- rmultinom(n = B, size = nMut, prob = q)
  AF <- apply(nPerCategory, 2, function(x) {
    sample(rep(alleleSet, times = x))
  })
  if (is.null(dim(AF))) {
    AF <- matrix(AF, ncol = B, nrow = nMut, )
  }
  return(lapply(1:B, function(kk) {
    x <- rbinom(nMut, size = m[, kk], prob = AF[, kk])
    if (onlyReturnData) cbind(nMutAllele = x, nReads = m[, 
                                                         kk]) else data.frame(AF = AF[, kk], nMutAllele = x, 
                                                                              nReads = m[, kk], obsAF = x/m[, kk])
  }))
}
