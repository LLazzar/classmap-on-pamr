daisy_vcr <- function(x, daisy.out = NULL, type = NULL,
                     weights = NULL, mybounds = NULL,
                     warnBin = FALSE, stand = TRUE) {
  #
  # To do:
  # - let daisy_vcr take X, Xnew and then compute the rectangular
  #   matrix of dissimilarities that we need, saving time and
  #   storage space. Jakob suggests that the data.matrix()
  #   mechanism needs to be replaced for this.
  #
  # Function that computes dissimilarities, based on cluster::daisy.
  # In case daisy.out is not NULL, x is considered to be new data,
  # and the information on the variables is extracted from daisy.out.
  # First, compatibility of x with daisy.out$x is checked.
  #
  # Type should be a list with elements:
  # $asymm, $symm, $ordratio, $logratio, $nominal.
  # Variable types not declared explicitly are assumed to be numeric.
  # [Peter] In fact, "$ordratio" is a strange name, it should just
  #   be "$ordinal". I modified the code to accept "$ordinal" too.
  #   Also, $logratio really just means take the logarithm, and this
  #   is superfluous in R where such a pretransformation is trivial.
  #   In general I think that type2 is much easier to work with than
  #   type, but we have to be backward compatible here. Maybe not in
  #   the side project later, that function will get a different
  #   name anyway.
  #
  # Options that were removed: metric (you can use dist() for
  # "euclidean" and "manhattan"), stand (you can standardize before
  # calling daisy_vcr), and most warn* (Jakob wants the warnings).
  # [Peter] I've put back one warning option (warnBin=F) as this
  #   warning would confuse some users, and it is a harmless
  #   occurrence since the result is the same anyway.
  #
  # Arguments:
  # x         : data matrix, in which the variables may be of
  #             mixed types.
  # daisy.out : if not NULL, the output of daisy_vcr on training
  #             data. If NULL, x is treated as training data.
  # type      : for training data, a list with variable types
  #             containing variable names (or their column numbers).
  #             If NULL, all variables are treated as numeric.
  #             For new data, the types of the training variables
  #             are used.
  # weights   : For training data, these are the variable weights,
  #             and by default they are all 1.
  #             For new data, the weights of the training variables
  #             are used.
  # mybounds  : a numeric vector c(a, b) of length two, with a < b.
  #             Only used for new data. After scaling a numeric
  #             variable (by the parameters it has in the training
  #             data), there may be values below 0 or above 1.
  #             These are then truncated to lie in the interval
  #             [a, b]. By default no truncation takes place.
  # warnBin   : if TRUE, a warning is given whenever a numeric
  #             variable takes on only 2 values.
  #
  if (length(dim(x)) != 2 && !(is.data.frame(x))) {
    stop("x is not a dataframe or a matrix.")
    }
  xorig <- x
  n <- nrow(x)
  p <- ncol(x)
  #
  if (is.null(daisy.out)) {
    # Checking the weights:
    if (is.null(weights)) weights <- 1
    if (length(weights) == 1) {
      weights <- rep.int(weights, p)
    } else if (length(weights) != p) {
      stop(paste0("'weights' must be of length ", p," (or 1)."))
    }
    if (sum(weights < 0) > 0) stop("There are negative weights.")
    if (is.null(type)) type <- list()
    #
  } else {
    # calculate dissimilarities on new data using
    # the result daisy.out from the training data:
    if (!is.null(weights)) warning(paste0(
      "\nFor new data, the weights of the training data are used.",
      "\nThe argument 'weights' in this call will not be used."))
    weights <- daisy.out$weights
    if (!is.null(type)) warning(paste0(
      "\nFor new data, the types of the training data are used.",
      "\nThe argument 'type' in this call will not be used."))
    type  <- daisy.out$type
    x_old <- daisy.out$x
    #
    # Check: is new data format compatible with training format?
    #
    allowNA <- TRUE # is NA allowed in test data?
    varsInModel <- which(daisy.out$weights > 0)
    checkXnew(x_old, x, allowNA = allowNA,
              varsInModel = varsInModel)
    # Match colnames to make sure they are in the same order, so
    # we can use the old weights, types, ranges, etc.:
    namesmatch <- match(colnames(x_old), colnames(x))
    x <- x[, namesmatch]
  }
  #
  # Now perform checks on the variable types:
  #
  varnms <- dimnames(x)[[2]]
  pColl <- function(n) paste(n, collapse = ", ")
  if (length(type)) {
    if (!is.list(type) || is.null(ntyp <- names(type)) ||
        any(ntyp == ""))
      stop(gettextf("invalid %s; must be named list",
                    sQuote("type")))
    typelist <- c("asymm", "symm", "nominal",
                 "logratio", "ordratio", "ordinal")
    for (nt in ntyp) { # e.g. nt = "nominal"
      if (!(nt %in% typelist)) stop(paste0(
        "type ",nt," is not allowed."))
      cvec <- type[[nt]]
      ct <- paste0("type$", nt)
      if (is.character(cvec)) {
        if (!is.null(varnms) && !all(cvec %in% varnms))
          stop(gettextf("%s has invalid column names",
                        ct))
      }
      else if (is.numeric(cvec)) {
        if (!all(1 <= cvec & cvec <= p))
          stop(gettextf("%s must be between 1 and ncol(x)",
                        ct))
      }
      else stop(gettextf("%s must contain column names or numbers",
                         ct))
    }
    tA <- type$asymm
    tS <- type$symm
    if (!is.null(tA) || !is.null(tS)) {
      d.bin <- cbind(as.data.frame(x[, tA, drop = FALSE]),
                     x[, tS, drop = FALSE])
      lenB <- sapply(lapply(d.bin, function(y) levels(as.factor(y))),
                     length)
      if (any(lenB > 2))
        stop("at least one binary variable has more than 2 levels.")
      if (!is.null(daisy.out) && any(lenB < 2))
        warning("at least one binary variable has < 2 levels.")
      if (any(is.f <- sapply(d.bin, is.factor)))
        d.bin[is.f] <- lapply(d.bin[is.f], function(f)
          as.integer(as.character(f)))
      if (!all(sapply(d.bin, function(y) is.logical(y) ||
                      all(sort(unique(as.numeric(y[!is.na(y)]))) %in%
                          0:1))))
        stop("at least one binary variable has values not in {0,1,NA}.")
    }
  }
  #
  # Convert the data to a numeric matrix, and create type2:
  #
  if (is.data.frame(x)) {
    type2 <- sapply(x, data.class)
    x <- data.matrix(x) # turns characters into numeric values.
  } else {
    type2 <- rep("numeric", p)
    names(type2) <- colnames(x)
  }
  if (length(type)) {
    tT <- type$ordratio
    if (is.null(tT)) tT <- type$ordinal # added
    x[, names(type2[tT])] <- unclass(as.ordered(x[, names(type2[tT])]))
    type2[tT] <- "T"
    tL <- type$logratio
    if (sum(as.vector(x[, names(type2[tL])]) <= 0) > 0) stop(
      "\nVariable(s) to be log-transformed are not always > 0.")
    x[, names(type2[tL])] <- log(x[, names(type2[tL])])
    type2[tL] <- "L"
    type2[tA] <- "A"
    type2[tS] <- "S"
    tN <- type$nominal
    type2[tN] <- "N"
  }
  type2[tI <- type2 %in% c("numeric", "integer")] <- "I"
  # here "I" stands for Interval-scaled.
  if (is.null(daisy.out) && warnBin && n > 9 && any(tI) &&
      any(iBin <- apply(x[, tI, drop = FALSE],
                        2, function(v) length(table(v)) == 2))) {
    warning(gettextf("\nNumeric variable(s) %s take only two values.",
                     pColl(which(tI)[iBin]))) }
  type2[type2 == "ordered"]   <- "O"
  type2[type2 == "factor"]    <- "N"
  type2[type2 == "character"] <- "N"
  if (any(ilog <- type2 == "logical")) {
    warning(sprintf(ngettext(
      sum(ilog),
      "setting 'logical' variable %s to type 'asymm'",
      "setting 'logical' variables %s to type 'asymm'"),
      pColl(which(ilog))), domain = NA)
    type2[ilog] <- "A"
  }
  #
  # Compute location and scale of the variables:
  #
  if (is.null(daisy.out)) {
    if (stand) {# classical Daisy standardization using the range
      colR   <- apply(x, 2, range, na.rm = TRUE)
      colmin <- colR[1, ]
    } else {# no standardization
      colR   <- rep(1, ncol(x))
      colmin <- rep(0, ncol(x))
    }

  } else {
    colR   <- daisy.out$colR
    colmin <- daisy.out$colmin
  }
  sx <- colR[2, ] - colmin
  if (any(sx == 0)) {sx[sx == 0] <- 1}
  #
  #
  for (j in seq_len(p)) { # iterate over the variables
    if (weights[j] > 0) { # skip variables with zero weight
      tempType <- type2[j]
      if (tempType %in% c("I", "T", "L")) {
        # Only Interval/Ordinal/Log variables are scaled:
        x[, j] <- scale(x[, j], center = colmin[j], scale = sx[j])

      } # ends if(type)
    } # ends if(weight > 0)
  } # ends loop over variables
  #
  # Actual distance computations:
  #
  nums   <- matrix(0, n, n) # numerator of the final dissim. matrix
  denoms <- matrix(0, n, n) # denominator of the final dissim. matrix
  for (j in seq_len(p)) { # iterate over the variables
    if (weights[j] > 0) { # skip variables with zero weight
      tempVar    <- x[, j]
      tempType   <- type2[j]
      denomstemp <- matrix(1, n, n) # by default, the denominator
      # is increased by 1 (we will multiply by weights[j] later).
      #
      # compute contribution to the numerator:
      if (tempType %in% c("A", "S", "N")) { # binary and nominal
        numtemp <- outer(tempVar, tempVar, "!=") + 0
      } else {# Interval/Ordinal/Log: Manhattan distance
        numtemp <- abs(outer(tempVar, tempVar, "-"))
      }
      # adjust denominator for asymmetric binary
      if (tempType == "A") {
        denomstemp <- (outer(tempVar, tempVar, "+") > 0) + 0
      }
      # handle NA's: the distance is 0, but the denominator does
      # not get incremented:
      numtemp[which(is.na(numtemp))] <- 0
      denomstemp[which(is.na(denomstemp))] <- 0
      #
      # use weights:
      numtemp    <- weights[j] * numtemp
      denomstemp <- weights[j] * denomstemp
      #
      # add to the overall distance matrix:
      nums   <- nums + numtemp
      denoms <- denoms + denomstemp
    }
  }
  disv <- as.dist(nums / denoms) # makes dissimilarity object
  #
  if (anyNA(disv)) attr(disv, "NA.message") <-
    "NA-values in the dissimilarity matrix !"
  class(disv) <-  c("dissimilarity","dist")
  attr(disv, "Labels") <- dimnames(x)[[1]]
  attr(disv, "Size")   <- n
  attr(disv, "Metric") <- "mixed"
  attr(disv, "Types")  <- type2
  #
  result <- list(disv = disv, # the dissimilarity object
                type = type, # types of the variables
                type2 = type2, # simpler version, as a check
                colR = colR, # column ranges (scales)
                colmin = colmin, # minimum (location) of each column
                weights = weights, # variable weights
                x = xorig)
  return(result)
}
