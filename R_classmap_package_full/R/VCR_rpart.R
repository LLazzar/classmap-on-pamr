vcr.rpart.train <- function(X, y, trainfit, type = list(),
                            k = 5, stand = TRUE) {
  #
  # Constructs vcr output from an rpart object.
  #
  # Arguments:
  # X        : A rectangular matrix or data frame, where the
  #            columns (bvariables) may be of mixed type and
  #            may contain NA's.
  # y        : factor with the given class labels.
  #            It is crucial that X and y are EXACTLY the same
  #            as in the call to rpart().
  #            Since rpart() allows y to contain NA's, we do too.
  #            If _all_ entries of y are NA we stop with a message,
  #            but in that situation rpart() crashes anyway.
  # k        : the number of nearest neighbors used in the
  #            farness computation.
  # trainfit : the output of an rpart() training cycle
  # other arguments: development
  #
  # Values include:
  #   yint        : given labels as integer 1, 2, 3, ...
  #   y           : given labels
  #   levels      : levels of y
  #   predint     : predicted class number. For each case this
  #                 is the class with the highest mixture density.
  #   pred        : predicted label of each object.
  #   altint      : competing class as integer, i.e. the non-given
  #                 class number with the highest mixture density.
  #   altlab      : competing class of each object.
  #   PAC         : PAC of each object
  #   farness     : farness of each object to its given class
  #   ofarness    : For each object i, its lowest farness to any
  #                 class (including its own).

  # Auxiliary function:
  #
  importanceWeights <- function(varimp, cnames) {
    # Takes variable importance from the package rpart or the
    # package randomForest.
    # In varimp the variables are listed from highest to lowest
    # importance. The function rearranges them in the order of
    # the original variables, where cnames contains the column
    # names of the training data.
    # Variables not occurring in varimp are assigned weight zero.
    #
    varimp <- pmax(varimp, 0)
    varimp <- drop(varimp / sum(varimp))
    namesimp <- names(varimp)
    nonzeroVars <- which(cnames %in% namesimp)
    varweights <- rep(0, length(cnames))
    for (j in nonzeroVars) { # j=2
      varnamej <- cnames[j]
      wgt <- varimp[which(namesimp == varnamej)]
      varweights[j] <- wgt
    }
    return(varweights)
  }

  # Here the main code starts
  if (nrow(X) == 1) X <- t(X)
  if (length(dim(X)) != 2 ||
      !(is.data.frame(X) || is.numeric(X)))
    stop("x is not a dataframe or a numeric matrix.")
  n <- nrow(X)
  if (n < 2) stop("The training data should have more than one case.")
  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = TRUE)
  #                      stop.training.if.a.response.is.NA = F)
  # If it did not stop: yint has length n, and its values are in
  # 1, ..., nlab without gaps. It is NA where y is: outside indsv.
  #
  lab2int <- checked$lab2int # is a function
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y)
  indsv   <- checked$indsv # cases in indsv can be visualized
  nvis    <- length(indsv)
  yintv   <- yint[indsv]
  #
  # As in vcr.svm.train, treating X as "new data" allows us to
  # obtain predictions also for cases whose y was NA:
  preds <- predict(object = trainfit, type = "class",
                   newdata = as.data.frame(X))
  predint <- lab2int(preds)
  #
  # From here on we only use the cases with available y
  #
  probs <- predict(object = trainfit, type = "prob",
                   newdata = as.data.frame(X))[indsv, ]
  probs <- probs[, order(lab2int(colnames(probs)))]
  # rowSums(probs) # all 1
  PAC <- altint <- rep(NA, n)
  # We can only compute PAC and altint where yint is not NA:
  altintv <- PACv <- rep(NA, nvis)
  for (i in seq_len(nvis)) {
    ptrue      <- probs[i, yintv[i]]
    others     <- (seq_len(nlab))[-yintv[i]]
    palt       <- max(probs[i, others])
    altintv[i] <- others[which.max(probs[i, others])]
    PACv[i]    <- palt / (ptrue + palt)
  }
  #
  PAC[indsv]    <- PACv    # the other entries of PAC stay NA
  altint[indsv] <- altintv # the other entries of altint stay NA
  if (nlab == 2) altint <- 3 - yint # for consistency with svm etc.
  if (is.null(k)) {k <- 5}
  #
  # Compute farness based on knn
  #
  ### Construct initial daisy matrix
  #
  varweights <- importanceWeights(trainfit$variable.importance,
                                  colnames(X))
  X <- X[indsv, ]
  daisy.out <- daisy_vcr(X, type = type, weights = varweights,
                         stand = stand)
  # daisy_vcr skips the variables with zero weight.
  dismat  <- as.matrix(daisy.out$disv)
  meandis <- mean(as.vector(dismat), na.rm = TRUE)
  dismat[which(is.na(dismat))] <- meandis
  #
  # Compute neighbors by sorting dissimilarities:
  sortNgb <- t(apply(dismat, 1, order))[, -1]
  sortDis <- t(apply(dismat, 1, sort))[, -1]
  #
  # Compute initial fig[i, g] from cases to classes:
  #
  fig <- matrix(rep(NA, nvis * nlab), ncol = nlab)
  for (i in seq_len(nvis)) { # loop over all cases in indsv
    for (g in seq_len(nlab)) { # loop over classes
      ngbg <- which(yintv[sortNgb[i, ]] == g)
      if (length(ngbg) > k) {ngbg <- ngbg[seq_len(k)]}
      fig[i, g] <- median(sortDis[i, ngbg])
    }
  }
  #
  # Compute farness:
  #
  farout <- compFarness(type = "knn", testdata = FALSE, yint = yintv,
                        nlab = nlab, X = NULL, fig = fig,
                        d = NULL, figparams = NULL)
  fig <- matrix(rep(0, n * nlab), ncol = nlab)
  fig[indsv, ] <- farout$fig
  farness <- ofarness <- rep(NA, n)
  farness[indsv] <- farout$farness
  ofarness[indsv] <- farout$ofarness
  figparams <- farout$figparams
  figparams$daisy.out <- daisy.out
  figparams$meandis <- meandis
  figparams$k       <- k
  #
  return(list(X         = X,
              yint      = yint,
              y         = levels[yint],
              levels    = levels,
              predint   = predint,
              pred      = levels[predint],
              altint    = altint,
              altlab    = levels[altint],
              PAC       = PAC,
              figparams = figparams,
              fig       = fig,
              farness   = farness,
              ofarness  = ofarness,
              trainfit  = trainfit))
}


vcr.rpart.newdata <- function(Xnew, ynew = NULL,
                              vcr.rpart.train.out,
                              LOO = FALSE) {
  #
  # Predicts class labels for new data by rpart, using the output
  # of vcr.svm.train() on the training data.
  # For new data cases whose label in ynew is non-missing, additional
  # output is produced for constructing graphical displays.
  #
  # Arguments:
  #   Xnew               : data matrix of the new data, with the same
  #                        number of columns d as in the training data.
  #                        Missing values are not allowed.
  #   ynew               : factor with class membership of each new
  #                        case. Can be NA for some or all cases.
  #                        If NULL, is assumed to be NA everywhere.
  #   vcr.rpart.train.out : output of vcr.rpart.train() on the training
  #                        data. This also contains nlab, and levels.
  #   LOO                : leave one out. Only used when testing this
  #                        function on a subset of the training data.
  #                        Default is LOO=F.
  #
  # Values include:
  #   yintnew   : given labels as integer 1, 2, 3, ..., nlab.
  #               Can be NA.
  #   ynew      : labels if yintnew is available, else NA.
  #   levels    : levels of the response, from vcr.rpart.train.out
  #   predint   : predicted class number, always exists.
  #   pred      : predicted label of each object
  #   altint    : alternative label as integer if yintnew was given,
  #               else NA
  #   altlab    : alternative label if yintnew was given, else NA
  #   PAC       : probability of alternative class for each object
  #               with non-NA yintnew
  #   farscales : (from training data) scales to use in fig
  #   fig       : distance of each object i from each class g.
  #               Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any
  #               class (including its own). Always exists.
  #
  X <- vcr.rpart.train.out$X # training data with all variables
  # d <- ncol(X)
  trainfit <- vcr.rpart.train.out$trainfit
  varsInModel <- attr(trainfit$terms, "term.labels")
  checkXnew(Xtrain = X, Xnew = Xnew, allowNA = TRUE,
            varsInModel = varsInModel)
  # From here on we know that the variable names match, but
  # they might be in a different order. To correct for this:
  namesmatch <- match(colnames(X), colnames(Xnew))
  Xnew <- Xnew[, namesmatch]
  #
  yint   <- vcr.rpart.train.out$yint # class membership in training data
  indold <- which(!is.na(yint))
  yint   <- yint[indold]
  ntrain <- length(indold)
  n      <- nrow(Xnew) # number of new cases
  levels <- vcr.rpart.train.out$levels # same as in training data
  nlab   <- length(levels) # number of classes
  #
  if (is.null(ynew)) ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = FALSE, levels)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv   # INDiceS of cases we can Visualize,
  #                         # can even be empty, is OK.
  nvis    <- length(indsv)   # can be zero.
  yintnew <- rep(NA, n)       # needed to overwrite "new" levels:
  yintnew[indsv] <- lab2int(ynew[indsv])
  # Any "new" levels are now set to NA.
  yintv   <- yintnew[indsv]  # entries of yintnew that can be visualized
  # From here on yintnew has length n, with values in seq_len(nlab).
  #
  newdata <- as.data.frame(Xnew)
  preds <- predict(trainfit, newdata, type = "class")
  predint <- lab2int(preds)
  #
  # From here on we only use the cases with available y
  #
  PAC <- altint <- rep(NA, n)
  probs <- predict(object = trainfit, newdata,
                  type = "prob")[indsv, ]
  probs <- probs[, order(lab2int(colnames(probs)))]
  # rowSums(probs) # all 1
  # We can only compute PAC and altint where yint is not NA:
  if (nvis > 0) { # if there are cases to visualize:
    altintv <- PACv <- rep(NA, nvis)
    for (i in seq_len(nvis)) {
      ptrue      <- probs[i, yintv[i]]
      others     <- (seq_len(nlab))[-yintv[i]]
      palt       <- max(probs[i, others])
      altintv[i] <- others[which.max(probs[i, others])]
      PACv[i]    <- palt / (ptrue + palt)
    }
    PAC[indsv]    <- PACv    # the other entries of PAC stay NA
    altint[indsv] <- altintv # the other entries of altint stay NA
  }
  if (nlab == 2) altint <- 3 - yintnew # for consistency with svm etc.
  #
  # Compute farness
  figparams <- vcr.rpart.train.out$figparams
  k         <- vcr.rpart.train.out$figparams$k
  #
  # Construct dissimilarity matrix
  #
  Xall <- rbind(X, Xnew) # has ntrain + n rows
  # Here we have used that the variables of Xnew are in the
  # same order as in X.
  indn <- (ntrain + 1):(ntrain + n) # range of new rows
  daisynew <- daisy_vcr(Xall, daisy.out = figparams$daisy.out)
  dismatnew <- as.matrix(daisynew$disv)[indn, -indn, drop = FALSE]
  dismatnew[is.na(dismatnew)] <- figparams$meandis
  sortNgb <- t(apply(dismatnew, 1, order))
  sortDis <- t(apply(dismatnew, 1, sort))
  if (LOO) sortNgb <- sortNgb[, -1] # leave one out, only for testing
  if (LOO) sortDis <- sortDis[, -1] # leave one out, only for testing
  #
  # Compute initial fig[i,g] from all new cases to classes:
  #
  fignew <- matrix(rep(NA, n * nlab), ncol = nlab)
  for (i in seq_len(n)) {
    for (g in seq_len(nlab)) {
      ngbg <- which(yint[sortNgb[i, ]] == g)
      if (length(ngbg) > k) {ngbg <- ngbg[seq_len(k)]}
      fignew[i, g] <- median(sortDis[i, ngbg])
    }
  }
  farout <- compFarness(type = "knn", testdata = TRUE, yint = yintnew,
                       nlab = nlab, X = NULL, fig = fignew,
                       figparams = vcr.rpart.train.out$figparams)

  return(list(yintnew = yintnew,
              ynew = levels[yintnew],
              levels = levels,
              predint = predint,
              pred = levels[predint],
              altint = altint,
              altlab = levels[altint],
              PAC = PAC,
              fig = farout$fig,
              farness = farout$farness,
              ofarness = farout$ofarness))
}
