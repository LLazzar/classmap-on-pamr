vcr.da.train <- function(X, y, rule = "QDA",
                         estmethod = "meancov") {
  #
  # Custom DA function which prepares for graphical displays.
  # The disciminant analysis itself is carried out by the maximum
  # a posteriori rule which maximizes the density of the mixture.
  #
  # Arguments:
  #   X         : data matrix. Missing values are not allowed.
  #   y         : vector with the given class labels.
  #               Missing y are not allowed, such cases should be
  #               removed first.
  #   rule      : either "QDA" for quadratic discriminant analysis
  #               or "LDA" for linear discriminant analysis.
  #   estmethod : function for location and covariance estimation.
  #               Should return a list with $m and $S. Can be "meancov"
  #               (classical mean and covariance matrix) or "DetMCD".
  #
  # Returns:
  #   yint      : given labels as integer 1, 2, 3, ...
  #   y         : given labels
  #   levels    : levels of y
  #   predint   : predicted class number. For each case this
  #               is the class with the highest mixture density.
  #   pred      : predicted label of each object.
  #   altint    : alternative class as integer, i.e. the non-given
  #               class number with the highest mixture density.
  #   altlab    : alternative class of each object.
  #   fig       : Mahalanobis distance of each object i to each class g
  #   classMS   : list with center and covariance matrix of each class
  #   lCurrent  : log(mixture density) in object's given class.
  #   lPred     : log(mixture density) in object's predicted class.
  #   lAlt      : log(mixture density) in object's alternative class.
  #   PAC       : probability of alternative class for each object
  #   figparams : parameters for computing fig, can be used for
  #               new data.
  #   fig       : distance of each object i from each class g.
  #   farness   : farness of each object to its given class.
  #   ofarness  : For each object i, its lowest farness to any class.
  #               If much above 1, we suspect i may be an outlier.
  #
  X <- as.matrix(X) # in case it is a data frame
  if (nrow(X) == 1) X <- t(X)
  if (sum(is.na(as.vector(X))) > 0) {
    stop("The data matrix X has NA's.")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2) stop("The training data should have more than one case.")
  if (!(rule %in% c("LDA", "QDA"))) {
    stop(" rule should be either LDA or QDA")
  }
  if (!is.function(estmethod)) {
    if (estmethod == "meancov") {
      estmethod <- meancov
    } else if (estmethod == "DetMCD") {
      estmethod <- DetMCD
    }
  }

  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = T)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y)
  yintv   <- yint[indsv]
  classSizes <- rep(0, nlab)
  for (g in seq_len(nlab)) {classSizes[g] <- sum(yintv == g)}
  # classSizes
  MD2s    <- matrix(0, n, nlab) # squared mahalanobis distances
  mixf    <- matrix(0, n, nlab) # component of mixture density
  lmixf   <- matrix(0, n, nlab) # the logarithm of that density
  classMS <- list()
  #
  if (rule == "QDA") {
    tinyclasses <- which(classSizes < (d + 1)) # was 2*(d+1)
    if (length(tinyclasses) > 0) {
      stop(paste0(" The following classes have too few members",
                  " to compute a covariance matrix:",
                  levels[tinyclasses]))
    }
    # Compute the center and covariance matrix of each class:
    iCN <- rep(NA, nlab) # for inverse condition numbers
    for (g in seq_len(nlab)) {
      clinds       <- indsv[which(yintv == g)]
      # class indices in 1,...,n
      class.out    <- estmethod(X[clinds, ]) # contains mu and Sigma
      iCN[g]       <- numericalCond(class.out$S)
      classMS[[g]] <- class.out # now has $m, $S
    }
    tinyCN <- which(iCN < 1e-06)
    if (length(tinyCN) > 0) {
      stop(paste0("\n The following classes have an ill-conditioned",
                  " covariance matrix: \n", levels[tinyCN],
                  "\n You may want to try LDA or another classifier"))
    } # From here onward we can use the inverse of each S.
  } # end of QDA
  #
  if (rule == "LDA") { # Then we compute a single S instead:
    Xv <- X[indsv, ] # only the rows with available y
    for (g in seq_len(nlab)) {
      clinds  <- which(yintv == g) # indices in 1, ..., indsv
      classMS[[g]] <- meancov(Xv[clinds, ]) # colMeans
      Xv[clinds, ] <- sweep(Xv[clinds, ], 2, classMS[[g]]$m) # centered
    }
    pool.out <- estmethod(Xv) # estimates m and S on pooled data.
    if (numericalCond(pool.out$S) < 1e-06) {
      stop(paste0("\n The pooled covariance matrix is ill-conditioned.",
                  "\n You may want to try a different classifier")) }
    for (g in seq_len(nlab)) {
      classMS[[g]]$m  <- classMS[[g]]$m + pool.out$m
      classMS[[g]]$S  <- pool.out$S
    }
  }
  # From here on everything is the same for LDA and QDA.
  #
  # compute MD2 of each object to all classes:
  for (g in seq_len(nlab)) { # g=1
    class.out <- classMS[[g]] # now nonsingular S
    MD2s[, g] <- mahalanobis(X, class.out$m, class.out$S)
  }
  # Compute prior probabilities:
  for (g in seq_len(nlab)) { # g=1
    clinds <- which(yintv == g) # indices in 1, ..., indsv
    classMS[[g]]$prior <- length(clinds) / length(indsv)
  }
  #
  for (g in seq_len(nlab)) { # g=1
    S          <- classMS[[g]]$S
    prior      <- classMS[[g]]$prior
    lmixf[, g] <- -0.5*(log(det(S)) + MD2s[, g]) + log(prior)
    mixf[, g]  <- exp(lmixf[, g])
  }
  # Compute predictions far all objects in 1, ..., n
  predint <- lPred <- rep(NA, n)
  lPred   <- apply(lmixf[, , drop = FALSE], 1, max)
  predint <- apply(lmixf[, , drop = FALSE], 1, which.max)
  #
  # Compute fCurrent and fAlt for all objects with available y:
  lCurrent <- fCurr <- fAlt <- lAlt <- altint <- rep(NA, n)
  for (g in seq_len(nlab)) { # g=1
    clinds <- indsv[which(yintv == g)] # indices in 1, ..., n
    others <- (seq_len(nlab))[-g] # alternative classes
    fCurr[clinds]    <- mixf[clinds, g]
    lCurrent[clinds] <- lmixf[clinds, g]
    fAlt[clinds]     <- apply(mixf[clinds, others, drop = FALSE],
                              1, max)
    lAlt[clinds]     <- apply(lmixf[clinds, others, drop = FALSE],
                              1, max)
    altint[clinds]   <- others[apply(lmixf[clinds, others, drop = FALSE],
                                     1, which.max)]
  }
  #
  # Compute PAC:
  #
  PAC <- rep(NA, n)
  PAC[indsv] <- fAlt[indsv] / (fCurr[indsv] + fAlt[indsv])
  #
  # Compute farness:
  #
  farout <- compFarness(type = "affine", testdata = F, yint = yint,
                        nlab = nlab, X = NULL, fig = sqrt(MD2s),
                        d = d, figparams = NULL)
  return(list(yint = yint,
              y = levels[yint],
              levels = levels,
              predint = predint,
              pred = levels[predint],
              altint = altint,
              altlab = levels[altint],
              PAC = PAC,
              figparams = farout$figparams,
              fig = farout$fig,
              farness = farout$farness,
              ofarness = farout$ofarness,
              classMS = classMS,
              lCurrent = lCurrent,
              lPred = lPred,
              lAlt = lAlt))
}


vcr.da.newdata <- function(Xnew, ynew = NULL, vcr.da.train.out) {
  #
  # Predicts class labels for new data by discriminant analysis,
  # using the output of vcr.da.train() on the training data.
  # For new data cases whose label in yintnew is non-missing,
  # additional output is produced for constructing graphical
  # displays.
  #
  # Arguments:
  #   Xnew             : data matrix of the new data, with the same
  #                      number of columns d as in the training data.
  #                      Missing values are not allowed.
  #   ynew             : factor with class membership of each new
  #                      case. Can be NA for some or all cases.
  #                      If NULL, is assumed to be NA everywhere.
  #   vcr.da.train.out : output of vcr.da.train() on the training
  #                      data. This also contains d, nlab, and levels.
  #
  # Returns:
  #   yintnew   : given labels as integer 1, 2, 3, ..., nlab.
  #               Can be NA.
  #   ynew      : labels if yintnew is available, else NA.
  #   levels    : levels of the response, from vcr.da.train.out
  #   predint   : predicted label as integer, always exists.
  #   pred      : predicted label of each object
  #   altint    : alternative label as integer if yintnew was given,
  #               else NA
  #   altlab    : alternative label if yintnew was given, else NA
  #   fig       : Mahalanobis distance of each object i to each class g
  #   classMS   : list with center and covariance matrix of each class,
  #               from vcr.da.train.out
  #   lCurrent  : log of density in object's given class, if yintnew
  #               was given, else NA
  #   lPred     : log of density in object's predicted class
  #   lAlt      : log of density in object's alternative class if
  #               yintnew was given, else NA
  #   PAC       : probability of alternative class for each object
  #               with non-NA yintnew
  #   figparams : (from training) parameters used to compute fig
  #   fig       : farness of each object i from each class g.
  #               Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any class.
  #               Always exists. If much above 1, we suspect i may be
  #               an outlier.
  #
  Xnew <- as.matrix(Xnew) # in case it is a data frame
  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  d <- length(vcr.da.train.out$classMS[[1]]$m) # = ncol(X)
  if (ncol(Xnew) != d) {
    stop(paste0("Xnew should have ", d,
                " columns, like the training data."))
  }
  ### We could also check whether the colnames are the same...
  if (sum(is.na(as.vector(Xnew))) > 0) {
    stop("The data matrix Xnew contains NA's.")
  }
  levels <- vcr.da.train.out$levels # same names as in training data
  nlab   <- length(levels) # number of classes
  #
  if (is.null(ynew)) ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = F, levels)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv   # INDiceS of cases we can Visualize,
  #                         # can even be empty, is OK.
  nvis    <- length(indsv)   # can be zero.
  yintnew <- rep(NA, n)       # needed to overwrite "new" levels.
  yintnew[indsv] <- lab2int(ynew[indsv])
  # Any "new" levels are now set to NA.
  yintv   <- yintnew[indsv]  # entries of yintnew that can be visualized
  #
  MD2s <- matrix(0, n, nlab) # for squared Mahalanobis distances
  mixf <- lmixf <- matrix(0, n, nlab) # mixture density and its log
  classMS <- vcr.da.train.out$classMS
  for (g in seq_len(nlab)) { # g=1
    M          <- classMS[[g]]$m
    S          <- classMS[[g]]$S
    prior      <- classMS[[g]]$prior
    MD2s[, g]  <- mahalanobis(Xnew, M, S)
    lmixf[, g] <- -0.5*(log(det(S)) + MD2s[, g]) + log(prior)
    mixf[, g]  <- exp(lmixf[, g])
  }
  # Compute predictions far all objects in 1, ..., n
  predint <- lPred <- rep(NA, n)
  lPred   <- apply(lmixf[, , drop = FALSE], 1, max)
  predint <- apply(lmixf[, , drop = FALSE], 1, which.max)
  #
  # For PAC:
  #
  lCurrent <- fCurr <- fAlt <- lAlt <- altint <- rep(NA, n)
  PAC <- rep(NA, n) # will remain NA for cases outside indsv
  # We can only compute PAC for cases with available yint.
  if (nvis > 0) { # if there are new data with labels
    yintv  <- yintnew[indsv] # yint of cases we will Visualize
    ayintv <- sort(unique(yintv)) # available yintv values
    # We cannot assume that ayintv = 1, ..., nlab for test data.
    for (g in ayintv) { # g=1 # g=3
      clinds <- indsv[which(yintv == g)] # indices in 1, ..., n
      others <- (seq_len(nlab))[-g]
      fCurr[clinds]    <- mixf[clinds, g]
      lCurrent[clinds] <- lmixf[clinds, g]
      fAlt[clinds]     <- apply(mixf[clinds, others, drop = FALSE],
                                1, max)
      lAlt[clinds]     <- apply(lmixf[clinds, others, drop = FALSE],
                                1, max)
      altint[clinds]   <- others[apply(lmixf[clinds, others, drop = FALSE],
                                       1, which.max)]
    }
    #
    # Compute PAC:
    #
    PAC[indsv] <- fAlt[indsv] / (fCurr[indsv] + fAlt[indsv])
  }
  #
  # Compute farness:
  #
  farout <- compFarness(type = "affine", testdata = T,
                        yint = yintnew,
                        nlab = nlab, X = NULL,
                        fig = sqrt(MD2s), d = d,
                        figparams = vcr.da.train.out$figparams)
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
              ofarness = farout$ofarness,
              classMS = classMS,
              lCurrent = lCurrent,
              lPred = lPred,
              lAlt = lAlt))
}
