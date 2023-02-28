vcr.svm.train <- function(X, y, svfit, ortho = FALSE) {
  #
  # Computes the quantities required for classmap().
  # Unfortunately X and y are not in the output of svm(),
  # so they must be put in explicitly here. Be careful that
  # X and yint are the same as what was used in svm() !
  #
  # Arguments:
  # X       : matrix of data coordinates, as used in svm().
  #           Missing values are not allowed.
  # y       : factor with the given (observed) class labels.
  #           It is crucial that X and y are EXACTLY the same
  #           as in the call to e1071:::svm().
  #           Since svm() allows y to contain NA's, we do too.
  #           If _all_ entries of y are NA we stop with a message,
  #           but in that situation svm() stops anyway.
  # svfit   : an object returned by svm(), called with EXACTLY
  #           the same X and y as above.
  # ortho   : TRUE will compute farness in the orthogonal
  #           complement of the vector beta given by svm.
  #           Is only possible for 2 classes, else there
  #           would be several beta vectors.
  #
  # Values:
  #   yint      : given class numbers of the cases, as it was input.
  #   levels    : labels (names) of the classes, as it was input.
  #   y         : given label of each case, as names.
  #   svfit     : as it was input, will be useful for new data.
  #   predint   : number of the predicted class, for each object.
  #   pred      : predicted label of each object.
  #   altint    : number of alternative class, i.e. the preferred
  #               non-self class.
  #   altlab    : alternative label of each object.
  #   PAC       : probability of alternative class for each object
  #   PCAfits   : PCA fits to each class, estimated from the
  #               training data but also useful for new data.
  #   figparams : parameters used in fig, can be used for new data.
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
  d <- ncol(X) # can be 1, is OK
  if (n < 2) stop(" X should have more than one row")
  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = TRUE)
  # If it did not stop: yint has length n, and its values are in
  # 1, ..., nlab without gaps. It is NA where y is: outside indsv.
  lab2int <- checked$lab2int # is a function
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y)
  indsv   <- checked$indsv
  nvis    <- length(indsv) # is >= 1 (else checkLabels has stopped)
  yintv   <- yint[indsv] # entries of yint that can be visualized
  #
  if (ortho == TRUE & nlab != 2) stop(
    " ortho=TRUE is only possible for 2 classes")
  if (is.null(svfit$compprob)) stop(" svfit$compprob is missing")
  if (svfit$compprob == FALSE) {
    stop(paste0("\n The object svfit needs to be obtained by ",
                "\n running e1071::svm with probability=TRUE"))
  }
  if (is.null(svfit$kernel)) {
    stop(" svfit$kernel is missing")
  } else {
    k <- svfit$kernel
  }
  if (!(k == 0 | k == "linear" | k == 1 | k == "polynomial" |
        k == 2 | k == "radial" | k == 3 | k == "sigmoid")) {
    stop(paste(" Unknown kernel type: ", k, sep = ""))
  }
  if (is.null(svfit$scaled)) {
    stop(" svfit$scaled is missing")
  } else {
    scaled <- svfit$scaled
  }
  if (is.null(svfit$degree)) {
    stop(" svfit$degree is missing")
  } else {
    degree <- svfit$degree
  }
  if (is.null(svfit$gamma)) {
    stop(" svfit$gamma is missing")
  } else {
    gamma <- svfit$gamma
  }
  if (is.null(svfit$coef0)) {
    stop(" svfit$coef0 is missing")
  } else {
    coef0 <- svfit$coef0
  }
  if (is.null(svfit$cost)) {
    stop(" svfit$cost is missing")
  } else {
    cost <- svfit$cost
  }
  if (is.null(svfit$fitted) | is.null(svfit$decision.values)) {
    # svfit is not an svm object, so we need to run svm() in
    # order to obtain both these outputs that we need.
    if (k == 0) {ktype <- "linear"}
    if (k == 1) {ktype <- "polynomial"}
    if (k == 2) {ktype <- "radial"}
    if (k == 3) {ktype <- "sigmoid"}
    svfit <- svm(X, y, scale = scaled, kernel = ktype,
                 degree = degree, gamma = gamma, coef0 = coef0,
                 cost = cost, probability = TRUE)
  }
  #
  # Compute predint, altint, and PAC
  #
  preds <- predict(object = svfit,
                   newdata = data.frame(X = X),
                   probability = TRUE)

  predint <- lab2int(preds) # for all cases, even those
  #                        # with missing y.
  altint <- PAC <- rep(NA, n)
  # altint and PAC will be NA outside on indsv, so we put:
  altintv <- PACv <- rep(NA, nvis)
  #
  if (nlab == 2) { # Binary classification. In this setting there
    # is only a single decision variable.
    altintv <- 3 - yintv # the alternative is always the other class
    #
    probs <- attr(preds, "probabilities")[indsv, ]
    probs <- probs[, order(lab2int(colnames(probs)))]
    for (i in seq_len(nvis)) PACv[i] <- probs[i, altintv[i]]
  }
  if (nlab > 2) { # More than 2 classes
    probs <- attr(preds, "probabilities")[indsv, ]
    probs <- probs[, order(lab2int(colnames(probs)))]
    # rowSums(probs) # all 1
    for (i in seq_len(nvis)) {
      ptrue      <- probs[i, yintv[i]]
      others     <- (seq_len(nlab))[-yintv[i]]
      palt       <- max(probs[i, others])
      altintv[i] <- others[which.max(probs[i, others])]
      PACv[i]    <- palt / (ptrue + palt)
    }
  } # ends nlab > 2
  altint[indsv] <- altintv # and stays NA outside indsv
  PAC[indsv]    <- PACv    # and stays NA outside indsv
  #
  # Compute farness
  #
  X <- X[indsv, ] # because this is what we will need to
  # keep for use with the new data
  Xloadings <- transfmat <- nbeta <- NULL
  if (k == 0) { # kernel was linear
    Xfar <- X
  } else {# kernel was nonlinear, so construct feature space:
    FVout <- makeFV(makeKernel(X, svfit = svfit), precS = 1e-12)
    Xfar  <- FVout$Xf # X in feature space
    FVout <- FVout$transfmat # to save space
    transfmat <- FVout; rm(FVout) # R's rename
    # We keep transfmat for new data.
  }

  if (k == 0) { # kernel was linear
    PCout <- truncPC(Xfar, ncomp = min(20, ncol(Xfar)),
                     center = FALSE, scale = FALSE)
    eigs  <- PCout$eigenvalues[seq_len(PCout$rank)]
    qkeep <- sum(eigs > 1e-8)
    Xloadings <- PCout$loadings[, seq_len(qkeep), drop = FALSE]
    Xfar <- Xfar %*% Xloadings
  } else {# kernel is nonlinear, so Xfar stems from makeFV
    qkeep <- min(20, ncol(Xfar))
    Xfar <- Xfar[, seq_len(qkeep)]
    attr(transfmat, "keep") <- qkeep
  }
  #
  if (ortho == TRUE) {
    # We would have stopped if nlab!=2, so nlab is 2.
    # In examples with spherical kernels (RBF, string)
    # the choice ortho=T made no visible difference.
    fitfs <- svm(Xfar, y, scale = FALSE, kernel = "linear",
                 cost = cost)
    beta <- coef(fitfs)[-1]
    beta <- matrix(beta, ncol = 1)

    # For orhogonal complement, first normalize beta:
    nbeta <- beta / sqrt(sum(beta * beta))
    # Project Xfar orthogonally to nbeta:
    Xfar <- Xfar - Xfar %*% nbeta %*% t(nbeta)
  }
  farout <- compFarness(type = "pca", testdata = FALSE, yint = yintv,
                        nlab = nlab, X = Xfar, fig = NULL, d = NULL,
                        figparams = NULL, PCAfits = NULL, keepPCA = TRUE)
  fig <- matrix(rep(0, n * nlab), ncol = nlab)
  fig[indsv,] <- farout$fig
  farness <- ofarness <- rep(NA, n)
  farness[indsv] <- farout$farness
  ofarness[indsv] <- farout$ofarness
  # Add all intermediate farness objects to figparams:
  figparams <- farout$figparams
  figparams$ortho <- ortho
  figparams$nbeta <- nbeta
  figparams$Xloadings <- Xloadings
  figparams$transfmat <- transfmat
  figparams$PCAfits <- farout$PCAfits
  return(list(yint = yint,
              y = levels[yint],
              levels = levels,
              predint = predint,
              pred = levels[predint],
              altint = altint,
              altlab = levels[altint],
              PAC = PAC,
              figparams = figparams,
              fig = fig,
              farness = farness,
              ofarness = ofarness,
              svfit = svfit,
              X = X))
}


vcr.svm.newdata <- function(Xnew, ynew = NULL, vcr.svm.train.out) {
  #
  # Predicts class labels for new data by SVM, using the output
  # of vcr.svm.train() on the training data.
  # For new data cases whose label in yintnew is non-missing,
  # additional output is produced for constructing graphical
  # displays.
  #
  # Arguments:
  #   Xnew              : data matrix of the new data, with the same
  #                       number of columns d as in the training data.
  #                       Missing values in Xnew are not allowed.
  #   ynew              : factor with class membership of each new
  #                       case. Can be NA for some or all cases.
  #                       If NULL, is assumed to be NA everywhere.
  #   vcr.svm.train.out : output of vcr.svm.train() on the training
  #                       data. This also contains d, nlab, and levels.
  #
  # Values:
  #   yintnew   : given labels as integer 1, 2, 3, ..., nlab.
  #               Can be NA.
  #   ynew      : labels if yintnew is available, else NA.
  #   levels    : levels of the response, from vcr.svm.train.out
  #   predint   : predicted class number, always exists.
  #   pred      : predicted label of each object
  #   altint    : alternative label as integer if yintnew was given,
  #               else NA
  #   altlab    : alternative label if yintnew was given, else NA
  #   PAC       : probability of alternative class for each object
  #               with non-NA yintnew
  #   PCAfits   : PCA fits to each class, from vcr.svm.train.out .
  #               This is required, else an error message is given.
  #   figparams : (from training data) parameters used for fig.
  #   fig       : farness of each object in the new data from each
  #               class g. Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any class.
  #               Always exists. If much above 1, we suspect i may be
  #               an outlier.
  #
  Xnew <- as.matrix(Xnew) # in case it is a data frame
  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  if (is.null(vcr.svm.train.out$figparams$PCAfits)) {
    stop(paste0("\nFor test data, the farness computation ",
                "requires the trained object", "\nto contain ",
                "the value figparams$PCAfits."))
  }
  PCAfits <- vcr.svm.train.out$figparams$PCAfits
  if (sum(is.na(as.vector(Xnew))) > 0) stop(
    "The data matrix Xnew contains NA's.")
  levels <- vcr.svm.train.out$levels # same as in training data
  nlab   <- length(levels) # number of classes
  #
  if (is.null(ynew)) ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = FALSE, levels)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv   # INDiceS of cases we can Visualize,
  #                         # can even be empty, is OK.
  nvis    <- length(indsv)   # can be zero.
  yintnew <- rep(NA, n)       # needed to overwrite "new" levels.
  yintnew[indsv] <- lab2int(ynew[indsv])
  # Any "new" levels are now set to NA.
  yintv   <- yintnew[indsv]  # entries of yintnew that can be visualized
  #
  # Make predictions for new data:
  #
  svfit <- vcr.svm.train.out$svfit # trained model
  k     <- svfit$kernel
  preds <- predict(object = svfit,
                   newdata = data.frame(X = Xnew),
                   decision.values = TRUE,
                   probability = TRUE)
  predint <- lab2int(preds)
  #
  # Compute altint and PAC
  #
  altint <- PAC <- rep(NA, n) # initializations
  # altint and PAC will be NA outside of indsv.
  if (nvis > 0) { # If there are cases to visualize
    altintv <- PACv <- rep(NA, nvis)
    #
    if (nlab == 2) { # Binary classification.
      altintv <- 3 - yintv # the alternative is always the other class
      #
      probs <- attr(preds, "probabilities")[indsv, ]
      probs <- probs[, order(lab2int(colnames(probs)))]
      for (i in seq_len(nvis)) PACv[i] <- probs[i, altintv[i]]
    }
    if (nlab > 2) { # More than 2 classes
      #
      probs <- attr(preds, "probabilities")[indsv, ]
      probs <- probs[, order(lab2int(colnames(probs)))]
      # nvis by nlab matrix with estimated class probabilities
      # rowSums(probs) # all 1
      for (i in seq_len(nvis)) {
        ptrue      <- probs[i, yintv[i]]
        others     <- (seq_len(nlab))[-yintv[i]]
        palt       <- max(probs[i, others])
        altintv[i] <- others[which.max(probs[i, others])]
        PACv[i]    <- palt / (ptrue + palt)
      }
    }
    PAC[indsv] <- PACv # in the other cases it remains NA
    altint[indsv] <- altintv # also here.
  }
  #
  # Compute farness
  #
  Zfar <- Xnew; rm(Xnew) # R's rename to save space
  if (k == 0) { # kernel was linear
    Xloadings <- vcr.svm.train.out$figparams$Xloadings # from training data
    if (!is.null(Xloadings)) { # (else Zfar remains unchanged)
      # We carried out an initial PCA on the training data, so
      # we have to replicate this here, else our PCAfits are
      # no longer in the same coordinate system.
      Zfar <- Zfar %*% Xloadings # to match the training data
    }
  } else {# kernel was not linear
    if (is.null(vcr.svm.train.out$figparams$transfmat)) {
      stop("vcr.svm.train.out$figparams$transfmat is missing")
    }
    # We had to store X in vcr.svm.train.out to be able to do:
    Zfar <- makeFV(makeKernel(Zfar, vcr.svm.train.out$X,
                              svfit = vcr.svm.train.out$svfit),
                   vcr.svm.train.out$figparams$transfmat)$Xf
  }
  # Zfar lives in the Xfar space from the training data.
  # dim(Zfar)
  if (vcr.svm.train.out$figparams$ortho == TRUE) {
    # We would have stopped if nlab!=2, so nlab is 2.
    nbeta <- vcr.svm.train.out$figparams$nbeta
    # nbeta is in the Xfar space too.
    # Project Zfar orthogonally to nbeta:
    Zfar <- Zfar - Zfar %*% nbeta %*% t(nbeta)
    # (Zfar%*%nbeta) # all 0: Zfar is on a hyperplane through 0.
  }
  #  dim(Zfar)
  farout <- compFarness(type = "pca", testdata = TRUE, yint = yintnew,
                        nlab = nlab, X = Zfar, fig = NULL,
                        d = NULL,
                        figparams = vcr.svm.train.out$figparams,
                        PCAfits = PCAfits, keepPCA = FALSE)
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
