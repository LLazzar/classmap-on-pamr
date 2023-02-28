
vcr.neural.train <- function(X, y, probs, estmethod = meancov) {
  #
  # Using the outputs of a neural network for classification
  # applied to the training data, this function prepares for
  # graphical displays.
  #
  # Arguments:
  #   X         : the coordinates of the n objects of the training
  #               data, in the layer chosen by the user. Missing
  #               values are not allowed.
  #   y         : factor with the given class labels of the objects.
  #               Make sure that the levels are in the same order as
  #               used in the neural net, i.e. the columns of its
  #               binary "once-hot-encoded" response vectors.
  #   probs     : posterior probabilities obtained by the neural
  #               net, e.g. in keras. For each case (row of X),
  #               the classes have probabilities that add up to 1.
  #               Each row of the matrix probs contains these
  #               probabilities. The columns of probs must be in
  #               the same order as the levels of y.
  #   estmethod : function for location and covariance estimation.
  #               Should return a list with $m and $S. Can be meancov
  #               (classical mean and covariance matrix) or DetMCD.
  #               If one or more classes have a singular covariance
  #               matrix, the function automatically switches to
  #               the PCA-based farness used in vcr.svm.train().
  #
  # Returns:
  #   yint      : given labels as integer 1, 2, 3, ...
  #   y         : given labels
  #   levels    : levels of y
  #   predint   : predicted class number. For each case this
  #               is the class with the highest probability,
  #               that is, which.max(probs[i,]).
  #   pred      : predicted label of each object.
  #   altint    : alternative class as integer, i.e. the non-given
  #               class number with the highest mixture density.
  #   altlab    : alternative class of each object.
  #   classMS   : list with center and covariance matrix of each class
  #   PAC       : probability of alternative class for each object
  #   PCAfits   : if not NULL, PCA fits to each class, estimated from
  #               the training data but also useful for new data.
  #   figparams : parameters for computing fig, can be used for
  #               new data.
  #   fig       : distance of each object i from each class g.
  #   farness   : farness of each object to its given class.
  #   ofarness  : For each object i, its lowest farness to any
  #               class, including its own. Always exists.
  #
  keepPCA <- TRUE; prec <- 1e-10

  X <- as.matrix(X) # in case it is a data frame
  if (nrow(X) == 1) X <- t(X)
  if (sum(is.na(as.vector(X))) > 0) {
    stop("The coordinate matrix X has NA's.")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2) stop("The training data should have more than one case.")
  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = TRUE)
  # If it did not stop: yint has length n, and its values are in
  # 1, ..., nlab without gaps. It is NA where y is: outside indsv.
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y)
  yintv   <- yint[indsv]
  classSizes <- rep(0, nlab)
  for (g in seq_len(nlab)) {classSizes[g] <- sum(yintv == g)}
  # classSizes
  #
  # Check matrix of posterior probabilities:
  #
  probs <- as.matrix(probs)
  if (length(dim(probs)) != 2) stop("probs should be a matrix.")
  if (nrow(probs) != n) stop(paste0(
    "The matrix probs should have ", n, " rows"))
  if (ncol(probs) != nlab) stop(paste0(
    "The matrix probs should have ", nlab, " columns"))
  if (any(is.na(probs))) stop("probs should not have any NA's.")
  #
  # Compute prediction for all objects in the training data:
  #
  predint <- apply(probs[, , drop = FALSE], 1, which.max)
  #
  # Compute ptrue and palt for all objects with available y:
  #
  ptrue <- palt <- altint <- PAC <- rep(NA, n)
  for (g in seq_len(nlab)) { # g=1
    clinds <- indsv[which(yintv == g)] # indices in 1, ..., n
    others <- (seq_len(nlab))[-g] # alternative classes
    ptrue[clinds]  <- probs[clinds, g]
    palt[clinds]   <- apply(probs[clinds, others, drop = FALSE], 1, max)
    altint[clinds] <- others[apply(probs[clinds, others, drop = FALSE],
                                   1, which.max)]
  }
  #
  # Compute PAC:
  #
  PAC[indsv] <- palt[indsv] / (ptrue[indsv] + palt[indsv])
  # (PAC and altint stay NA outside indsv)
  #
  # Compute farness:
  #
  computeMD <- TRUE
  tinyclasses <- which(classSizes < (d + 1))
  if (length(tinyclasses) > 0) {
    wnq(paste0("\nThere is at least one class with fewer",
               " than the ", d + 1, "\nmembers required to compute a ",
               "nonsingular covariance matrix, ", "\nso farness ",
               "will be computed from PCA."))
    computeMD <- FALSE
  } else {
    classMS <- list()
    # Compute the center and covariance matrix of each class:
    iCN <- rep(NA, nlab) # for numerical conditions
    for (g in seq_len(nlab)) {
      clinds       <- indsv[which(yintv == g)]
      # class indices in 1, ..., n
      class.out    <- estmethod(X[clinds, ]) # contains mu and Sigma
      iCN[g]       <- numericalCond(class.out$S)
      classMS[[g]] <- class.out # now has $m, $S
    }
    tinyCN <- which(iCN < prec) # also ran with 1e-10 to test PCA
    if (length(tinyCN) > 0) {
      wnq(paste0("\nThere is at least one class with ",
                 "(near-) singular covariance matrix, ",
                 "\nso farness will be computed from PCA."))
      computeMD <- FALSE
    }
  }
  if (computeMD == TRUE) {
    # compute MD2 of each object to all classes:
    initfig <- matrix(0, n, nlab) # for initial farness
    for (g in seq_len(nlab)) { # g=1
      M <- classMS[[g]]$m
      S <- classMS[[g]]$S # now S is nonsingular
      initfig[, g]  <- sqrt(mahalanobis(X, M, S))
      # Computed for all objects, even when response y is NA.
    }
    farout <- compFarness(type = "affine", testdata = FALSE, yint = yint,
                          nlab = nlab, X = NULL, fig = initfig,
                          d = d, figparams = NULL)
  } else {
    classMS <- NULL
    farout <- compFarness(type = "pca", testdata = FALSE, yint = yint,
                          nlab = nlab, X = X, fig = NULL,
                          d = NULL, figparams = NULL,
                          PCAfits = NULL, keepPCA = keepPCA)
  }

  figparams <- farout$figparams
  figparams$ncolX <- d
  figparams$computeMD <- computeMD
  figparams$classMS <- classMS
  figparams$PCAfits <- farout$PCAfits
  return(list(X = X,
              yint = yint,
              y = levels[yint],
              levels = levels,
              predint = predint,
              pred = levels[predint],
              altint = altint,
              altlab = levels[altint],
              PAC = PAC,
              figparams = figparams,
              fig = farout$fig,
              farness = farout$farness,
              ofarness = farout$ofarness))
}


vcr.neural.newdata <- function(Xnew, ynew = NULL, probs,
                               vcr.neural.train.out){
  #
  # Prepares graphical display of new data fitted by a neural
  # net that was modeled on the training data, using the output
  # of vcr.neural.train() on the training data.
  #
  # Arguments:
  #   Xnew                 : data matrix of the new data, with the
  #                          same number of columns d as in the
  #                          training data.
  #                          Missing values are not allowed.
  #   ynew                 : factor with class membership of each new
  #                          case. Can be NA for some or all cases.
  #                          If NULL, is assumed to be NA everywhere.
  #   probs                : posterior probabilities obtained by
  #                          running the neural net on the new data.
  #   vcr.neural.train.out : output of vcr.neural.train() on the
  #                          training data.
  #
  # Returns:
  #   yintnew   : given labels as integers 1, 2, 3, ..., nlab.
  #               Can have NA's.
  #   ynew      : labels if yintnew is available, else NA.
  #   levels    : levels of the response, from vcr.neural.train.out
  #   predint   : predicted label as integer, always exists.
  #   pred      : predicted label of each object
  #   altint    : alternative label as integer if yintnew was given,
  #               else NA.
  #   altlab    : alternative label if yintnew was given, else NA.
  #   classMS   : list with center and covariance matrix of each class,
  #               from vcr.neural.train.out
  #   PAC       : probability of alternative class for each object
  #               with non-NA yintnew.
  #   figparams : (from training) parameters used to compute fig
  #   fig       : farness of each object i from each class g.
  #               Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any
  #               class, including its own. Always exists.
  #
  Xnew <- as.matrix(Xnew) # in case it is a data frame
  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  d <- vcr.neural.train.out$figparams$ncolX
  if (ncol(Xnew) != d) {
    stop(paste0("Xnew should have ", d,
                " columns, like the training data."))
  }
  if (sum(is.na(as.vector(Xnew))) > 0) {
    stop("The coordinate matrix Xnew contains NA's.")
  }
  levels <- vcr.neural.train.out$levels # as in training data
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
  #
  probs <- as.matrix(probs)
  if (length(dim(probs)) != 2) stop("probs should be a matrix.")
  if (ncol(probs) == 1) probs <- t(probs) # if new data is 1 object
  if (ncol(probs) != nlab) stop(paste0(
    "The matrix probs should have ", nlab, " columns"))
  if (any(is.na(probs))) stop("probs should not have any NA's.")
  #
  # Compute prediction for all objects in the new data:
  #
  predint <- apply(probs[, , drop = FALSE], 1, which.max)
  #
  # Compute PAC for all objects with available y:
  #
  ptrue <- palt <- altint <- PAC <- rep(NA, n)
  if (nvis > 0) { # if there are new data with labels
    yintv  <- yintnew[indsv] # yint of cases we will Visualize
    ayintv <- sort(unique(yintv)) # available yintv values
    # ayintv can be a proper subset of 1, ..., nlab for new data.
    for (g in ayintv) {
      clinds <- indsv[which(yintv == g)] # indices in 1, ..., n
      others <- (seq_len(nlab))[-g] # non-self classes
      ptrue[clinds]  <- probs[clinds, g]
      palt[clinds]   <- apply(probs[clinds, others, drop  = FALSE], 1, max)
      altint[clinds] <- others[apply(probs[clinds, others, drop = FALSE],
                                     1, which.max)]
    }
    PAC[indsv] <- palt[indsv] / (ptrue[indsv] + palt[indsv])
    # (PAC and altint stay NA outside indsv)
  }
  #
  # Compute farness:
  #
  if (vcr.neural.train.out$figparams$computeMD == TRUE) {
    # Compute MD2S
    initfig <- matrix(0, n, nlab) # for Mahalanobis distances
    classMS <- vcr.neural.train.out$figparams$classMS
    for (g in seq_len(nlab)) {
      M <- classMS[[g]]$m
      S <- classMS[[g]]$S
      initfig[, g] <- sqrt(mahalanobis(Xnew, M, S))
      # Computed for all objects, even when ynew is NA.
    }
    farout <- compFarness(type = "affine", testdata = TRUE,
                          yint = yintnew, nlab = nlab, X = NULL,
                          fig = initfig, d = d,
                          figparams = vcr.neural.train.out$figparams)
  } else {
    farout <- compFarness(type = "pca", testdata = TRUE, yint = yintnew,
                          nlab = nlab, X = Xnew, fig = NULL, d = NULL,
                          figparams = vcr.neural.train.out$figparams,
                          PCAfits = vcr.neural.train.out$figparams$PCAfits,
                          keepPCA = FALSE)
  }
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
