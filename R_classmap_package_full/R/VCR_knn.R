
vcr.knn.train <- function(X, y, k){
  #
  # Carries out a k-nearest neighbor classification on the
  # training data. Each case is predicted from the labels of
  # all cases except its own.
  # Various additional output is produced for the purpose
  # of constructing graphical displays.
  #
  # X        : This can be a rectangular matrix or data frame
  #            of (already standardized) measurements, or a
  #            dist object obtained from dist() or daisy().
  #            Missing values are not allowed.
  # y        : factor with the given (observed) class labels.
  #            There need to be non-missing y in order to be
  #            able to train the classifier.
  # k        : the number of nearest neighbors used. It can
  #            be selected by running cross-validation using
  #            a different package.
  #
  # Values include:
  #   yint      : given labels as integer 1, 2, 3, ...
  #   y         : given labels
  #   levels    : levels of y
  #   predint   : predicted class number. For each case this
  #               is the class with the highest mixture density.
  #   pred      : predicted label of each object.
  #   altint    : alternative class as integer, i.e. the non-given
  #               class number with the highest mixture density.
  #   altlab    : alternative class of each object.
  #   PAC       : probability of alternative class for each object
  #   fig       : distance of each object i from each class g.
  #   farness   : farness of each object to its given class.
  #   ofarness  : For each object i, its lowest farness to any class.
  #               If much above 1, we suspect i may be an outlier.

  # Auxiliary function for kNN:
  #
  predDis <- function(nlab, labsi, dissi) {
    # Computes the kNN prediction and its average distance.
    # When the highest frequency is tied, the predicted class
    # is the one with lowest average distance from i.
    freqz <- rep(0, times = nlab)
    for (j in seq_len(nlab)) { # freqz = frequencies of labels
      freqz[j] <- sum(labsi == j, na.rm = TRUE)
    }
    ranking <- order(freqz, decreasing = TRUE)
    # positions of most to least frequent labels
    sfreq <- freqz[ranking] # frequencies in that order
    predi <- NA
    if (sfreq[1] > sfreq[2]) { # frequency has only 1 maximizer
      predi <- ranking[1] # predicted class
      averdis <- sum(dissi[labsi == predi]) / sfreq[1]
      seconddis <- NA
    } else {# frequency has several maximizers
      multipreds <- ranking[which(sfreq == sfreq[1])]
      multi <- length(multipreds); multi
      avdis <- rep(0, times = multi)
      for (j in seq_len(multi)) { # adds at each tied solution
        avdis[j] <- sum(dissi[labsi == multipreds[j]]) / sfreq[1]
      }
      bestInMulti <- which.min(avdis) # avdis minimizer
      averdis <- avdis[bestInMulti] # smallest average distance
      # compute second smallest average distance:
      seconddis <- min(avdis[-bestInMulti])
      predi <- ranking[bestInMulti] # predicted class
    }
    freq <- sfreq[1]
    return(list(predi = predi, averdis = averdis, freq = freq,
                seconddis = seconddis))
  }

  # Here the main function starts:
  if (inherits(X, "dist")) { # input is a dist object
    n <- cluster::sizeDiss(X); n # 150
    # turn it into a symmetric square matrix:
    X <- as.matrix(X); dim(X)
    dismat <- X; rm(X) # R's renaming, to save space
    if (sum(is.na(as.vector(dismat))) > 0) {
      stop("The input dissimilarity matrix X has NA's.")
    }
    X <- NULL # we don't need to return this
  } else {
    # Assumes X is already standardized
    X <- as.matrix(X) # in case it is a data frame
    if (nrow(X) == 1) X <- t(X)
    if (sum(is.na(as.vector(X))) > 0) {
      stop(paste0("The data matrix X has NA's. You could either",
                  " remove them", "\nor make a dissimilarity",
                  " matrix without NA's, e.g. with daisy()."))
    }
    n <- nrow(X)
    dismat <- as.matrix(dist(X, method = "euclidean"))
    # If you want a different metric, call dist() or daisy() first
    # and then input the result into vcr.knn.train instead of
    # the coordinate matrix.
  }
  if (n < 2) stop("The training data should have more than one case.")
  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = TRUE)
  # If it did not stop: yint has length n, and its values are in
  # 1, ..., nlab without gaps. It is NA where y is: outside indsv.
  lab2int <- checked$lab2int
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y)
  indsv   <- checked$indsv # cases in indsv can be visualized
  nvis    <- length(indsv)
  yintv   <- yint[indsv]
  #
  # Only neighbors with known class are useful for training:
  # dim(dismat)
  dismat <- dismat[, indsv] # has the right colnames
  # The k nearest neighbors of each i, as indices in 1, ..., n are:
  sortNgb <- t(apply(dismat, 1, order))[, -1] # excludes i itself,
  # this is the usual convention for knn on training data.
  # To make the entries case indices in 1, ..., n we do:
  sortNgb <- apply(sortNgb, 2, function(x) indsv[x])
  dismat <- t(apply(dismat, 1, sort))[, -1] # excludes d(i, i)=0
  sortDis <- dismat; rm(dismat) # R's renaming, to save space

  #
  # Make neighborhoods:
  #
  # First take care of situations where the k-th distance
  # is tied with (at least) the (k+1)-th. For such an object
  # we include all those tied distances.
  prec <- 1e-12
  ktrues <- rep(NA, times = n)
  for (i in seq_len(n)) { # i=1
    rowdis <- sortDis[i, ]
    ktrues[i] <- length(which(rowdis < rowdis[k] + prec))
  }
  #
  # Make predictions:
  #
  preds  <- matrix(rep(NA, times = 2 * n), ncol = 2)
  freqs  <- matrix(rep(0, times = 3 * n), ncol = 3) # no NA's, OK
  avedis <- matrix(rep(NA, times = 3 * n), ncol = 3)
  #
  for (i in seq_len(n)) { # loop over cases to predict # i=7
    #
    # Compute statistics for the true label, if it exists:
    #
    ktrue <- ktrues[i]
    labsi <- yint[sortNgb[i, seq_len(ktrue)]]
    # labels of neighbors exist
    dissi <- sortDis[i, seq_len(ktrue)] # dissimilarities from i
    if (!is.na(yint[i])) {
      freqs[i, 3] <- sum(labsi == yint[i])
      # Note that freqs[i, 3] can be zero, but not NA
      if (freqs[i, 3] > 0) {
        avedis[i, 3] <- sum(dissi[labsi == yint[i]]) / freqs[i, 3]
      }
    }
    #
    # Compute predicted label:
    #
    out1       <- predDis(nlab, labsi, dissi)
    # calls the auxiliary function predDis
    preds[i, 1] <- out1$predi
    freqs[i, 1] <- out1$freq
    if (!is.na(out1$averdis)) {
      avedis[i, 1] <- out1$averdis
    } else {
      avedis[i, 1] <- min(sortDis[i, yint[sortNgb[i, ]] == out1$predi])
    }
    #
    # Compute second best prediction if there is one:
    #
    if (out1$freq < ktrues[i]) { # else nothing to do
      remNgb      <- which(labsi != out1$predi) # remaining neighbors
      labsi       <- labsi[remNgb]
      dissi       <- dissi[remNgb]
      out2        <- predDis(nlab, labsi, dissi)
      preds[i, 2]  <- out2$predi   # was initialized as NA
      freqs[i, 2]  <- out2$freq    # was initialized as 0
      avedis[i, 2] <- out2$averdis # was initialized as NA
    }
  } # ends loop over cases to predict
  #
  counts <- freqs[, c(3, 1, 2)]
  # given class, predicted class, alternative class
  predint <- preds[, 1] # integer predictions for all cases,
  #                   # even those with missing y.
  freq1v <- freqs[indsv, 1] # highest frequency (= of prediction)
  freq2v <- freqs[indsv, 2] # second highest frequency
  freq3v <- freqs[indsv, 3] # frequency of self
  pred1v <- preds[indsv, 1] # prediction
  pred2v <- preds[indsv, 2] # best non-prediction
  #
  # Compute PAC:
  #
  PAC <- altint <- rep(NA, n)
  freqalt <- freq1v
  altintv <- pred1v # preliminary
  matched <- which(pred1v == yintv) # prediction == selfclass
  freqalt[matched] <- freq2v[matched]
  altintv[matched] <- pred2v[matched]
  PACv <- freqalt / (freq3v + freqalt)
  PAC[indsv] <- PACv
  altint[indsv] <- altintv
  #
  ### Compute initial fig
  #
  fig <- matrix(rep(NA, nvis * nlab), ncol = nlab)
  # compute fig[i, g] of object i to class g:
  for (i in seq_len(nvis)) { # loop over all cases in indsv
    for (g in seq_len(nlab)) { # loop over classes
      ngbg <- which(yint[sortNgb[indsv[i], ]] == g)
      if (length(ngbg) > k) {
        ngbg <- ngbg[seq_len(k)]
      }
      fig[i, g] <- median(sortDis[indsv[i], ngbg])
    }
  }
  #
  # Compute farness:
  #
  farout <- compFarness(type = "knn", testdata = FALSE, yint = yintv,
                        nlab = nlab, X = NULL, fig = fig, d = NULL,
                        figparams = NULL)
  fig <- matrix(rep(0, n * nlab), ncol = nlab)
  fig[indsv, ] <- farout$fig
  farness <- ofarness <- rep(NA, n)
  farness[indsv] <- farout$farness
  ofarness[indsv] <- farout$ofarness
  if (nlab == 2) altint <- 3 - yint # for consistency with svm etc.
  # This way altint exists, even if it has frequency zero.
  return(list(yint = yint,
              y = levels[yint],
              levels = levels,
              predint = predint,
              pred = levels[predint],
              altint = altint,
              altlab = levels[altint],
              PAC = PAC,
              figparams = farout$figparams,
              fig = fig,
              farness = farness,
              ofarness = ofarness,
              k = k,
              ktrues = ktrues,
              counts = counts,
              X = X))
}

vcr.knn.newdata <- function(Xnew, ynew = NULL,
                            vcr.knn.train.out, LOO = FALSE) {
  #
  # Predicts class labels for new data by k nearest neighbors,
  # using the output of vcr.knn.train() on the training data.
  # For new data with available labels in ynew, additional output
  # is produced for constructing graphical displays.
  #
  #   Xnew              : If the training data was a matrix of
  #                       coordinates, Xnew must be such a matrix
  #                       with the same number of columns.
  #                       If the training data was a set of
  #                       dissimilarities, Xnew must be a rectangular
  #                       matrix of dissimilarities, with each row
  #                       containing the dissmilarities of a new
  #                       case to all training cases.
  #                       Missing values are not allowed.
  #   ynew              : factor with class membership of each new
  #                       case. Can be NA for some or all cases.
  #                       If NULL, is assumed to be NA everywhere.
  #   vcr.knn.train.out : output of vcr.knn.train() on the training
  #                       data. This contains nlab, levels and k.
  #   LOO               : leave one out. Only used when testing this
  #                       function on a subset of the training data.
  #                       Default is LOO=FALSE.
  #
  # Values include:
  #   yintnew   : given labels as integer 1, 2, 3, ..., nlab.
  #               Can be NA.
  #   ynew      : labels if yintnew is available, else NA.
  #   levels    : levels of the response, from vcr.da.train.out
  #   predint   : predicted class number, always exists.
  #   pred      : predicted label of each object
  #   altint    : alternative label as integer if yintnew was given,
  #               else NA
  #   altlab    : alternative label if yintnew was given, else NA
  #   PAC       : probability of alternative class for each object
  #               with non-NA yintnew
  #   fig       : distance of each object i from each class g.
  #               Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any class.
  #               Always exists. If much above 1, we suspect i may be
  #               an outlier.

  # Auxiliary function for kNN:
  #
  predDis <- function(nlab, labsi, dissi) {
    # Computes kNN prediction and its average distance.
    # This works even if the data contains < nlab labels.
    freqz <- rep(0, times = nlab)
    for (j in seq_len(nlab)) { # freqz = frequencies of labels
      freqz[j] <- sum(labsi == j, na.rm = TRUE)
    }
    ranking <- order(freqz, decreasing = TRUE)
    # positions of most to least frequent labels
    sfreq <- freqz[ranking] # frequencies in that order
    predi <- NA
    if (sfreq[1] > sfreq[2]) { # frequency has only 1 maximizer
      predi <- ranking[1] # predicted class
      averdis <- sum(dissi[labsi == predi]) / sfreq[1]
      seconddis <- NA
    } else {# frequency has several maximizers
      multipreds <- ranking[which(sfreq == sfreq[1])]
      multi <- length(multipreds); multi
      avdis <- rep(0, times = multi)
      for (j in seq_len(multi)) { # adds at each tied solution
        avdis[j] <- sum(dissi[labsi == multipreds[j]]) / sfreq[1]
      }
      bestInMulti <- which.min(avdis) # avdis minimizer
      averdis <- avdis[bestInMulti] # smallest average distance
      # compute second smallest average distance:
      seconddis <- min(avdis[-bestInMulti])
      predi <- ranking[bestInMulti] # predicted class
    }
    freq <- sfreq[1]
    return(list(predi = predi, averdis = averdis, freq = freq,
                seconddis = seconddis))
  }

  # Here the main function starts:
  Xnew <- as.matrix(Xnew)
  if (is.null(vcr.knn.train.out$X)) {
    # The training data were given as a dissimilarity matrix.
    # Therefore Xnew must be a dissimilarity matrix too.
    if (is.vector(Xnew)) {Xnew <- matrix(Xnew, nrow = 1)}
    Xnew <- as.matrix(Xnew)
    if (sum(is.na(as.vector(Xnew))) > 0) {
      stop("The dissimilarity matrix Xnew contains NA's.")
    }
    dismat <- Xnew; rm(Xnew) # R's version of rename, to save space.
    ntrain <- length(vcr.knn.train.out$yint)
    if (ncol(dismat) != ntrain) {
      stop(paste0("\nThe dissimilarity matrix Xnew should have ", ntrain,
                  "\ncolumns, as many as there were training cases."))
    }
    n <- nrow(dismat) # number of new cases
  } else {
    # The training data was a matrix of coordinates.
    Xnew <- as.matrix(Xnew) # in case it is a data frame
    if (nrow(Xnew) == 1) Xnew <- t(Xnew)
    d <- ncol(vcr.knn.train.out$X)
    if (ncol(Xnew) != d) {
      stop(paste0("Xnew should have ", d,
                  " columns, like the training data."))
    }
    # Should check column names here (also in da)
    # first vcr.da.train should output $colnames
    if (sum(is.na(as.vector(Xnew))) > 0) {
      stop("The data matrix Xnew contains NA's.")
    }
    ntrain <- nrow(vcr.knn.train.out$X)
    n      <- nrow(Xnew) # number of new cases
    Xall   <- rbind(vcr.knn.train.out$X, Xnew) # has ntrain + n rows
    indn   <- (ntrain + 1):(ntrain + n) # range of new rows
    dismat <- as.matrix(dist(Xall, method = "euclidean")
    )[indn, seq_len(ntrain), drop = FALSE]
  }
  levels <- vcr.knn.train.out$levels # same names as in training data
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
  yint <- vcr.knn.train.out$yint # class membership in training data
  k    <- vcr.knn.train.out$k
  # For each new case, find its neighbors in the training data:
  sortNgb <- t(apply(dismat, 1, order))
  if (LOO) sortNgb <- sortNgb[, -1] # leave one out, only for testing
  colnames(sortNgb) <- NULL # colnames were irrelevant
  # Obtain their dissimilarities:
  dismat  <- t(apply(dismat, 1, sort)) # to save space
  if (LOO) dismat <- dismat[, -1] # leave one out, only for testing
  sortDis <- dismat; rm(dismat) # to save space
  #
  # Make neighborhoods:
  #
  # First take care of situations where the k-th distance
  # is tied with (at least) the (k+1)-th. For such an object
  # we include all those tied distances.
  prec <- 1e-12
  ktrues <- rep(NA, n)
  for (i in seq_len(n)) { # i=1
    rowdis <- sortDis[i, ]
    ktrues[i] <- length(which(rowdis < rowdis[k] + prec))
  }
  ktrues
  #
  # Make predictions:
  #
  preds  <- matrix(rep(NA, times = 2 * n), ncol = 2)
  freqs  <- matrix(rep(0, times = 3 * n), ncol = 3)
  avedis <- matrix(rep(NA, times = 3 * n), ncol = 3)
  #
  for (i in seq_len(n)) { # loop over cases to predict
    #
    # Compute statistics for the true label, if it exists:
    #
    ktrue <- ktrues[i]
    labsi <- yint[sortNgb[i, seq_len(ktrue)]]
    # labels in training data exist
    dissi <- sortDis[i, seq_len(ktrue)]
    # dissim. between i and training data
    if (is.na(yintnew[i])) { # if i has no given label
      freqs[i, 3] <- 0
    } else {# if i does have a given label
      freqs[i, 3] <- sum(labsi == yintnew[i])
    }
    if (freqs[i, 3] > 0) {
      avedis[i, 3] <- sum(dissi[labsi == yintnew[i]]) / freqs[i, 3]
    } else {
      avedis[i, 3] <- NA
    }
    #
    # Compute predicted label:
    #
    out1       <- predDis(nlab, labsi, dissi)
    preds[i, 1] <- out1$predi
    freqs[i, 1] <- out1$freq
    if (!is.na(out1$averdis)) {
      avedis[i, 1] <- out1$averdis } else {
        avedis[i, 1] <- min(sortDis[i, yint[sortNgb[i, ]] == out1$predi]) }
    #
    # Compute second best prediction if there is one:
    #
    if (out1$freq < ktrues[i]) { # else nothing to do
      remNgb      <- which(labsi != out1$predi) # remaining neighbors
      labsi       <- labsi[remNgb]
      dissi       <- dissi[remNgb]
      out2        <- predDis(nlab, labsi, dissi)
      preds[i, 2]  <- out2$predi   # was initialized as NA
      freqs[i, 2]  <- out2$freq    # was initialized as 0
      avedis[i, 2] <- out2$averdis # was initialized as NA
    }
  } # ends loop over cases to predict
  #
  counts <- freqs[, c(3, 1, 2)]
  # given class, predicted class, alternative class
  predint <- preds[, 1] # integer predictions for all cases,
  #                   # even those with missing ynew.
  #
  PAC <- altint <- rep(NA, n)
  if (nvis > 0) { # there are new data with labels
    # We can only compute PAC etc. for cases with available ynew.
    freq1v <- freqs[indsv, 1]
    freq2v <- freqs[indsv, 2] # can have (many) zeroes, no NAs
    freq3v <- freqs[indsv, 3]
    pred1v <- preds[indsv, 1] # prediction
    pred2v <- preds[indsv, 2] # best non-prediction
    #
    # Compute PAC:
    #
    freqalt <- freq1v # frequency of prediction
    altintv <- pred1v # preliminary
    matched <- which(pred1v == yintv) # prediction == selfclass
    freqalt[matched] <- freq2v[matched] # frequency of alternative
    altintv[matched] <- pred2v[matched]
    PACv <- freqalt / (freq3v + freqalt)
    PAC[indsv] <- PACv # the other entries of PAC remain NA
    altint[indsv] <- altintv
    if (nlab == 2) altint <- 3 - yintnew # for consistency with svm etc.
    # When ynew is NA, also altint will be
  }
  #
  # Compute initial fig[i, g] from new cases to training classes:
  #
  fig <- matrix(rep(NA, n * nlab), ncol = nlab)
  for (i in seq_len(n)) { # loop over all cases in newdata
    for (g in seq_len(nlab)) { # loop over classes in training data
      ngbg <- which(yint[sortNgb[i, ]] == g)
      if (length(ngbg) > k) {
        ngbg <- ngbg[seq_len(k)]
      }
      fig[i, g] <- median(sortDis[i, ngbg])
    }
  }
  #
  # Compute farness:
  #
  farout <- compFarness(type = "knn", testdata = TRUE, yint = yintnew,
                        nlab = nlab, X = NULL, fig = fig, d = NULL,
                        figparams = vcr.knn.train.out$figparams)
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
              k = k,
              ktrues = ktrues,
              counts = counts))
}
