vcr.svm.train = function (X, y, svfit, ortho = FALSE){
  X <- as.matrix(X)
  if (nrow(X) == 1)
    X <- t(X)
  if (sum(is.na(as.vector(X))) > 0) {
    stop("The data matrix X has NA's.")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2)
    stop(" X should have more than one row")
  checked <- checkLabels(y, n, training = TRUE)
  lab2int <- checked$lab2int
  levels <- checked$levels
  nlab <- length(levels)
  yint <- lab2int(y)
  indsv <- checked$indsv
  nvis <- length(indsv)
  yintv <- yint[indsv]
  if (ortho == TRUE & nlab != 2)
    stop(" ortho=TRUE is only possible for 2 classes")
  if (is.null(svfit$compprob))
    stop(" svfit$compprob is missing")
  if (svfit$compprob == FALSE) {
    stop(paste0("\n The object svfit needs to be obtained by ",
                "\n running e1071::svm with probability=TRUE"))
  }
  if (is.null(svfit$kernel)) {
    stop(" svfit$kernel is missing")
  }
  else {
    k <- svfit$kernel
  }
  if (!(k == 0 | k == "linear" | k == 1 | k == "polynomial" |
        k == 2 | k == "radial" | k == 3 | k == "sigmoid")) {
    stop(paste(" Unknown kernel type: ", k, sep = ""))
  }
  if (is.null(svfit$scaled)) {
    stop(" svfit$scaled is missing")
  }
  else {
    scaled <- svfit$scaled
  }
  if (is.null(svfit$degree)) {
    stop(" svfit$degree is missing")
  }
  else {
    degree <- svfit$degree
  }
  if (is.null(svfit$gamma)) {
    stop(" svfit$gamma is missing")
  }
  else {
    gamma <- svfit$gamma
  }
  if (is.null(svfit$coef0)) {
    stop(" svfit$coef0 is missing")
  }
  else {
    coef0 <- svfit$coef0
  }
  if (is.null(svfit$cost)) {
    stop(" svfit$cost is missing")
  }
  else {
    cost <- svfit$cost
  }
  if (is.null(svfit$fitted) | is.null(svfit$decision.values)) {
    if (k == 0) {
      ktype <- "linear"
    }
    if (k == 1) {
      ktype <- "polynomial"
    }
    if (k == 2) {
      ktype <- "radial"
    }
    if (k == 3) {
      ktype <- "sigmoid"
    }
    svfit <- svm(X, y, scale = scaled, kernel = ktype, degree = degree,
                 gamma = gamma, coef0 = coef0, cost = cost, probability = TRUE)
  }
  preds <- predict(object = svfit, newdata = data.frame(X = X),
                   probability = TRUE)
  predint <- lab2int(preds)
  altint <- PAC <- rep(NA, n)
  altintv <- PACv <- rep(NA, nvis)
  if (nlab == 2) {
    altintv <- 3 - yintv
    probs <- attr(preds, "probabilities")[indsv, ]
    probs <- probs[, order(lab2int(colnames(probs)))]
    for (i in seq_len(nvis)) PACv[i] <- probs[i, altintv[i]]
  }
  if (nlab > 2) {
    probs <- attr(preds, "probabilities")[indsv, ]
    probs <- probs[, order(lab2int(colnames(probs)))]
    for (i in seq_len(nvis)) {
      ptrue <- probs[i, yintv[i]]
      others <- (seq_len(nlab))[-yintv[i]]
      palt <- max(probs[i, others])
      altintv[i] <- others[which.max(probs[i, others])]
      PACv[i] <- palt/(ptrue + palt)
    }
  }
  altint[indsv] <- altintv
  PAC[indsv] <- PACv
  X <- X[indsv, ]
  Xloadings <- transfmat <- nbeta <- NULL
  if (k == 0) {
    Xfar <- X
    PCout <- cellWise::truncPC(Xfar, ncomp = min(20, ncol(Xfar)), center = FALSE,
                               scale = FALSE)
    eigs <- PCout$eigenvalues[seq_len(PCout$rank)]
    qkeep <- sum(eigs > 1e-08)
    Xloadings <- PCout$loadings[, seq_len(qkeep), drop = FALSE]
    Xfar <- Xfar %*% Xloadings
  }
  else { # k >= 1
    Xfar <- makeKernel(X, svfit = svfit)
    transfmat <- makeFV(Xfar, precS = 1e-12)$transfmat
    if(ncol(transfmat) > 20) transfmat = transfmat[, seq_len(20)]
    Xfar <- Xfar %*% transfmat
  }
  if (ortho == TRUE) {
    fitfs <- svm(Xfar, y, scale = FALSE, kernel = "linear",
                 cost = cost)
    beta <- coef(fitfs)[-1]
    beta <- matrix(beta, ncol = 1)
    nbeta <- beta/sqrt(sum(beta * beta))
    Xfar <- Xfar - Xfar %*% nbeta %*% t(nbeta)
  }
  farout <- compFarness(type = "pca", testdata = FALSE,
                         yint = yintv, nlab = nlab, X = Xfar, fig = NULL, d = NULL,
                         figparams = NULL, PCAfits = NULL, keepPCA = TRUE)
  fig <- matrix(rep(0, n * nlab), ncol = nlab)
  fig[indsv, ] <- farout$fig
  farness <- ofarness <- rep(NA, n)
  farness[indsv] <- farout$farness
  ofarness[indsv] <- farout$ofarness
  figparams <- farout$figparams
  figparams$ortho <- ortho
  figparams$nbeta <- nbeta
  figparams$Xloadings <- Xloadings
  figparams$transfmat <- transfmat
  figparams$PCAfits <- farout$PCAfits
  return(list(yint = yint, y = levels[yint], levels = levels,
              predint = predint, pred = levels[predint], altint = altint,
              altlab = levels[altint], PAC = PAC, figparams = figparams,
              fig = fig, farness = farness, ofarness = ofarness, svfit = svfit,
              X = X))
}


vcr.svm.newdata = function (Xnew, ynew = NULL, vcr.svm.train.out){
  Xnew <- as.matrix(Xnew)
  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  if (is.null(vcr.svm.train.out$figparams$PCAfits)) {
    stop(paste0("\nFor test data, the farness computation ",
                "requires the trained object", "\nto contain ",
                "the value figparams$PCAfits."))
  }
  PCAfits <- vcr.svm.train.out$figparams$PCAfits
  if (sum(is.na(as.vector(Xnew))) > 0)
    stop("The data matrix Xnew contains NA's.")
  levels <- vcr.svm.train.out$levels
  nlab <- length(levels)
  if (is.null(ynew))
    ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = FALSE, levels)
  lab2int <- checked$lab2int
  indsv <- checked$indsv
  nvis <- length(indsv)
  yintnew <- rep(NA, n)
  yintnew[indsv] <- lab2int(ynew[indsv])
  yintv <- yintnew[indsv]
  svfit <- vcr.svm.train.out$svfit
  k <- svfit$kernel
  preds <- predict(object = svfit, newdata = data.frame(X = Xnew),
                   decision.values = TRUE, probability = TRUE)
  predint <- lab2int(preds)
  altint <- PAC <- rep(NA, n)
  if (nvis > 0) {
    altintv <- PACv <- rep(NA, nvis)
    if (nlab == 2) {
      altintv <- 3 - yintv
      probs <- attr(preds, "probabilities")[indsv,
      ]
      probs <- probs[, order(lab2int(colnames(probs)))]
      for (i in seq_len(nvis)) PACv[i] <- probs[i, altintv[i]]
    }
    if (nlab > 2) {
      probs <- attr(preds, "probabilities")[indsv,
      ]
      probs <- probs[, order(lab2int(colnames(probs)))]
      for (i in seq_len(nvis)) {
        ptrue <- probs[i, yintv[i]]
        others <- (seq_len(nlab))[-yintv[i]]
        palt <- max(probs[i, others])
        altintv[i] <- others[which.max(probs[i, others])]
        PACv[i] <- palt/(ptrue + palt)
      }
    }
    PAC[indsv] <- PACv
    altint[indsv] <- altintv
  }
  Zfar <- Xnew
  rm(Xnew)
  # dim(Zfar) # 75 2
  if (k == 0) {
    Xloadings <- vcr.svm.train.out$figparams$Xloadings
    if (!is.null(Xloadings)) {
      Zfar <- Zfar %*% Xloadings
    }
  }
  else {
    if (is.null(vcr.svm.train.out$figparams$transfmat)) {
      stop("vcr.svm.train.out$figparams$transfmat is missing")
    }
    Zfar <- makeFV(makeKernel(Zfar, vcr.svm.train.out$X,
                               svfit = vcr.svm.train.out$svfit), vcr.svm.train.out$figparams$transfmat)$Xf
    # dim(Zfar) # 75 196
  }
  if (vcr.svm.train.out$figparams$ortho == TRUE) {
    nbeta <- vcr.svm.train.out$figparams$nbeta
    Zfar <- Zfar - Zfar %*% nbeta %*% t(nbeta)
  }
  farout <- compFarness(type = "pca", testdata = TRUE,
                         yint = yintnew, nlab = nlab, X = Zfar, fig = NULL, d = NULL,
                         figparams = vcr.svm.train.out$figparams, PCAfits = PCAfits,
                         keepPCA = FALSE)
  return(list(yintnew = yintnew, ynew = levels[yintnew], levels = levels,
              predint = predint, pred = levels[predint], altint = altint,
              altlab = levels[altint], PAC = PAC, fig = farout$fig,
              farness = farout$farness, ofarness = farout$ofarness))
}

