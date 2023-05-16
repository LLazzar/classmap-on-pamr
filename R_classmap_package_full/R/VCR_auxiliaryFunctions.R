
#######################
# Auxiliary functions #
#######################


# Hidden from the user:
wnq <- function(string, qwrite = TRUE) { # auxiliary function
  # writes a line without quotes
  if (qwrite) write(noquote(string), file = "", ncolumns = 100)
}


# Hidden from the user:
pnq <- function(string, qwrite = TRUE) { # auxiliary function
  # prints without quotes
  if (qwrite) print(noquote(string))
}


# Hidden from the user:
checkLabels <- function(y, n, training = TRUE, levels = NULL) {
  if (training == TRUE) { # for training data
    if (!is.null(levels)) {
      stop("For training data, the levels argument is not used.")
    }
    whichy <- "y"
  } else {# for new data
    if (is.null(levels)) {
      stop("For new data, the argument levels is required")
    }
    whichy <- "ynew"
  }
  if (!(is.factor(y))) {
    stop(paste0("\n The response ", whichy, " is not a factor."))
  }
  if (is.null(attr(y, "levels"))) {
    stop(paste0("\n The response factor ", whichy, " has no levels.",
                "\n Please use as.factor() or factor() first."))
  }
  levelsattr <- attr(y, "levels")
  if (sum(is.na(levelsattr)) > 0) {
    stop(paste0("\n The levels attribute of ", whichy,
                " has at least one NA."))
  }
  dupls <- duplicated(levelsattr)
  if (sum(dupls) > 0) stop(paste0(
    "\n The label ", levelsattr[dupls == TRUE], " occurs more than ",
    " once in the levels attribute of ", whichy))
  if (length(y) != n) {
    stop(paste0("\n The response y should have length ", n,
                ", the number of cases."))
  }
  yv <- as.character(y) # this yv is not a factor any more.
  # When the factor y is malformed, the function stops here.
  # This happens when not all entries of y belong to levelsattr.
  indsv <- which(!is.na(yv)) # INDiceS of y that can be Visualized
  if (training == TRUE) { # for training data:
    if (length(indsv) == 0) {
      stop(paste0("The response factor y only contains NA's,"
                  , "\n so training is not possible."))
    }
    yv <- yv[indsv] # thus has at least 1 entry
    uniqy <- sort(unique(yv))
    for (g in seq_len(length(uniqy))) { # g = 1
      if (!(uniqy[g] %in% levelsattr)) {
        stop(paste0("\n The response label ", uniqy[g], " does not",
                    " appear in the levels attribute of y."))
      }
    }
    # From here on we know that all the yv are in levelsattr.
    if (length(uniqy) < 2) {
      stop(paste0("\n Only a single level occurs ( ", uniqy[1],
                  " ) so training is not possible."))
    }
    for (g in seq_len(length(levelsattr))) { # g = 1
      if (!(levelsattr[g] %in% uniqy)) {
        wnq(paste0("\nThe level ", levelsattr[g], " does not occur",
                   " among the training cases. A model trained",
                   "\non these data will not be able to",
                   " assign any objects to this class.\n"))
      }
    }
    # Create array "levels" with the labels that actually occur:
    levels <- levelsattr[which(levelsattr %in% uniqy)]; levels
    #
  } else {# for new data:
    if (length(indsv) > 0) {
      yv <- yv[indsv] # thus has at least 1 entry
      uniqy <- sort(unique(yv))
      for (g in seq_len(length(uniqy))) { # g = 1
        if (!(uniqy[g] %in% levelsattr)) {
          stop(paste0("\n The response label ", uniqy[g], " does not",
                      " appear in the levels attribute of ynew."))
        }
      }
      # From here on we know that all the yv are in levelsattr.
      badlevels <- uniqy[which(!(uniqy %in% levels))]; badlevels
      if (length(badlevels) > 0) { wnq(paste0(
        "\n The level ", badlevels, " occurs in ynew but",
        " was not in the training data.",
        "\n Such cases will be treated as if their ynew is NA",
        " so their response", "\n can be predicted, but",
        " these cases will not be in the class map.\n"))
        indsv <- indsv[which(!(yv %in% badlevels))]
        # the resulting indsv may again be empty
      }
    }
  }
  xlevelz <- levels # for use in the following function:
  lab2int <- function(labs){
    # labs is a vector of labels, that may contain NA's
    ints <- rep(NA, length(labs))
    indv <- which(!is.na(labs))
    labv <- labs[indv] # uses the final indsv above
    for (g in seq_len(length(xlevelz))) { # g = 1
      clinds <- indv[which(labv == xlevelz[g])] # in all of labs
      ints[clinds] <- g
    }
    ints
  }
  return(list(lab2int = lab2int, levels = levels, indsv = indsv))
}


# Hidden from the user:
meancov <- function(X) {
  # Wrapper around colMeans and classical covariance matrix
  return(list(m = colMeans(X), S = cov(X)))
}


# Hidden from the user:
DetMCD <- function(X) {
  # Wrapper around robustbase::covMCD
  detmcd.out <- robustbase::covMcd(X, nsamp = "deterministic", alpha = 0.5)
  return(list(m = detmcd.out$center, S = detmcd.out$cov))
}


# Hidden from the user:
numericalCond <- function(S) {
  # Computes the inverse condition number of a matrix, which
  # is the ratio: smallest eigenvalue/largest eigenvalue .
  # Unlike the CN, this avoids division by (near) zero.
  # I call this the "numerical condition" (NC).
  #
  S <- as.matrix(S)
  if (length(which(is.na(S))) > 0) stop(" S contains NA's.")
  d <- nrow(S)
  if (ncol(S) != d) stop(" S is not square.")
  maxS <- max(abs(as.vector(S)))
  if (!isSymmetric(S, tol = 1e-10 * maxS)) stop(" S is not symmetric.")
  eigvals <- eigen(S, only.values = TRUE)$values
  eigvals[d] / eigvals[1]
}


# Visible to the user (is called in example script):
makeKernel <- function(X1, X2 = NULL, svfit) {
  # Computes kernel value or kernel matrix, where the
  # kernel type is extracted from svfit.
  #
  # Arguments:
  # X1    : first data (coordinate) matrix or vector.
  # X2    : if not NULL, second data matrix or vector.
  # svfit : an object returned by e1071::svm(), or a list
  #         with the members $scaled, $kernel, $degree,
  #         $gamma, and $coef0 .
  #         Here $kernel may be in c(0, 1, 2, 3) or in
  #         c("linear", "polynomial", "radial", "sigmoid").
  #
  # A heuristic choice for gamma in radial and sigmoid is
  #    0.5/median(as.vector(as.matrix(dist(Xmat)))^2)
  #
  # Value:
  # kmat  : the kernel matrix or kernel value.
  #
  if (is.vector(X1)) {
    X1 <- t(as.matrix(X1))
  } else {
    X1 <- as.matrix(X1)
  }
  d <- ncol(X1)
  if (!is.null(X2)) {
    if (is.vector(X2)) {
      X2 <- t(as.matrix(X2))
    } else {
      X2 <- as.matrix(X2)
    }
    if (ncol(X2) != d) stop(" X1 and X2 have different dimensions")
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

  #

  if (length(scaled) != d) {
    stop(paste(" svfit$scaled should have length ", d, sep = ""))
  }
  # Now we scale the corresponding columns:
  sind <- which(scaled == TRUE)
  if (length(sind) > 0) {
    if (is.null(svfit$x.scale)) {
      stop(" svfit$x.scale is missing")
    } else {
      ctrs <- svfit$x.scale[[1]]
      scls <- svfit$x.scale[[2]]
      X1[, sind] <- scale(X1[, sind], center = ctrs, scale = scls)
      if (!is.null(X2)) {
        X2[, sind] <- scale(X2[, sind], center = ctrs, scale = scls)
      }
    }
  }
  #
  # kernel definitions
  #
  if (k == 0 | k == "linear") {
    kern <- function(xi, yj) {sum(xi * yj)}
  }
  if (k == 1 | k == "polynomial") {
    kern <- function(xi, yj) {(gamma * sum(xi * yj) + coef0)^degree}
  }
  if (k == 2 | k == "radial") {
    kern <- function(xi, yj) {exp(-gamma * sum((xi - yj) * (xi - yj)))}
  }
  if (k == 3 | k == "sigmoid") {
    kern <- function(xi, yj) {tanh(gamma * sum(xi * yj) + coef0)}
  }
  n1 <- nrow(X1)
  if (is.null(X2)) { # compute n1 by n1 matrix kmat
    kmat <- matrix(rep(0, n1 * n1), nrow = n1)
    dim(kmat)
    # Only fill triangular matrix and then put in mirror
    # image to obtain whole symmetric matrix:
    for (i in seq_len(n1)) {
      for (j in seq_len(i)) {
        kmat[i, j] <- kern(X1[i, ], X1[j, ])
      }
    }
    kmat <- kmat + t(kmat)
    diag(kmat) <- 0.5 * diag(kmat)
  }
  if (!is.null(X2)) { # compute n1 by n2 matrix kmat
    n2 <- nrow(X2)
    kmat <- matrix(rep(0, n1 * n2), nrow = n1)
    for (i in seq_len(n1)) {
      for (j in seq_len(n2)) {kmat[i, j] <- kern(X1[i, ], X2[j, ])
      }
    }
  }
  if (nrow(kmat) == 1 & ncol(kmat) == 1) kmat <- kmat[1, 1]
  return(kmat)
}


# Visible to the user:
makeFV <- function(kmat, transfmat = NULL, precS = 1e-12)
{
  kmat <- as.matrix(kmat)
  if (length(which(is.na(kmat))) > 0) {
    stop(" The kernel matrix kmat contains NA's.")
  }
  n <- nrow(kmat)
  maxK <- max(abs(as.vector(kmat)))
  if (is.null(transfmat)) { # so kmat contains training data
    if (ncol(kmat) != n) {
      stop(" The training kernel matrix kmat is not square.")
    }
    if (maxK < 0.001) {
      stop(" The training kernel matrix is too close to zero.")
    }
    if (max(abs(as.vector(kmat - t(kmat)))) > 1e-10 * maxK) {
      stop(" The training kernel matrix kmat is not symmetric.")
    }
    defaultPrecS <- 1e-16
    if (is.null(precS)) {
      precS <- defaultPrecS
    }
    if (precS < defaultPrecS) {
      precS <- defaultPrecS
    }
    eig <- eigen(kmat)
    eigvals <- pmax(eig$values, 0)
    eig <- eig$vectors
    V <- eig; rm(eig) # R's rename, to save memory
    keep <- length(which(eigvals > precS))
    if (keep < ncol(kmat)) {
      wnq(paste0(" The kernel matrix has ", ncol(kmat) -
                              keep, " eigenvalues below precS = ", precS,
                            "."))
      wnq(" Will take this into account.")
    }
    V <- V[,seq_len(keep)]
    colnames(V) <- NULL
    flipcolumn <- function(x) {
      if (x[which.max(abs(x))] < 0) {
        -x
      }
      else {
        x
      }
    }
    V <- apply(V, 2L, FUN = flipcolumn)
    scals <- sqrt(eigvals[seq_len(keep)])
    # Are the root sum of squares of the feature variables.
    Xf <- V %*% diag(scals) # This is safe because we do not
    # have to invert scals.
    V <- V %*% diag(1/scals)
    transfmat <- V; rm(V) # R's rename, to save memory
  }
  else { # so kmat contains new data
    dim <- nrow(transfmat)
    if (ncol(kmat) != dim) {
      stop(paste(" The matrix kmat has ", ncol(kmat), " columns",
                 " but the argument transfmat was made for ",
                 dim, " columns.", sep = ""))
    }
    Xf <- kmat %*% transfmat
    transfmat <- NULL # no need to output it again
  }
  return(list(Xf = Xf, transfmat = transfmat))
}

# Hidden from the user:
compFarness = function (type = "affine", testdata = FALSE, yint, nlab,
                         X = NULL, fig = NULL, d = NULL, PCAfits = NULL, keepPCA = FALSE,
                         figparams = NULL){
  far2probability <- function(farness, trainInds = NULL) {
    YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
      indlow <- which(y < 0)
      indhigh <- which(y >= 0)
      if (lambda != 0) {
        y[indhigh] <- ((1 + y[indhigh])^(lambda) - 1)/lambda
      }
      else {
        y[indhigh] <- log(1 + y[indhigh])
      }
      if (lambda != 2) {
        y[indlow] <- -((1 - y[indlow])^(2 - lambda) -
                         1)/(2 - lambda)
      }
      else {
        y[indlow] <- -log(1 - y[indlow])
      }
      if (stdToo) {
        if (length(y) > 1) {
          locScale <- cellWise::estLocScale(matrix(y,
                                                   ncol = 1), type = "hubhub")
          zt <- (y - locScale$loc)/locScale$scale
        }
        else {
          zt <- y
        }
      }
      else {
        zt <- NULL
      }
      return(list(yt = y, zt = zt))
    }
    origfarness <- farness
    if (!is.null(trainInds)) {
      farness <- farness[trainInds]
    }
    indnz <- which(farness > 1e-10)
    farnz <- farness[indnz]
    farloc <- median(farnz, na.rm = TRUE)
    farsca <- mad(farnz, na.rm = TRUE)
    if (farsca < 1e-10) {
      farsca <- sd(farnz, na.rm = TRUE)
    }
    sfar <- scale(farnz, center = farloc, scale = farsca)
    YJ.out <- cellWise::transfo(X = sfar, robust = TRUE, standardize = FALSE,
                                checkPars = list(silent = TRUE))
    xt <- YJ.out$Y
    tfarloc <- median(xt, na.rm = TRUE)
    tfarsca <- mad(xt, na.rm = TRUE)
    lambda <- YJ.out$lambdahats
    origIndnz <- which(origfarness > 1e-10)
    origFarnz <- origfarness[origIndnz]
    zt <- scale(YJ(scale(origFarnz, farloc, farsca), lambda = lambda,
                   stdToo = FALSE)$yt, tfarloc, tfarsca)
    probs <- rep(0, length(origfarness))
    probs[origIndnz] <- pnorm(zt)
    tfunc <- function(qs) {
      qsIndnz <- which(qs > 1e-10)
      qsnz <- qs[qsIndnz]
      zn <- scale(YJ(scale(qsnz, farloc, farsca), lambda = lambda,
                     stdToo = FALSE)$yt, tfarloc, tfarsca)
      pn <- rep(0, length(qs))
      pn[qsIndnz] <- pnorm(zn)
      return(pn)
    }
    return(list(probs = probs, far2prob = tfunc))
  }

  transformYJ <- function(X, lambda) {
    if (!is.vector(X)) {
      stop("transformYJ: the data should be a vector or a number")
    }
    if (sum(is.na(X)) > 0)
      stop("transformYJ: the data contains NA's")
    Xt <- rep(NA, length(X))
    indlow <- which(X < 0)
    indhigh <- which(X >= 0)
    if (lambda != 0) {
      Xt[indhigh] <- ((1 + X[indhigh])^(lambda) - 1)/lambda
    }
    else {
      Xt[indhigh] <- log(1 + X[indhigh])
    }
    if (lambda != 2) {
      Xt[indlow] <- -((1 - X[indlow])^(2 - lambda) - 1)/(2 -
                                                           lambda)
    }
    else {
      Xt[indlow] <- -log(1 - X[indlow])
    }
    return(Xt)
  }
  mednz <- function(x, prec = 1e-12) {
    if (!(prec > 0))
      stop("mednz: prec = ", prec, " should be > 0")
    xx <- x[which(!is.na(x))]
    xx <- xx[which(xx > prec)]
    if (length(xx) == 0) {
      mednzx <- prec
    }
    else {
      mednzx <- median(xx)
    }
    mednzx
  }

  if (!(type %in% c("affine", "knn", "pca"))) {
    stop(" type must be either \"affine\" or \"knn\" or \"pca\".")
  }
  n <- length(yint)
  indsv <- which(!is.na(yint))
  yintv <- yint[indsv]
  ayint <- sort(unique(yint[indsv]))
  farness <- rep(NA, n)
  classSizes <- rep(0, nlab)
  for (g in seq_len(nlab)) {
    classSizes[g] <- sum(yintv == g)
  }
  medOD <- medSD <- SDs <- ODs <- NULL
  if (type == "affine") {
    if (is.null(fig)) {
      stop("argument fig is missing")
    }
    if (testdata == FALSE) {
      for (g in ayint) {
        clinds <- which(yint == g)
        farness[clinds] <- fig[clinds, g]
      }
      omedian <- mednz(farness, prec = 1e-08)
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds <- which(yint == g)
        gfarness <- farness[clinds]
        farscales[g] <- mednz(gfarness, prec = 1e-08)/omedian
        farness[clinds] <- gfarness/farscales[g]
      }
      mfar <- median(farness, na.rm = TRUE)
      sfar <- max(mad(farness, na.rm = TRUE), 1e-08)
      farness <- as.vector(scale(farness, center = mfar,
                                 scale = sfar))
      tout <- cellWise::transfo(farness, type = "YJ", robust = TRUE,
                                standardize = FALSE,
                                checkPars = list(silent = TRUE))
      lambda <- tout$lambdahats[1]
      muhat <- tout$muhat
      sigmahat <- tout$sigmahat
      figparams <- list(farscales = farscales, mfar = mfar,
                        sfar = sfar, lambda = lambda, muhat = muhat,
                        sigmahat = sigmahat)
    }
    else {
      if (is.null(figparams))
        stop("Argument figparams is missing")
      farscales <- figparams$farscales
      if (is.null(farscales))
        stop("figparams$farscales is missing")
      mfar <- figparams$mfar
      sfar <- figparams$sfar
      lambda <- figparams$lambda
      muhat <- figparams$muhat
      sigmahat <- figparams$sigmahat
    }
    fig <- scale(fig, center = FALSE, scale = farscales)
    for (g in seq_len(nlab)) {
      fig[, g] <- as.vector(scale(fig[, g], center = mfar,
                                  scale = sfar))
      fig[, g] <- transformYJ(fig[, g], lambda)
      fig[, g] <- as.vector(scale(fig[, g], center = muhat,
                                  scale = sigmahat))
      fig[, g] <- pnorm(fig[, g])
    }
    for (g in ayint) {
      clinds <- which(yint == g)
      farness[clinds] <- fig[clinds, g]
    }
  }# ends "affine"
  if (type == "knn") {
    if (is.null(fig)) {
      stop("argument fig is missing")
    }
    if (testdata == FALSE) {
      for (g in ayint) {
        clinds <- which(yint == g)
        farness[clinds] <- fig[clinds, g]
      }
      omedian <- mednz(farness, prec = 1e-08)
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds <- which(yint == g)
        gfarness <- farness[clinds]
        farscales[g] <- mednz(gfarness, prec = 1e-08)/omedian
      }
      fig <- scale(fig, center = FALSE, scale = farscales)
      farness <- rep(NA, n)
      far2probFuncs <- list()
      probs <- matrix(NA, nrow = nrow(fig), ncol = ncol(fig))
      for (g in seq_len(nlab)) {
        clinds <- which(yint == g)
        ftemp <- fig[, g]
        far2prob.out <- far2probability(ftemp, trainInds = clinds)
        fig[, g] <- far2prob.out$probs
        farness[clinds] <- fig[clinds, g]
        far2probFuncs[[g]] <- far2prob.out$far2prob
      }
      figparams <- list(farscales = farscales, far2probFuncs = far2probFuncs)
    }
    else {
      if (is.null(figparams))
        stop("Argument figparams is missing")
      farscales <- figparams$farscales
      if (is.null(farscales))
        stop("figparams$farscales is missing")
      fig <- scale(fig, center = FALSE, scale = farscales)
      far2probFuncs <- figparams$far2probFuncs
      for (g in seq_len(nlab)) {
        far2probTemp <- far2probFuncs[[g]]
        fig[, g] <- far2probTemp(fig[, g])
      }
      for (g in ayint) {
        clinds <- which(yint == g)
        farness[clinds] <- fig[clinds, g]
      }
    }
  } # ends "knn"
  if (type == "pca") {
    if (is.null(X)) {
      stop("argument X is missing")
    }
    if (is.vector(X)) {
      X <- matrix(X, ncol = 1)
    }
    X <- as.matrix(X)
    dim(X) # 75 196
    d <- ncol(X)
    if (nrow(X) != n) {
      stop(" X and yint have incompatible sizes")
    }
    fig <- SDs <- ODs <- matrix(rep(NA, n * nlab), ncol = nlab)
    if (testdata == FALSE) {
      if (n < 2)
        stop(" X should have more than one row")
      PCAfits <- NULL
      if (keepPCA)
        PCAfits <- list()
      medSD <- medOD <- rep(NA, nlab)
      medscores <- madscores <- list()
      for (j in seq_len(nlab)) {
        clinds <- which(yint == j)
        Xj <- X[clinds, , drop = FALSE]
        nXju <- sum(duplicated(Xj) == FALSE)
        cntr <- colMeans(Xj)
        outpca <- prcomp(sweep(Xj, 2, cntr), scale = FALSE)
        if (nXju == 1) {
          outpca$rotation <- outpca$rotation * 0
        }
        sdev <- outpca$sdev
        tokeep <- length(which(sdev > 1e-08))
        ncomps <- max(tokeep, 1)
        sdev <- sdev[seq_len(ncomps)]
        loadings <- outpca$rotation[, seq_len(ncomps),
                                    drop = FALSE]
        if (keepPCA)
          PCAfits[[j]] <- list(lj = loadings, cj = cntr,
                               sj = sdev)
        dim(Xj)
        length(cntr)
        dim(loadings)
        pX <- sweep(X, 2, cntr) %*% loadings
        pX <- matrix(pX, nrow = n)
        dim(pX)
        if (tokeep == d) {
          ODs[, j] <- rep(0, n)
        }
        else {
          Xdiffs <- sweep(X, 2, cntr) - pX %*% t(loadings)
          ODs[, j] <- sqrt(rowSums(Xdiffs^2))
          ODs[clinds, j] <- rep(0, length(clinds))
        }
        medOD[j] <- mednz(ODs[setdiff(indsv, clinds),
                              j], prec = 1e-08)
        ODs[, j] <- ODs[, j]/medOD[j]
        scores <- pX[clinds, , drop = FALSE]
        dim(scores)
        medscores[[j]] <- apply(scores, 2, median)
        madscores[[j]] <- apply(scores, 2, mad)
        madscores[[j]] <- pmax(madscores[[j]], 1e-08)
        scscores <- scale(pX, center = medscores[[j]],
                          scale = madscores[[j]])
        SDs[, j] <- sqrt(rowSums(scscores^2))
        medSD[j] <- mednz(SDs[clinds, j], prec = 1e-08)
        SDs[, j] <- SDs[, j]/medSD[j]
        if (tokeep == 0)
          SDs[, j] <- rep(1, length(SDs[, j]))
        fig[, j] <- sqrt(SDs[, j]^2 + ODs[, j]^2)
      }
      for (g in seq_len(nlab)) {
        clinds <- which(yint == g)
        farness[clinds] <- fig[clinds, g]
      }
      omedian <- mednz(farness[indsv], prec = 1e-08)
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds <- which(yint == g)
        gfarness <- farness[clinds]
        farscales[g] <- mednz(gfarness, prec = 1e-08)/omedian
        farness[clinds] <- gfarness/farscales[g]
      }
      mfar <- median(farness[indsv], na.rm = TRUE)
      sfar <- mad(farness[indsv], na.rm = TRUE)
      sfar <- max(sfar, 1e-08)
      farness <- as.vector(scale(farness, center = mfar,
                                 scale = sfar))
      tout <- cellWise::transfo(farness[indsv], type = "YJ", robust = TRUE,
                                standardize = FALSE, checkPars = list(silent = TRUE))
      lambda <- tout$lambdahats[1]
      muhat <- tout$muhat
      sigmahat <- tout$sigmahat
      figparams <- list(medOD = medOD, medscores = medscores,
                        madscores = madscores, medSD = medSD, farscales = farscales,
                        mfar = mfar, sfar = sfar, lambda = lambda, muhat = muhat,
                        sigmahat = sigmahat)
    } # ends if(training)
    else { # if(newdata)
      if (is.null(PCAfits)) {
        stop(paste0("\nFor test data, the farness computation ",
                    "requires the trained object", "\nto contain ",
                    "the value $PCAfits. Rerun the training ",
                    "with keepPCA = TRUE"))
      }
      if (is.null(figparams))
        stop("Argument figparams is missing")
      medOD <- figparams$medOD
      madscores <- figparams$madscores
      medscores <- figparams$medscores
      medSD <- figparams$medSD
      farscales <- figparams$farscales
      for (j in seq_len(nlab)) { # j=1
        loadings <- PCAfits[[j]]$lj
        dim(loadings) # 20 20
        cntr <- PCAfits[[j]]$cj
        sdev <- PCAfits[[j]]$sj
        tokeep <- length(which(sdev > 1e-08))
        ncomps <- ncol(loadings)
        clinds <- which(yint == j)
        Xj = X
        if(nrow(loadings) < ncol(Xj)) {
          Xj = Xj[, seq_len(nrow(loadings))]
        }
        pX <- sweep(Xj, 2, cntr) %*% loadings
        Xdiffs <- sweep(Xj, 2, cntr) - pX %*% t(loadings)
        ODs[, j] <- sqrt(rowSums(Xdiffs^2))
        ODs[, j] <- ODs[, j]/medOD[j]
        scscores <- scale(pX, center = medscores[[j]],
                          scale = madscores[[j]])
        SDs[, j] <- sqrt(rowSums(scscores^2))
        SDs[, j] <- SDs[, j]/medSD[j]
        if (tokeep == 0)
          SDs[, j] <- rep(1, length(SDs[, j]))
        fig[, j] <- sqrt(SDs[, j]^2 + ODs[, j]^2)
      }
    }
    fig <- scale(fig, center = FALSE, scale = farscales)
    mfar <- figparams$mfar
    sfar <- figparams$sfar
    lambda <- figparams$lambda
    muhat <- figparams$muhat
    sigmahat <- figparams$sigmahat
    for (g in seq_len(nlab)) {
      fig[, g] <- as.vector(scale(fig[, g], center = mfar,
                                  scale = sfar))
      fig[, g] <- transformYJ(fig[, g], lambda)
      fig[, g] <- as.vector(scale(fig[, g], center = muhat,
                                  scale = sigmahat))
      fig[, g] <- pnorm(fig[, g])
    }
    for (g in ayint) {
      clinds <- which(yint == g)
      farness[clinds] <- fig[clinds, g]
    }
  }
  ofarness <- apply(fig, 1, min)
  if (!keepPCA)
    PCAfits <- NULL
  return(list(fig = fig, farness = farness, ofarness = ofarness,
              figparams = figparams, PCAfits = PCAfits, medSD = medSD,
              medOD = medOD, SDs = SDs, ODs = ODs))
}

###############################################################################


checkXnew <- function(Xtrain, Xnew, varsInModel = NULL,
                      allowNA = FALSE) {
  # This function checks whether Xnew follows the same format
  # as Xtrain. It checks variable names, variable types,
  # data.frame vs. matrix, and factor levels.
  #
  # varsInModel is a vector of indices indicating which
  # variables are used in the model, so these should be
  # present in Xnew in the same format as in the training data.
  # If NULL, it is assumed that all variables should be checked.
  #
  #
  # Check object type: data frame vs. matrix
  if (is.data.frame(Xtrain)) {
    if (!is.data.frame(Xnew)) stop(paste0(
      "\nThe test data should be a data frame, ",
      " like the training data."))
  }
  if (is.matrix(Xtrain)) {
    if (!is.matrix(Xnew)) stop(paste0(
      "\nThe test data should be a matrix, ",
      " like the training data."))
  }
  #
  # Check number of variables:
  if (is.null(varsInModel) | identical(varsInModel,
                                       seq_len(ncol(Xtrain)))) {
    varsInModel <- seq_len(ncol(Xtrain))
    if (ncol(Xtrain) != ncol(Xnew)) stop(paste0(
      "\nThe test data should have ", ncol(Xtrain), " variables, ",
      " like the training data."))
  }
  #
  # Only compare the variables actually used in the model:
  Xtrain <- Xtrain[, varsInModel]
  #
  # Check variable names:
  #
  namesmatch <- match(colnames(Xtrain), colnames(Xnew))
  # namesmatch
  if (anyNA(namesmatch)) {
    message("\nThe following variable name(s) are missing in Xnew:")
    print(colnames(Xtrain)[is.na(namesmatch)])
    stop(paste0("\nThe variable names of the test data should",
                "\nmatch those of the training data."))
  } else {
    Xnew <- Xnew[, namesmatch] # puts the columns of Xnew in the
    # same order, and only keeps those columns we actually need.
  }
  # From here on the variable names match, and are in the same order.
  #
  if (!allowNA) {
    if (anyNA(Xnew)) stop(paste0(
      "\nXnew contains at least one NA in a",
      " variable used in the model."))
  }
  #
  # Check variable types:
  #
  if (is.data.frame(Xtrain)) {
    # if (!identical(sapply(Xtrain, data.class),
    #                sapply(Xnew, data.class) )) { }
    differ <- which(sapply(Xtrain, data.class) !=
                      sapply(Xnew, data.class))
    if (length(differ) > 0) {
      message("\nThe type of the following variable(s) in Xnew:")
      print(colnames(Xnew)[differ])
      message("differs from their type in the training data.")
      stop(paste0("\nThe variable types of the test data should",
                  "\nmatch those of the training data."))
    }
  }
  if (is.matrix(Xtrain)) {
    if (mode(Xnew) != mode(Xtrain)) stop(paste0(
      "\nThe mode of the test data matrix is ", mode(Xnew),
      ", \nbut that of the training data was ", mode(Xtrain), "."))
  }
  #
  # Check levels of factors and of character variables:
  #
  if (is.data.frame(Xtrain)) {
    varTypes <- sapply(Xtrain, data.class)
    #
    # First the factors:
    #
    facVars  <- which(varTypes == "factor")
    for (j in seq_len(length(facVars))) {
      charvec <- as.character(Xnew[, facVars[j]])
      # if the test factor is malformed, this throws an error.
      # It is not clear how to make the code report which
      # factor is responsible before the error is thrown.

      # check whether there are any new levels in the test factor
      levelsTrain <- attr(Xtrain[, facVars[j]], "levels")
      levelsTest  <- attr(Xnew[, facVars[j]], "levels")
      # check for new levels in the test variable:
      if (any(!(levelsTest %in% levelsTrain))) stop(
        paste0("\nThe factor \"", colnames(Xnew)[facVars[j]],
               "\" of the test data has at least one level",
               "\nthat does not appear in the training data."))
      # check for duplicate levels in the test variable:
      if (any(duplicated(levelsTest))) stop(
        paste0("\nThe factor \"", colnames(Xnew)[facVars[j]],
               "\" of the test data has duplicate levels."))
    }
    #
    # Now for the character variables. Note that this function
    # should be called with the variables actually used in the
    # modeling, and not e.g. "name" variables.
    #
    charVars <- which(varTypes == "character")
    for (j in seq_len(length(charVars))) {
      # Check whether there are any new values in the factor:
      valsTrain <- unique(Xtrain[, charVars[j]])
      valsTest  <- unique(Xnew[, charVars[j]])
      if (allowNA) {valsTrain <- c(valsTrain, NA)}
      #
      # check for new levels in the test variable
      if (any(!(valsTest %in% valsTrain))) stop(paste0(
        "\nThe character variable \"", colnames(Xnew)[charVars[j]],
        "\" of the test data has", "\nat least one value",
        " that does not appear in the training data."))
    }
  }
  #
  if (is.matrix(Xtrain)) {
    if ((mode(Xtrain) == "character")) {
      # a matrix of character variables
      for (j in seq_len(ncol(Xtrain))) {
        # Check whether there are any new values in the factor:
        valsTrain <- unique(Xtrain[, j])
        valsTest  <- unique(Xnew[, j])
        if (allowNA) {
          valsTrain <- c(valsTrain, NA)
        }
        # check for new levels in test variable
        if (any(!(valsTest %in% valsTrain))) stop(paste0(
          "\nThe character variable \"", colnames(Xnew)[j],
          "\" of the test data has", "\nat least one value",
          " that does not appear in the training data."))
      }
    }
  }
}
