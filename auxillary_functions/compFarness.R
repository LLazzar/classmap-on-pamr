compFarness <- function(type = "affine", testdata = FALSE, yint, nlab,
                        X = NULL, fig = NULL, d = NULL, PCAfits = NULL,
                        keepPCA = FALSE, figparams = NULL) {
  #
  # Computes farness measures of different types.
  #
  # Arguments:
  #   type  : when "affine" this will standardize the Mahalanobis
  #           distances relative to each class.
  #           When "pca" it will compute the score distances and
  #           orthogonal distances to the PCA of each class.
  #   yint  : the class of each object, coded as 1, 2, ...
  #   X     : if not NULL, a data matrix with the x_i as rows.
  #   fig   : if not NULL, initial farness of each object i
  #           to each class g.
  #   d     : if not NULL, the dimension of the x_i.
  
  # Auxiliary functions:
  
  far2probability <- function(farness, trainInds = NULL) {
    #
    # converts farness values into probabilities
    # args:
    #   farness: a vector of farness values
    #   trainInds: if not NULL, these will be used in the estimation
    # returns:
    #   probs: the farness values converted to probabilities
    #   tfunc: a function which transforms the farness according to
    #   the estimated transformation (including all standardizations)
    #
    
    YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
      indlow <- which(y < 0)
      indhigh <- which(y >= 0)
      if (lambda != 0) {
        y[indhigh] <- ((1 + y[indhigh])^(lambda) - 1) / lambda
      }
      else {
        y[indhigh] <- log(1 + y[indhigh])
      }
      if (lambda != 2) {
        y[indlow] <- -((1 - y[indlow])^(2 - lambda) - 1) / (2 - lambda)
      }
      else {
        y[indlow] <- -log(1 - y[indlow])
      }
      if (stdToo) {
        if (length(y) > 1) {
          locScale <- cellWise::estLocScale(matrix(y, ncol = 1),
                                            type = "hubhub")
          zt <- (y - locScale$loc) / locScale$scale
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
    indnz  <- which(farness > 1e-10)
    farnz  <- farness[indnz]
    farloc <- median(farnz, na.rm = TRUE)
    farsca <- mad(farnz, na.rm = TRUE)
    if (farsca < 1e-10) {farsca <- sd(farnz, na.rm = TRUE)}
    sfar <- scale(farnz, center = farloc, scale = farsca)
    # Fit a Yeo-Johnson transformation:
    YJ.out  <- cellWise::transfo(X = sfar, robust = TRUE,
                                 prestandardize = FALSE,
                                 checkPars = list(silent = TRUE))
    xt      <- YJ.out$Xt
    tfarloc <- median(xt, na.rm = TRUE)
    tfarsca <- mad(xt, na.rm = TRUE)
    lambda  <- YJ.out$lambdahats
    #
    # Now transform origfarness (to be used in fig):
    origIndnz <- which(origfarness > 1e-10)
    origFarnz <- origfarness[origIndnz]
    zt        <- scale(YJ(scale(origFarnz, farloc, farsca),
                          lambda = lambda, stdToo = FALSE)$yt,
                       tfarloc, tfarsca)
    probs     <- rep(0,length(origfarness))
    probs[origIndnz] <- pnorm(zt)
    tfunc <- function(qs) {
      #
      # In Titanic data, 52% of casualties had duplicated features.
      # Therefore, we take out zero farness values in YJ computation:
      #
      qsIndnz <- which(qs > 1e-10)
      qsnz <- qs[qsIndnz]
      zn <- scale(YJ(scale(qsnz, farloc, farsca),
                     lambda = lambda, stdToo = FALSE)$yt,
                  tfarloc, tfarsca)
      pn <- rep(0,length(qs))
      pn[qsIndnz] <- pnorm(zn)
      return(pn)
    }
    return(list(probs = probs, far2prob = tfunc))
  }
  
  
  transformYJ <- function(X, lambda) {
    # Yeo-Johnson transformation with given parameter lambda.
    # The argument must be a vector or a number, without NA's.
    #
    if (!is.vector(X)) {
      stop("transformYJ: the data should be a vector or a number")
    }
    if (sum(is.na(X)) > 0) stop("transformYJ: the data contains NA's")
    Xt <- rep(NA, length(X))
    indlow  <- which(X < 0)
    indhigh <- which(X >= 0)
    if (lambda != 0) {
      Xt[indhigh] <- ((1 + X[indhigh])^(lambda) - 1) / lambda
    } else {
      Xt[indhigh] <- log(1 + X[indhigh])
    }
    if (lambda != 2) {
      Xt[indlow] <- -((1 - X[indlow])^(2 - lambda) - 1) / (2 - lambda)
    } else {
      Xt[indlow] <- -log(1 - X[indlow])
    }
    return(Xt)
  }
  
  
  
  mednz <- function(x, prec = 1e-12) {
    # Median of values > 0. Can be useful in denominators.
    #
    if (!(prec > 0)) stop("mednz: prec = ", prec, " should be > 0")
    xx <- x[which(!is.na(x))]
    xx <- xx[which(xx > prec)]
    if (length(xx) == 0) {mednzx <- prec} else {mednzx <- median(xx)}
    # By construction, always mednz >= prec .
    mednzx
  }
  
  
  # Here the main function compFarness starts:
  #
  if (!(type %in% c("affine", "knn", "pca"))) {
    stop(' type must be either "affine" or "knn" or "pca".')
  }
  n          <- length(yint)
  indsv      <- which(!is.na(yint))
  yintv      <- yint[indsv] # yint of cases that can be visualized
  ayint      <- sort(unique(yint[indsv])) # the values it takes
  farness    <- rep(NA, n) # will remain NA for cases outside indsv
  classSizes <- rep(0, nlab)
  for (g in seq_len(nlab)) {classSizes[g] <- sum(yintv == g)} #get class sizes
  medOD <- medSD <- SDs <- ODs <- NULL
  #
  if (type == "affine") {
    if (is.null(fig)) {stop("argument fig is missing")
    }
    #
    if (testdata == FALSE) { # this is for training data
      for (g in ayint) {
        clinds <- which(yint == g) # indices in all of 1, ..., n, label g index
        # Indices for which yint is NA do not get in clinds.
        farness[clinds] <- fig[clinds, g] #it's farness to true label
      }
      omedian <- mednz(farness, prec = 1e-08) # overall median
      #
      # Renormalize the classes:
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds          <- which(yint == g) # in 1, ..., n
        gfarness        <- farness[clinds] # all non-NA
        farscales[g]    <- mednz(gfarness, prec = 1e-08) / omedian
        farness[clinds] <- gfarness / farscales[g]
      }
      #
      # Fit a distribution to the homogenized farness measures:
      mfar    <- median(farness, na.rm = TRUE)
      sfar    <- max(mad(farness, na.rm = TRUE), 1e-08) # no division by 0
      farness <- as.vector(scale(farness, center = mfar, scale = sfar))
      tout    <- transfo(farness, type = "YJ", robust = TRUE,
                         prestandardize = FALSE, checkPars = list(silent = TRUE))
      lambda    <- tout$lambdahat[1]
      muhat     <- tout$muhat
      sigmahat  <- tout$sigmahat
      figparams <- list(farscales = farscales,
                        mfar = mfar,
                        sfar = sfar,
                        lambda = lambda,
                        muhat = muhat,
                        sigmahat = sigmahat)
    } else {# This is for new data
      # for new data, figparams should be given:
      if (is.null(figparams)) stop("Argument figparams is missing")
      farscales <- figparams$farscales
      if (is.null(farscales)) stop("figparams$farscales is missing")
      mfar     <- figparams$mfar
      sfar     <- figparams$sfar
      lambda   <- figparams$lambda
      muhat    <- figparams$muhat
      sigmahat <- figparams$sigmahat
    }
    fig <- scale(fig, center = FALSE, scale = farscales) # for training and newdata
    #
    # Now convert these distances to probabilities:
    for (g in seq_len(nlab)) { # cycle over columns of fig[, ]
      fig[, g] <- as.vector(scale(fig[, g], center = mfar, scale = sfar))
      # transform by YJ using lambda:
      fig[, g] <- transformYJ(fig[, g], lambda)
      fig[, g] <- as.vector(scale(fig[, g], center = muhat, scale = sigmahat)) #standardizing to get the normal distr
      fig[, g] <- pnorm(fig[, g])
    }
    for (g in ayint) { # for training and newdata
      clinds <- which(yint == g)
      farness[clinds] <- fig[clinds, g]
      # so it remains NA outside of indsv
    }
  } # end of type == "affine"
  #
  if (type == "knn") {
    if (is.null(fig)) {stop("argument fig is missing")
    }
    #
    if (testdata == FALSE) { # this is for training data
      # compute initial farness and overall median
      # farness = rep(NA, n) # already defined above
      for (g in ayint) {
        clinds <- which(yint == g) # indices in 1, ..., n
        farness[clinds] <- fig[clinds, g]
      }
      omedian <- mednz(farness, prec = 1e-08) # overall median
      #
      # Renormalize the classes:
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds       <- which(yint == g) # in all of 1, ..., n
        gfarness     <- farness[clinds] # all non-NA
        farscales[g] <- mednz(gfarness, prec = 1e-08) / omedian
        # farness[clinds] <- gfarness / farscales[g]
      }
      # wnq("farscales:")
      # pnq(farscales)
      fig <- scale(fig, center = FALSE, scale = farscales)
      farness <- rep(NA, n)
      far2probFuncs <- list()
      # list of estimated functions to transform raw farness values
      probs <- matrix(NA, nrow = nrow(fig), ncol = ncol(fig))
      for (g in seq_len(nlab)) {
        clinds             <- which(yint == g) # class indices
        ftemp              <- fig[, g]
        far2prob.out       <- far2probability(ftemp, trainInds = clinds)
        fig[, g]           <- far2prob.out$probs
        farness[clinds]    <- fig[clinds, g]
        far2probFuncs[[g]] <- far2prob.out$far2prob
      }
      figparams <- list(farscales = farscales,
                        far2probFuncs = far2probFuncs)
    } else {# this is for new data
      # for new data, figparams should be given:
      if (is.null(figparams)) stop("Argument figparams is missing")
      farscales <- figparams$farscales
      if (is.null(farscales)) stop("figparams$farscales is missing")
      fig <- scale(fig, center = FALSE, scale = farscales)
      #
      # transform the fig values using the far2prob functions:
      far2probFuncs <- figparams$far2probFuncs
      for (g in seq_len(nlab)) {
        far2probTemp <- far2probFuncs[[g]]
        fig[, g]      <- far2probTemp(fig[, g])
      }
      # now compute farness for all new data in indsv:
      for (g in ayint) { # g=3 # loop over classes in new data
        clinds <- which(yint == g)
        farness[clinds] <- fig[clinds, g]
        # so it remains NA outside of indsv
      }
    } # end of new data
  } # end of type == "knn"
  #
  if (type == "pca") {
    if (is.null(X)) {stop("argument X is missing")}
    if (is.vector(X)) {X <- matrix(X, ncol = 1)}
    X <- as.matrix(X)
    d <- ncol(X)
    if (nrow(X) != n) {stop(" X and yint have incompatible sizes")}
    fig <- SDs <- ODs <- matrix(rep(NA, n * nlab), ncol = nlab)
    #
    if (testdata == FALSE) { # this is for training data
      if (n < 2) stop(" X should have more than one row")
      PCAfits <- NULL
      if (keepPCA) PCAfits <- list()
      medSD <- medOD <- rep(NA, nlab)
      medscores <- madscores <- list()
      for (j in seq_len(nlab)) { # Loop over the classes # j=1
        clinds <- which(yint == j) # class indices
        # Cases with missing yint do not enter any clinds.
        Xj <- X[clinds, , drop = FALSE] # class j might contain only 1 case
        nXju <- sum(duplicated(Xj) == FALSE)
        cntr <- colMeans(Xj) # full dimension
        outpca <- prcomp(sweep(Xj, 2, cntr), scale = FALSE)
        if (nXju == 1) {outpca$rotation <- outpca$rotation * 0}
        sdev <- outpca$sdev # square roots of eigenvalues
        tokeep <- length(which(sdev > 1e-8))
        ncomps <- max(tokeep, 1)
        qmax <- NULL
        if (!is.null(qmax)) ncomps <- min(ncomps, qmax)
        # wnq(paste(" Keeping ", ncomps,
        #           " principal components in class ", j, sep=""))
        sdev <- sdev[seq_len(ncomps)]
        loadings <- outpca$rotation[, seq_len(ncomps), drop = FALSE]
        # columns are eigenvectors
        # store loadings, cntr, and sdev per class:
        if (keepPCA) PCAfits[[j]] <- list(lj = loadings, cj = cntr, sj = sdev)
        #
        dim(Xj); length(cntr); dim(loadings)
        pX <- sweep(X, 2, cntr) %*% loadings # projects all cases of X
        pX <- matrix(pX, nrow = n) # so it stays a matrix when ncomps=1
        dim(pX)
        #
        ## For orthogonal distances:
        if (tokeep == d) { # then no orthogonal distances exist
          ODs[, j] <- rep(0, n) # for cases in all classes
        } else {
          Xdiffs <- sweep(X, 2, cntr) - pX %*% t(loadings)
          # differences in full-dimensional space, for all cases
          ODs[, j] <- sqrt(rowSums(Xdiffs^2)) # Euclidean
          ODs[clinds, j] <- rep(0, length(clinds))
          # cases have no orthogonal distance to their own class
        }
        medOD[j] <- mednz(ODs[setdiff(indsv, clinds), j], prec = 1e-8)
        ODs[, j] <- ODs[, j] / medOD[j]
        # for all cases, but =zero in class j
        #
        ## For score distances:
        scores <- pX[clinds, , drop = FALSE]
        # projected, ONLY for cases in class j
        dim(scores)
        # distances from robustly scaled scores:
        medscores[[j]] <- apply(scores, 2, median)
        madscores[[j]] <- apply(scores, 2, mad)
        madscores[[j]] <- pmax(madscores[[j]], 1e-8)
        # clips at same precision as above, avoids dividing by 0
        scscores <- scale(pX, center = medscores[[j]], scale = madscores[[j]])
        # we will use this for all cases, not only in class j
        SDs[, j]  <- sqrt(rowSums(scscores^2))
        medSD[j] <- mednz(SDs[clinds, j], prec = 1e-8)
        # only over cases IN class j
        SDs[, j] <- SDs[, j] / medSD[j] # again for all cases
        # Special situation where sdev[1] ~ 0:
        if (tokeep == 0) SDs[, j] <- rep(1, length(SDs[, j]))
        ## Now combine SD with OD:
        fig[, j] <- sqrt(SDs[, j]^2 + ODs[, j]^2)
      } # ends loop over classes j
      #
      for (g in seq_len(nlab)) { # g=1
        clinds <- which(yint == g) # indices in 1, ..., n
        farness[clinds] <- fig[clinds, g] # only for available yint
      }
      omedian <- mednz(farness[indsv], prec = 1e-08) # overall median
      # omedian # 1 # This is because the farness is SDs^2 at this
      # stage, as ODs=0 in each class, and the SDs have already
      # been standardized to median 1 in each class. If no
      # farness is below 1e-08, mednz is the median, and if all
      # classes have odd cardinality we get omedian=1.
      #
      # Renormalize the classes to bring farness on the same scale:
      farscales <- rep(NA, nlab)
      for (g in ayint) {
        clinds       <- which(yint == g) # in 1, ..., n
        gfarness     <- farness[clinds] # all non-NA
        farscales[g] <- mednz(gfarness, prec = 1e-08) / omedian
        # These will also be 1 when no gfarness is below 1e-08,
        # otherwise not.
        farness[clinds] <- gfarness / farscales[g]
      }
      # Fit a distribution to the homogenized farness measures:
      mfar      <- median(farness[indsv], na.rm = TRUE)
      sfar      <- mad(farness[indsv], na.rm = TRUE)
      sfar      <- max(sfar, 1e-08) # avoids division by 0
      farness   <- as.vector(scale(farness, center = mfar, scale = sfar))
      tout      <- transfo(farness[indsv], type = "YJ", robust = TRUE,
                           prestandardize = FALSE,
                           checkPars = list(silent = TRUE))
      lambda    <- tout$lambdahat[1]
      muhat     <- tout$muhat
      sigmahat  <- tout$sigmahat
      figparams <- list(medOD = medOD,
                        medscores = medscores,
                        madscores = madscores,
                        medSD = medSD,
                        farscales = farscales,
                        mfar = mfar,
                        sfar = sfar,
                        lambda = lambda,
                        muhat = muhat,
                        sigmahat = sigmahat)
      # end of computing initial fig on training data
    } else {# test data, so we use the PCA fits of the training data
      #
      if (is.null(PCAfits)) {
        stop(paste0("\nFor test data, the farness computation ",
                    "requires the trained object","\nto contain ",
                    "the value $PCAfits. Rerun the training ",
                    "with keepPCA = TRUE"))
      }
      if (is.null(figparams)) stop("Argument figparams is missing")
      # if (is.null(farscales)) stop("Argument farscales is missing")
      medOD     <- figparams$medOD
      madscores <- figparams$madscores
      medscores <- figparams$medscores
      medSD     <- figparams$medSD
      farscales <- figparams$farscales
      #
      # Computations for new data apply the existing figparams:
      #
      for (j in seq_len(nlab)) { # Loop over the classes # j=1
        loadings  <- PCAfits[[j]]$lj
        cntr      <- PCAfits[[j]]$cj
        sdev      <- PCAfits[[j]]$sj
        tokeep    <- length(which(sdev > 1e-8))
        # tokeep can still be zero when ncomps is 1
        ncomps <- ncol(loadings)
        # wnq(paste0(" Using ", ncomps,
        #            " principal components in class ", j))
        clinds <- which(yint == j) # class indices
        Xj <- X[clinds, , drop = FALSE] # class j can contain only 1 case
        dim(Xj); length(cntr); dim(loadings)
        pX <- sweep(X, 2, cntr) %*% loadings # projects all cases of X
        scores <- pX[clinds, ] # projected, ONLY for cases in class j
        dim(pX); dim(scores)
        Xdiffs <- sweep(X, 2, cntr) - pX %*% t(loadings) # for ALL cases
        # differences in full-dimensional space
        ODs[, j] <- sqrt(rowSums(Xdiffs^2))
        ODs[, j] <- ODs[, j] / medOD[j] # for all cases, but =zero in class j
        scscores <- scale(pX, center = medscores[[j]],
                          scale = madscores[[j]])
        SDs[, j] <- sqrt(rowSums(scscores^2))
        SDs[, j] <- SDs[, j] / medSD[j] # again for all cases
        if (tokeep == 0) SDs[, j] <- rep(1, length(SDs[, j]))
        fig[, j] <- sqrt(SDs[, j]^2 + ODs[, j]^2)
      } # ends loop over classes j
    } # end of computing initial fig on test data
    #
    # The following is for both training data and test data:
    #
    fig <- scale(fig, center = FALSE, scale = farscales)
    #
    mfar      <- figparams$mfar
    sfar      <- figparams$sfar
    lambda    <- figparams$lambda
    muhat     <- figparams$muhat
    sigmahat  <- figparams$sigmahat
    # Now convert these distances to probabilities:
    for (g in seq_len(nlab)) { # cycle over columns of fig[, ]
      fig[, g] <- as.vector(scale(fig[, g], center = mfar, scale = sfar))
      # transform by YJ using lambda:
      fig[, g] <- transformYJ(fig[, g], lambda)
      fig[, g] <- as.vector(scale(fig[, g], center = muhat, scale = sigmahat))
      fig[, g] <- pnorm(fig[, g])
    }
    for (g in ayint) { # g=1
      clinds <- which(yint == g)
      farness[clinds] <- fig[clinds, g]
      # so it remains NA outside of indsv
    }
  } # end of type == "pca"
  #
  ofarness <- apply(fig, 1, min) # ofarness even exists for cases
  # with missing yint, for which farness cannot be defined.
  if (!keepPCA) PCAfits <- NULL
  return(list(fig = fig,
              farness = farness,
              ofarness = ofarness,
              figparams = figparams,
              PCAfits = PCAfits,
              medSD = medSD,
              medOD = medOD,
              SDs = SDs,
              ODs = ODs))
}