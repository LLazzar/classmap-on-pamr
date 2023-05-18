######################### (code that will not be needed in case of classmap merge)
#to import functions needed for VCR_pamr (especially for checkLabels, computeFarness..)
library(cellWise) #for transfo function used in Comp fareness in VCR_auxillaryFunctions.R
source("R_classmap_package_full/R/VCR_auxiliaryFunctions.R") #importing auxillary functions needed
                                                             # this script is available in classmap package
                                                             # so in case of integration of VCR_pamr this import would be useless


library(pamr) #to get pamr.predict for newdata function, in package can just import that with import pamr.predict from pamr
#########################
source("feature_code/R/VCR_auxiliaryFunctions_alt.R") # an alternative version of the file used for some experiments
                                                      # for example i tried to removed the division by omedian in compFareness

vcr.pamr.train <- function(data, pamrfit, pamrfitcv=NULL, threshold) {
  #
  # Using the outputs of just pamr.train (or of pamr.cv) for classification
  # applied to the training data, this function prepares for
  # graphical displays.
  #
  ##########
  # TO DO LIST:
  # MISSING vcr.pamr.newdata
  #########
  #
  # putted pamrfit object just like forest.vcr takes forest fit
  #
  #
  # Arguments:
  #   data      : the same input data to pamr (same format) that is a
  #               list with components: x- an expression genes in the rows,
  #               samples in the columns), and y- a vector of the class labels
  #               for each sample. Optional components- genenames, a vector
  #               of gene names, and geneid- a vector of gene identifiers.
  #   pamrfit   : the result of a call to pamr.train or pamr.cv (PAMR.CV TO BE TESTED)
  #   threshold : the desired threshold value
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
  #
  # Subsetting to the same subset (of variables and observation) on which pamr fit works on.

  #if (!is.null(pamrfit$gene.subset)) {

  data$x=data$x[pamrfit$gene.subset,pamrfit$sample.subset] #can subset for both genes and samples #PROBABLY SUBSETTING GENE NOT USEFUL IN OURCASE
  #} else  { #because pamr.cv does not have directly gene.subset in the output
    #data$x=data$x[,pamrfit$sample.subset] ###PROBLEM TO GET GENE SUBSET but problably not useful pamrcv already inherits gene.subset
                                            # problem in matrix moltiplication in DD function if there subset of gene (t(X)%*%centroids)
                                            #centroids vector is shrinked into gene sub dimension, t(x) is not
  #}
  #
  #
  X <- as.matrix(t(data$x)) # in case it is a data frame
                            # also transpose back since pamr takes rows as variables and columns as observation
                            # IS THIS NECESSARY??
  if (nrow(X) == 1) X <- t(X)
  if (sum(is.na(as.vector(X))) > 0) {
    stop("The coordinate matrix X has NA's.") #it's ok to leave that because pamr don't fit with NAs
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2) stop("The training data should have more than one case.")

  y=as.factor(data$y) #FACTORIZE IT OR NOT?

  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = TRUE) #PROBABLY SHOULD RE DIG DEEP TO UNDERSTAND THIS FUNCTION
  # If it did not stop: yint has length n, and its values are in
  # 1, ..., nlab without gaps. It is NA where y is: outside indsv.
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv
  levels  <- checked$levels
  nlab    <- length(levels)
  yint    <- lab2int(y) #given label (true) as integer
  yintv   <- yint[indsv]
  #
  #
  #
  # Getting threshold index from inputted threshold value (idea of this code/logic from pamr.confusion)
  #
  #

  if (is.null(pamrfitcv)) {
    ii <- (1:length(pamrfit$threshold))[pamrfit$threshold >= threshold] ##ADD STOP IF THRESHOLD VALUE IS OUTSIDE
    ii <- ii[1] #taking the first in the list
  } else {
    ii <- (1:length(pamrfitcv$threshold))[pamrfit$threshold >= threshold] ##ADD STOP IF THRESHOLD VALUE IS OUTSIDE
    ii <- ii[1] #taking the first in the list
  }

  #
  # Check matrix of posterior probabilities:
  #

  if (is.null(pamrfitcv)) {
  probs <- as.matrix(pamrfit$prob[,,ii]) #prob object in pamr output is consistent for our purpose

  } else {
    probs <- as.matrix(pamrfitcv$prob[,,ii]) #if a crossvalidated object is given, we want to evaluate cv results
  }
  if (length(dim(probs)) != 2) stop("probs should be a matrix.")
  if (nrow(probs) != n) stop(paste0(
    "The matrix probs should have ", n, " rows"))
  if (ncol(probs) != nlab) stop(paste0(
    "The matrix probs should have ", nlab, " columns"))
  if (any(is.na(probs))) stop("probs should not have any NA's.")
  #
  # Compute prediction for all objects in the training data:
  #
  # MAYBE SHOULD ADD LINE 89 VCR_FOREST (CHECK labels switching)
  #
  predint <- apply(probs[, , drop = FALSE], 1, which.max) #should be ok but check on pamr if this value corresponds to the yhat
  #                                                       # CAN BE PROBABLY REDUCE LIKE VCR_FOREST
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
  # GETTING DISCRIMINANT SCORES (they are our distance measures on which we build farness)
  #
  #

  #### Auxillary functions needed: ####

  soft.shrink <-function(delta, threshold) {
    dif <- abs(delta) - threshold
    delta <- sign(delta) * dif * (dif > 0)
    nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
    attr(delta, "nonzero") <- nonzero
    delta
  }

  diag.disc.original <-function(x, centroids, prior, weight) {
    ### Computes the class discriminant functions assuming scaled x and centroids
    if(! missing(weight)) {
      posid <- (weight > 0)
      if(any(posid)) {
        weight <- sqrt(weight[posid])
        centroids <- centroids[posid,  , drop = FALSE] * weight
        x <- x[posid,  , drop = FALSE] * weight
      }
      else {
        mat <- outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    dd <- t(x) %*% centroids
    dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
    names(dd0) <- NULL
    scale(dd, dd0, FALSE)
  }

  mdS2 <-function(x, centroids, prior, weight) {
    ### Computes the class discriminant functions assuming scaled x and centroids
    if(! missing(weight)) {
      posid <- (weight > 0)
      if(any(posid)) {
        weight <- sqrt(weight[posid])
        centroids <- centroids[posid,  , drop = FALSE] * weight
        x <- x[posid,  , drop = FALSE] * weight
      }
      else {
        mat <- outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    p=ncol(t(x))
    n=nrow(t(x))
    k=ncol(centroids)
    dd=matrix(NA, nrow=n, ncol=k)
    for (k in 1:ncol(centroids)){
      dd[,k]=mahalanobis(t(x),centroids[,k],cov=diag(p))
    }
    dd
  }
  ###################

  #actually reconstrunctanting dd
  #getting quantities needed from pamrfit

  if (is.null(pamrfitcv)) { #it means we have a classic pamr.train object


  centroids=pamrfit$centroids #centroids per variable per class
  centroid.overall=pamrfit$centroid.overall
  sd=pamrfit$sd
  threshold=pamrfit$threshold[ii]
  se.scale=pamrfit$se.scale
  threshold.scale=pamrfit$threshold.scale
  prior=pamrfit$prior
  nonzero=pamrfit$nonzero[ii]
  K=length(prior)

  #getting deltas (dik)
  delta <- (centroids - centroid.overall)/sd
  delta <- scale(delta, FALSE, threshold.scale * se.scale) #gives division by mk
  #getting the shrunken ones (d'ik)
  delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
  #getting d'ik*mk
  delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))
  nonzero_check <- attr(delta.shrunk, "nonzero")
  if(!nonzero_check==nonzero){
    stop(nonzero_check)
  }
  posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
  xtest<-data$x #NBNB we evaluate directly to train set here #in VCR.pamr.newdata or pamrcv probably needs to be modified
  dd <- mdS2((xtest - centroid.overall)/sd, delta.shrunk, prior, weight = posid)

  } else { #it means we have a pamr.cv, so different process to reconstruct dd

    dd=matrix(NA, nrow=n , ncol=nlab) #defining the matrix contain dd (n x k)
    for (nf in 1:length(pamrfitcv$folds)) {

      pamrfits=pamrfitcv$cv.objects[[nf]] #retrieve training object for the given fold in the loop
      centroids=pamrfits$centroids #centroids per variable per class
      centroid.overall=pamrfits$centroid.overall
      sd=pamrfits$sd
      threshold=pamrfits$threshold[ii]
      se.scale=pamrfits$se.scale
      threshold.scale=pamrfits$threshold.scale
      prior=pamrfits$prior
      nonzero=pamrfits$nonzero[ii]
      K=length(prior)

      #getting deltas (dik)
      delta <- (centroids - centroid.overall)/sd
      delta <- scale(delta, FALSE, threshold.scale * se.scale) #gives division by mk
      #getting the shrunken ones (d'ik)
      delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
      #getting d'ik*mk
      delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))
      nonzero_check <- attr(delta.shrunk, "nonzero")
      if(!nonzero_check==nonzero){
        stop(nonzero_check)
      }
      posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      xtest<-data$x[,pamrfitcv$folds[[nf]]] #I want dd only in the obvs considered as test in that fold
                                         # this is done following the logic inside nsccv ( called in pamr.cv) where for each fold nsc is called with x (obvs in train) and xtest the obvs considered as test
                                         # then in nsc xtest in called in dd (as here below)
      dd[pamrfitcv$folds[[nf]],] <- mdS2((xtest - centroid.overall)/sd, delta.shrunk, prior, weight = posid)


    }
  }

  if (any(is.na(dd))) { ##CHECK PUT FOR DEVELOPING REASONS PROBABLY COULD REMOVE
    stop(dd)
  }

  initfig<-(-dd) #minus because wrt to paper is with opposite sign here

  #if (min(initfig)<0){ #shift to all positive values to try box cox transformation
  #  initfig=initfig-(min(initfig))+3000
  #}

  #initfig=initfig+10 #TRYING TRANSFORMING DD TO SEE RESULTS
  #initfig=log(initfig)

  rd=pamrfit$nonzero[ii] #getting reducted dimension
  farout <- compFarness(type = "affine", testdata = FALSE, yint = yint,
                        nlab = nlab, X = NULL, fig = initfig,
                        d = NULL, figparams = NULL) #DON'T REALLYU KNOW IF TO FEED OR NOT D

  figparams <- farout$figparams
  figparams$ncolX <- d

  ####################################################
  # calculating pairwise distances, for additional visualization feature

  pw_mdS2 <-function(x, sd, prior, weight) { #pairwise mahalobis squared
    if(! missing(weight)) {
      posid <- (weight > 0)
      if(any(posid)) {
        weight <- sqrt(weight[posid])
        x <- x[posid,  , drop = FALSE] * weight #get only positions non zero positions
      }
      else {
        mat <- outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    p=ncol(t(x))
    n=nrow(t(x))
    pwd=matrix(NA, nrow=n, ncol=n)
    sd=sd[posid]
    for (i in 1:n){
      pwd[,i]=mahalanobis(t(x),t(x)[i,],cov=diag(sd^2))
    }
    pwd
  }

  #pwd=pw_mdS2(xtest, sd, weight=posid)
  ###############################################

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
              ofarness = farout$ofarnes,
              pamrfit = pamrfit, #needed for computations in vcr.pamr.newdata
              pamrfitcv = pamrfitcv, #needed to test in vcr.pamr.newdata whether vcr.pamr.train out is right
              threshold = threshold, #effective threshold used
              ii=ii, #number of threshold selected
              initfig=initfig, #INITFIG ADDED FOR TEST
              posid=posid, #added for test
              sd=sd#pwd=pwd #pairwise matrix distance for mds
              ))
}

vcr.pamr.newdata <- function(newdata, vcr.pamr.train.out, threshold=vcr.pamr.train.out$threshold, prior=NULL, threshold.scale=NULL){
  #
  # Prepares graphical display of new data fitted by a neural
  # net that was modeled on the training data, using the output
  # of vcr.pamr.train() on the training data.
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
  # subsetting to same subset of gene

  #newdata$x=newdata$x[vcr.pamr.train.out$pamrfit$gene.subset , ] #subsetting not needed, not done in pamr.predict

  Xnew <- as.matrix(t(newdata$x))

  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  d <- vcr.pamr.train.out$figparams$ncolX
  if (ncol(Xnew) != d) {
    stop(paste0("newdata$x should have ", d,
                " columns, like the training data$x."))
  }
  if (sum(is.na(as.vector(Xnew))) > 0) {
    stop("The coordinate matrix newdata$x contains NA's.")
  }
  levels <- vcr.pamr.train.out$levels # as in training data
  nlab   <- length(levels) # number of classes

  ynew=newdata$y

  if (is.null(ynew)) ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = FALSE, levels)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv   # INDiceS of cases we can Visualize,
  #                         # can even be empty, is OK.
  nvis    <- length(indsv)   # can be zero.
  yintnew <- rep(NA, n)       # needed to overwrite "new" levels.
  yintnew[indsv] <- lab2int(ynew[indsv])
  #
  #
  # computing posterior using function of pamr called pamr.predict
  probs=pamr.predict(vcr.pamr.train.out$pamrfit, newx=newdata$x, threshold=threshold, type = c("posterior"))

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
  soft.shrink <-function(delta, threshold) {
    dif <- abs(delta) - threshold
    delta <- sign(delta) * dif * (dif > 0)
    nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
    attr(delta, "nonzero") <- nonzero
    delta
  }

  diag.disc.original <-function(x, centroids, prior, weight) {
    ### Computes the class discriminant functions assuming scaled x and centroids
    if(! missing(weight)) {
      posid <- (weight > 0)
      if(any(posid)) {
        weight <- sqrt(weight[posid])
        centroids <- centroids[posid,  , drop = FALSE] * weight
        x <- x[posid,  , drop = FALSE] * weight
      }
      else {
        mat <- outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    dd <- t(x) %*% centroids
    dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
    names(dd0) <- NULL
    scale(dd, dd0, FALSE)
  }

  mdS2 <-function(x, centroids, prior, weight) {
    ### Computes the class discriminant functions assuming scaled x and centroids
    if(! missing(weight)) {
      posid <- (weight > 0)
      if(any(posid)) {
        weight <- sqrt(weight[posid])
        centroids <- centroids[posid,  , drop = FALSE] * weight
        x <- x[posid,  , drop = FALSE] * weight
      }
      else {
        mat <- outer(rep(1, ncol(x)), log(prior), "*")
        dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
        return(mat)
      }
    }
    p=ncol(t(x))
    n=nrow(t(x))
    k=ncol(centroids)
    dd=matrix(NA, nrow=n, ncol=k)
    for (k in 1:ncol(centroids)){
      dd[,k]=mahalanobis(t(x),centroids[,k],cov=diag(p))
    }
    dd
  }
  ###################

  #actually reconstrunctanting dd
  #getting quantities needed from pamrfit
  pamrfit=vcr.pamr.train.out$pamrfit
  ii=vcr.pamr.train.out$ii

  if (!is.null(vcr.pamr.train.out$pamrfitcv)){
    stop("The vcr.pamr.train.out feeded is calculated on pamr.cv object and not on pamr.train object")
  }


  centroids=pamrfit$centroids #centroids per variable per class
  centroid.overall=pamrfit$centroid.overall
  sd=pamrfit$sd
  threshold=pamrfit$threshold[ii]
  se.scale=pamrfit$se.scale
  threshold.scale=pamrfit$threshold.scale
  prior=pamrfit$prior
  nonzero=pamrfit$nonzero[ii]
  K=length(prior)

  #getting deltas (dik)
  delta <- (centroids - centroid.overall)/sd
  delta <- scale(delta, FALSE, threshold.scale * se.scale) #gives division by mk
  #getting the shrunken ones (d'ik)
  delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
  #getting d'ik*mk
  delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))
  nonzero_check <- attr(delta.shrunk, "nonzero")
  if(!nonzero_check==nonzero){
    stop(nonzero_check)
  }
  posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
  newDistToClass <- mdS2((newdata$x - centroid.overall)/sd, delta.shrunk, prior, weight = posid) #since the fucntion assumes scaled x, xnew is scale before


  if (any(is.na(newDistToClass))) { ##CHECK PUT FOR DEVELOPING REASONS PROBABLY COULD REMOVE
    stop("newDistToClass has NAs")
  }

  newDistToclass=-newDistToClass #minus because wrt to paper is with opposite sign here

  rd=pamrfit$nonzero[ii] #getting reducted dimension

  farout <- compFarness(type = "affine", testdata = TRUE, yint = yintnew,
                        nlab = nlab, X = NULL, fig = newDistToclass,
                        d = NULL, figparams = vcr.pamr.train.out$figparams) #DON'T REALLYU KNOW IF TO FEED OR NOT D



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
              threshold=threshold #effective threshold used in pamr.prediction
              ))
}


