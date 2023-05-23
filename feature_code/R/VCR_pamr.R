######################### (code that will not be needed in case of classmap merge)
#to import functions needed for VCR_pamr (especially for checkLabels, computeFarness..)
library(cellWise) #for transfo function used in Comp fareness in VCR_auxillaryFunctions.R
source("R_classmap_package_full/R/VCR_auxiliaryFunctions.R") #importing auxillary functions needed
                                                             # this script is available in classmap package
                                                             # so in case of integration of VCR_pamr this import would be useless


library(pamr) #to get pamr.predict for newdata function, in package can just import that with import pamr.predict from pamr, also to get softshirnki etc
#########################

# functions here devised in same style as they are for other classifiers in classmap package
# particularly here of inspiration was vcr.forest.train that takes in a forestfit object of particular package
# here we take in pamr.train object form pamr package

vcr.pamr.train <- function(data, pamrfit, pamrfitcv=NULL, threshold) {
  #
  # Using the outputs of just pamr.train (or also in addition pamr.cv) for classification
  # applied to the training/cv data, this function prepares for
  # graphical displays.
  #
  #
  # putted pamrfit object just like vcr.forest takes forest fit
  #
  #
  # Arguments:
  #   data      : exactly the same input data used in pamr.train that is a
  #               list with components. -x an expression genes in the rows,
  #               samples in the columns). y- a factor with the class labels
  #               for each sample. Here data$y should be a factor and possibly you should be factorize and the begininning
  #               of the processing because the levels should be in the same order as
  #               used by pamr.train. Additional components: -genenames, a vector
  #               of gene names, and -geneid a vector of gene identifiers.
  #   pamrfit   : the result of a call to pamr.train.
  #   pamrfitcv : optional. the result of a call to pamr.cv, then pamrfit should be
  #               the same as the one given to pamr.cv and the crossvalidated results
  #               are considered
  #   threshold : the desired threshold value at which to evaluate. Then is taken the first theshold
  #               value in the evaluated grid by pamr greater then theee threshold inserted
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

  if (!is.factor(data$y)) {
    stop("data$y is not a factor")
  }
  y=data$y #FACTORIZE IT OR NOT?

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
  # Build D(x,g) taking as reference discriminant score (and then do as classmap do for DA beacuse pamr can be seen as particular case of DA)
  #
  #

  #### Auxillary functions needed: ####

  #function that allows shrinkcage #internal function of pamr package
  soft.shrink <-function(delta, threshold) {
    dif <- abs(delta) - threshold
    delta <- sign(delta) * dif * (dif > 0)
    nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
    attr(delta, "nonzero") <- nonzero
    delta
  }

  #function that calculates diagonal discriminant used to then assign and
  #compute posteriors (here only kept to have inspiration to compute measure D(i,g) and build function just below)
  #internal function of pamr package
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

  #function actually used to compute D(i,g) distance of an observation to the class
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

  #actually getting all the quantities needed to compute D(i,g)
  #similar to what pamr.predict does inside pamr package

  norm.cen <- pamrfit$norm.cen  #to handle if hetero is specified in pamr fitting
  if(!is.null(norm.cen)) {
    data$x <- abs(t(scale(t(data$x), center = norm.cen, scale = FALSE)))
  }

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

  #handling sign contrast
  if(pamrfit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(pamrfit$sign.contrast=="negative"){delta <- delta*(delta<0)}

  #getting the shrunken ones (d'ik)
  delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
  #getting d'ik*mk
  delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))


  nonzero_check <- attr(delta.shrunk, "nonzero") #check added for debugging purposing
  if(!nonzero_check==nonzero){
    stop(nonzero_check)
  }

  posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0 #to store non zero gene/variable positions

  distToClass <- mdS2((data$x - centroid.overall)/sd, delta.shrunk, prior, weight = posid)

  } else { #it means we have a pamr.cv, so different process to reconstruct dd

    distToClass=matrix(NA, nrow=n , ncol=nlab) #defining the matrix contain dd (n x k)
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

      #handling sign contrast
      if(pamrfit$sign.contrast=="positive"){delta <- delta*(delta>0)}
      if(pamrfit$sign.contrast=="negative"){delta <- delta*(delta<0)}

      #getting deltas (dik)
      delta <- (centroids - centroid.overall)/sd
      delta <- scale(delta, FALSE, threshold.scale * se.scale) #gives division by mk
      #getting the shrunken ones (d'ik)
      delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
      #getting d'ik*mk
      delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))

      nonzero_check <- attr(delta.shrunk, "nonzero") #check for debug
      if(!nonzero_check==nonzero){
        stop(nonzero_check)
      }

      posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
      xfold<-data$x[,pamrfitcv$folds[[nf]]] #I want dd only in the obvs considered as test in that fold
                                         # this is done following the logic inside nsccv ( called in pamr.cv) where for each fold nsc is called with x (obvs in train) and xtest the obvs considered as test
                                         # then in nsc xtest in called in dd (as here below)
      distToClass[pamrfitcv$folds[[nf]],] <- mdS2((xfold - centroid.overall)/sd, delta.shrunk, prior, weight = posid)
    }
  }

  if (any(is.na(distToClass))) { # check for debug
    stop("Calculated distToClass has NAs")
  }

  distToClass<-(-distToClass) #minus because wrt to paper sense is with opposite sign here

  rd=pamrfit$nonzero[ii] #getting reduced dimension

  farout <- compFarness(type = "affine", testdata = FALSE, yint = yint, #with affine we can feed our matrix D(i,g) and estimation process according to paper is brought up
                        nlab = nlab, X = NULL, fig = distToClass,
                        d = NULL, figparams = NULL)

  figparams <- farout$figparams
  figparams$ncolX <- d

  ####################################################
  # calculating pairwise distances, for additional visualization feature MDScolorScape
  # thios could probably slow in computation with low thresholds (all genes)

  pw_mdS2 <-function(x, sd, prior, weight) { #pairwise mahalanbis squared
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

  #pwd=pw_mdS2(xtest, sd, weight=posid) #for now it is inactive

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
              ofarness = farout$ofarness,
              pamrfit = pamrfit, #needed for computations in vcr.pamr.newdata
              pamrfitcv = pamrfitcv, #needed to test in vcr.pamr.newdata whether vcr.pamr.train out is right
              threshold = threshold, #effective threshold used
              ii=ii #number of threshold selected
              #distToClass=distToClass, #added for debug
              #posid=posid, #posid and sd could be added for outside pairwise dissimilarity computation
              #sd=sd
              ))
}

vcr.pamr.newdata <- function(newdata, vcr.pamr.train.out, prior=NULL){ #threshold scale and specific threshold (in pamr.predict available) not feeded here because then it would imply to recalculate all figparams, in vcrpamrout figparams are computed on the specific threshold
  #
  # Prepares graphical display of new data fitted by pamr that was modeled on the training data,
  # using the output of vcr.pamr.train() on the training data and predicting newdata at the same threshold.
  #
  # Arguments:
  #   newdata              :exactly the same input data used in pamr.train that is a
  #                         list with components. -x an expression genes in the rows,
  #                         samples in the columns). y- a factor with the class labels
  #                         for each sample. NB #check label switching # check also what checklabels does for test if it reorders based on the string
  #                        .Additional components: -genenames, a vector
  #                         of gene names, and -geneid a vector of gene identifiers.
  #   vcr.pamr.train.out   :output of vcr.pamr.train on the training data (with pamr.cv=NULL)
  #   prior                :a different prior to be used in the preidiction, as in pamr.predict allows
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
  #   fig       : farness of each object i from each class g.
  #               Always exists.
  #   farness   : farness of each object to its given class, for
  #               objects with non-NA yintnew.
  #   ofarness  : For each object i, its lowest farness to any
  #               class, including its own. Always exists.
  #
  #

  #subsetting to same subset of gene
  #newdata$x=newdata$x[vcr.pamr.train.out$pamrfit$gene.subset , ] #subsetting not needed, not done in pamr.predict, should be already inherited thogh pamrfit from vcr.pamr.out

  Xnew <- as.matrix(t(newdata$x))

  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  n <- nrow(Xnew)
  d <- vcr.pamr.train.out$figparams$ncolX
  if (ncol(Xnew) != d) {
    stop(paste0("newdata$x should have ", d,
                " columns, like the training data($x) feeded in vcr.pamr.train"))
  }
  if (sum(is.na(as.vector(Xnew))) > 0) {
    stop("The coordinate matrix newdata$x contains NA's.")
  }
  levels <- vcr.pamr.train.out$levels # as in training data
  nlab   <- length(levels) # number of classes

  if (!is.factor(newdata$y)){
    stop("newdata$y is not a factor. Remember it is better to factorize all together data$y before beginning of pipeline and before splitting to test train.")
  }
  ynew=newdata$y

  # checking labels and producing related quantities

  if (is.null(ynew)) ynew <- factor(rep(NA, n), levels = levels)
  checked <- checkLabels(ynew, n, training = FALSE, levels)
  lab2int <- checked$lab2int # is a function
  indsv   <- checked$indsv   # INDiceS of cases we can Visualize,
  #                         # can even be empty, is OK.
  nvis    <- length(indsv)   # can be zero.
  yintnew <- rep(NA, n)       # needed to overwrite "new" levels.
  yintnew[indsv] <- lab2int(ynew[indsv])
  #
  #getting the threshold the effective one used in vcr.pamr.train (since here after modle fitted we can predicted on any given punctual threshold value)
  threshold=vcr.pamr.train.out$threshold
  #
  #
  # computing posterior using function of pamr called pamr.predict
  if (is.null(prior)) { #if prior not specified we do noot feed it, leve default
    probs=pamr.predict(vcr.pamr.train.out$pamrfit, newx=newdata$x, threshold=threshold, type = c("posterior"))
  }
  else { #prior is specified
    probs=pamr.predict(vcr.pamr.train.out$pamrfit, newx=newdata$x, threshold=threshold, type = c("posterior"), prior=prior)
  }

  #internal check of probs matrix for debugging
  if (length(dim(probs)) != 2) stop("probs should be a matrix.")
  if (ncol(probs) == 1) probs <- t(probs) # if new data is 1 object
  if (ncol(probs) != nlab) stop(paste0(
    "The matrix probs should have ", nlab, " columns"))
  if (any(is.na(probs))) stop("probs should not have any NA's.")

  # Compute prediction for all objects in the new data:
  #
  predint <- apply(probs[, , drop = FALSE], 1, which.max)
  #
  #internal check for debug if max posterior is consistent with prediction
  if (!identical(as.numeric(predint),as.numeric(pamr.predict(vcr.pamr.train.out$pamrfit, newx=newdata$x, threshold=threshold, type = c("class"))))) {
    stop("predint (by max posterior) and class prediction (from pamr predict) do not match")
  }

  # Compute PAC for all objects with available ynew:
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
  #
  #
  #
  # Compute farness:
  #
  # Build D(x,g) taking as reference discriminant score (and then do as classmap do for DA beacuse pamr can be seen as particular case of DA)
  #

  #### Auxillary functions needed: ###############

  #function that allows shrinkcage #internal function of pamr package
  soft.shrink <-function(delta, threshold) {
    dif <- abs(delta) - threshold
    delta <- sign(delta) * dif * (dif > 0)
    nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
    attr(delta, "nonzero") <- nonzero
    delta
  }

  #function that calculates diagonal discriminant used to then assign and
  #compute posteriors (here only kept to have inspiration to compute measure D(i,g) and build function just below)
  #internal function of pamr package
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

  #function actually used to compute D(i,g) distance of an observation to the class
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

  ################### end of auxillary functions needed

  #actually getting all the quantities needed to compute D(i,g)
  #similar to what pamr.predict does inside pamr package

  pamrfit=vcr.pamr.train.out$pamrfit #setting pamrfit object
  ii=vcr.pamr.train.out$ii

  norm.cen <- pamrfit$norm.cen  #to handle if hetero is specified in pamr fitting
  if(!is.null(norm.cen)) {
    data$x <- abs(t(scale(t(data$x), center = norm.cen, scale = FALSE)))
  }

  #check that vcr.pamr.out is done on pamr.train and not on pamr.cv (could also make separate vcr.pamr.cv function)
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

  #handling sign contrast
  if(pamrfit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(pamrfit$sign.contrast=="negative"){delta <- delta*(delta<0)}

  #getting the shrunken ones (d'ik)
  delta.shrunk=soft.shrink(delta,threshold) #we have a problem here, all zero
  #getting d'ik*mk
  delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * se.scale))

  nonzero_check <- attr(delta.shrunk, "nonzero") #check for debug reasons
  if(!nonzero_check==nonzero){
    stop(nonzero_check)
  }

  posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
  newDistToClass <- mdS2((newdata$x - centroid.overall)/sd, delta.shrunk, prior, weight = posid) #since the fucntion assumes scaled x, xnew is scale before

   if (any(is.na(newDistToClass))) { #check for debug reasons
    stop("newDistToClass has NAs")
  }

  newDistToclass=-newDistToClass #minus because wrt to paper is with opposite sign here

  rd=pamrfit$nonzero[ii] #getting reduced dimension

  farout <- compFarness(type = "affine", testdata = TRUE, yint = yintnew,
                        nlab = nlab, X = NULL, fig = newDistToclass,
                        d = NULL, figparams = vcr.pamr.train.out$figparams)



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


