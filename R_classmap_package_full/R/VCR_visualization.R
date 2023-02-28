

confmat.vcr <- function(vcrout, cutoff = 0.99,
                        showClassNumbers = FALSE,
                        showOutliers = TRUE, silent = FALSE) {
  # Builds a confusion matrix from a vcr object.
  #
  # Arguments:
  # vcrout           : output of vcr.*.train() or vcr.*.newdata().
  # cutoff           : the objects with overall farness
  #                    ofarness > cutoff are flagged as outliers.
  # showClassNumbers : if TRUE, the row and column names are the
  #                    number of each level instead of the level
  #                    itself. Useful for long level names.
  # showOutliers     : if TRUE and some points were flagged as
  #                    outliers, it adds an extra column on the
  #                    right of the confusion matrix for these
  #                    outliers, with label "outl".
  # silent           : if FALSE, the confusion matrix and
  #                    accuracy are shown on the screen.
  #
  # Returns:
  # confMat          : the confusion matrix.
  #
  levels <- vcrout$levels
  nlab   <- length(levels)
  if (is.null(vcrout$yintnew)) {
    if (is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew")
    }
    given <- vcrout$yint
    training <- TRUE
  } else {
    given <- vcrout$yintnew
    training <- FALSE
  }
  whichyint <- if (training) {"yint"} else {"yintnew"}
  indsv <- which(!is.na(given))
  if (length(indsv) == 0) {
    stop(paste0(whichyint, " has no available values, ",
                " so no confusion matrix can be made."))
  }
  given <- given[indsv]
  if (sum(!(given %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ", whichyint, "[indsv] lie in 1:", nlab))
    }
  avobs <- sort(unique(given)) # available given values
  predicted <- vcrout$predint[indsv] # corresponding predictions
  if (sum(is.na(predicted)) > 0) {
    stop("\n predint[indsv] has missing values.")}
  if (sum(!(predicted %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all predint[indsv] lie in 1:", nlab))
    }
  # compute the accuracy before taking outliers into account:
  accuracy <- sum(given == predicted) / length(indsv)
  myrownames <- levels[avobs]
  ofarness <- vcrout$ofarness[indsv]
  outliers <- (ofarness > cutoff)
  if (showOutliers & sum(outliers) > 0) {
    predicted[which(outliers)] <- nlab + 1
    confMat <- table(given, predicted)
    mycolints <- sort(unique(predicted))
    mycolnames <- c(levels, "outl")[mycolints]
    mycolints[which(mycolints == (nlab + 1))] <- "outl"
  } else {
    confMat <- table(given, predicted)
    mycolints <- sort(unique(predicted))
    mycolnames <- levels[mycolints]
  }
  if (showClassNumbers) {
    rownames(confMat) <- avobs
    colnames(confMat) <- mycolints
  } else {
    rownames(confMat) <- myrownames
    colnames(confMat) <- mycolnames
  }
  if (silent != TRUE) {
    wnq(paste0("\nConfusion matrix:"))
    pnq(confMat)
    wnq(paste0("\nThe accuracy is ", round(100 * accuracy, 2), "%."))
  }
  invisible(confMat)
}



stackedplot <- function(vcrout, cutoff = 0.99, classCols = NULL,
                        classLabels = NULL, separSize = 1,
                        minSize = 1.5, showOutliers = TRUE,
                        showLegend = FALSE, main = NULL,
                        htitle = NULL, vtitle = NULL)
{
  # Argument added:
  # main              : title for the plot
  #
  if (is.null(vcrout$yintnew)) {
    if (is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew")
    }
    whichdata <- "[training]"
  }
  else {
    whichdata <- "[newdata]"
  }
  if (is.null(htitle))
    htitle <- paste0(whichdata, " given class")
  if (is.null(vtitle))
    vtitle <- paste0(whichdata, " predicted class")
  nlab <- length(vcrout$levels)
  if (is.null(classLabels)) {
    lvls <- vcrout$levels
    if ("outl" %in% lvls)
      stop("Please rename the level currently called 'outl'.")
  }
  else {
    if (!is.vector(classLabels)) {
      stop("\n classLabels should be a vector")
    }
    lvls <- classLabels
    if ("outl" %in% lvls)
      stop("Do not use 'outl' as a class label.")
    if (length(lvls) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ", nlab, "."))
    }
  }
  if (is.null(classCols)) {
    mosaicColors <- rainbow(nlab)
  }
  else {
    if (!is.vector(classCols)) {
      stop("\n classCols should be a vector")
    }
    if (length(classCols) < nlab) {
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ", nlab,
                  "."))
    }
    mosaicColors <- classCols[seq_len(nlab)]
  }
  confusi <- confmat.vcr(vcrout, cutoff = cutoff, showOutliers = showOutliers,
                         showClassNumbers = TRUE, silent = TRUE)
  pconfusi <- confusi / rowSums(confusi) * 100
  ayint <- as.numeric(rownames(confusi))
  nlaba <- length(ayint)
  if (nlaba < nlab) {
    wnq(paste0("\nNot all classes occur in these data.",
               " The classes to plot are:"))
    pnq(ayint)
  }
  colints <- colnames(confusi)
  colints[which(colints == "outl")] <- nlab + 1
  for (j in seq_len(nlaba)) {
    confvec <- pconfusi[j, ]
    tooSmall <- which((confvec > 0) & (confvec < minSize))
    if (length(tooSmall) > 0) {
      confvec[tooSmall] <- minSize
      confvec <- confvec/sum(confvec) * 100
      pconfusi[j, ] <- confvec
    }
  }
  hasOutliers <- 0
  if (showOutliers) {
    if ("outl" %in% colnames(confusi)) {
      lvls <- c(lvls, "outliers")
      hasOutliers <- 1
      mosaicColors <- c(mosaicColors[seq_len(nlab)], "gray20")
    }
  }
  confRow <- rep(0, length(lvls))
  colints <- as.numeric(colints)
  confRow[colints] <- pconfusi[1, ]
  df <- data.frame(g1_fill = factor(lvls[seq_len(nlab + hasOutliers)],
                                    lvls), g1_count = confRow)
  if (nlaba > 1) {
    for (j in 2:nlaba) {
      g <- ayint[j]
      confRow <- rep(0, length(lvls))
      colints <- as.numeric(colints)
      confRow[colints] <- pconfusi[j, ]
      fillname <- paste0("g", j, "_fill")
      perm <- c(g, (seq_len(nlab + hasOutliers))[-g])
      df[[fillname]] <- factor(lvls[perm], lvls[perm])
      countname <- paste0("g", j, "_count")
      df[[countname]] <- confRow[perm]
    }
  }
  widths <- table(vcrout$yint)/length(vcrout$yint)
  pos <- 0.5 * (cumsum(widths) + cumsum(c(0, widths[-length(widths)])))
  gg <- ggplot(df)
  for (i in seq_len(nlaba)) {
    fillname <- paste0("g", i, "_fill")
    countname <- paste0("g", i, "_count")
    gg <- gg + geom_bar(mapping = aes_string(x = pos[i], y = countname,
                                             fill = fillname),
                        width = widths[i], stat = "identity",
                        size = separSize, colour = "white",
                        show.legend = showLegend)
  }
  gg <- gg + scale_x_continuous(htitle, breaks = pos,
                                labels = lvls[ayint],
                                expand = c(0.02, 0))
  gg <- gg + scale_y_continuous(vtitle, trans = "reverse",
                                expand = c(0.02, 0))
  gg <- gg + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank())
  if (showLegend) {
    gg <- gg + labs(fill = "Classes")
  }
  gg <- gg + theme(axis.ticks.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())
  gg <- gg + scale_fill_manual(values = mosaicColors)
  if (!is.null(main)) {
    gg <- gg + labs(title = paste0(main))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5)) # center title
  }
  return(gg)
}




classmap <- function(vcrout, whichclass, classLabels = NULL,
                     classCols = NULL, main = NULL, cutoff = 0.99,
                     plotcutoff = TRUE, identify = FALSE, cex = 1,
                     cex.main = 1.2, cex.lab = NULL, cex.axis = NULL,
                     opacity = 1, squareplot = TRUE,
                     maxprob = NULL, maxfactor = NULL) {
  #
  # Draws class maps to visualize classification results.
  #
  # Arguments:
  # vcrout       : the result of a function vcr.*.*()
  # whichclass   : the class to be displayed. Required. The colors
  #                of the points are those of their predicted
  #                labels. Required.
  # classLabels  : the labels (levels) of the classes. If NULL, they
  #                are taken from vcrout.
  # classCols    : a list of colors for the class labels 1,... There
  #                should be at least as many as there are levels.
  #                If NULL the classCols are taken as 2,...
  # main         : title for the plot.
  # cutoff       : objects with overall farness > cutoff will be
  #                flagged as outliers and plotted with a black
  #                border around them.
  # plotcutoff   : if TRUE, a dashed vertical line is drawn at
  #                the cutoff.
  # identify     : if TRUE, left-click on a point to get its number,
  #                then ESC to exit.
  # cex          : passed on to plot().
  # cex.main     : same, for title.
  # cex.lab      : same, for labels on horizontal and vertical axes.
  # cex.axis     : same, for axes.
  # opacity      : determines opacity of plotted dots.
  #                Value between 0 and 1, where 0 is transparent
  #                and 1 is opaque.
  # squareplot   : if TRUE, the horizontal and vertical axis will
  #                get the same length.
  # maxprob      : draws the farness axis at least upto probability
  #                maxprob.
  #                If NULL, the limits are obtained automatically.
  # maxfactor    : if not NULL, a number slightly higher than 1
  #                to increase the space at the right hand side of
  #                the plot, to make room for marking points.
  #
  # Returns:
  # coordinates  : a matrix with 2 columns containing the
  #                coordinates of the plotted points. The
  #                first coordinate is the quantile of the
  #                farness probability.

  # Auxiliary function:
  qfunc <- function(probs) {
    # quantile function of N(0,1) restricted to [0,4].
    if (min(probs) < 0) stop("There are negative probs")
    if (max(probs) > 1) stop("there are probs > 1")
    qnorm(probs * (pnorm(4) - 0.5) + 0.5)
  }

  # The main function starts here:
  nlab <- length(vcrout$levels)
  if (is.null(classLabels)) { # no classLabels are given
    levels <- vcrout$levels
  } else {# classLabels are given
    if (!is.vector(classLabels)) {
      stop("\n classlabels should be a vector")
    }
    levels <- classLabels
    if (length(levels) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ", nlab, "."))
    }
  }
  if (is.null(classCols)) {
    classCols <- seq_len(nlab) + 1 # red, green, blue,...
  } else {
    if (!is.vector(classCols)) {
      stop("\n classCols should be a vector")
    }
    if (length(classCols) < nlab) {
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ", nlab, "."))
    }
    classCols <- classCols[seq_len(nlab)]
  }
  if (is.null(vcrout$yintnew)) {
    if (is.null(vcrout$yint)) stop("vcrout$yint is missing")
    yint <- vcrout$yint
    training <- TRUE
  } else {
    yint <- vcrout$yintnew
    training <- FALSE
  }
  whichdata <- if (training) {"[training]"} else {"[newdata]"}
  whichyint <- if (training) {"yint"} else {"yintnew"}
  indsv <- which(!is.na(yint))
  if (length(indsv) == 0) {
    stop(paste0(whichyint, " has no available values, ",
                " so the plot cannot be made.")) }
  yint     <- vcrout$yint[indsv]
  predint  <- vcrout$predint[indsv]
  PAC      <- vcrout$PAC[indsv]
  ofarness <- vcrout$ofarness[indsv]
  outliers <- (ofarness > cutoff)
  # In the classmap we will plot farness at quantile positions:
  qfarness <- qfunc(vcrout$farness[indsv])
  coordinates <- cbind(qfarness, PAC)
  rownames(coordinates) <- indsv # case numbers in entire set.
  #
  fardashed <- qfunc(cutoff)
  if (!is.null(maxprob)) maxfar <- min(maxprob, 1)
  if (is.null(maxfactor)) maxfactor <- 1
  maxfar <- max(max(qfarness), fardashed)
  if (sum(!(yint %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ", whichyint, " lie in 1:", nlab)) }
  ayint <- sort(unique(yint)) # available yint values
  #
  if (sum(!(predint %in% seq_len(nlab))) > 0) {
    stop(paste0(whichdata, " Not all predint lie in 1:", nlab))
  }
  #
  if (is.numeric(whichclass)) {
    if (!(whichclass %in% seq_len(nlab))) stop(paste0(
      "whichclass must be in 1:", nlab, ", the number of classes."))
  } else {#
    classnumber <- match(whichclass, vcrout$levels)
    if (is.na(classnumber)) stop(paste0(
      'level "', whichclass, '" does not exist.'))
    whichclass <- classnumber
  }
  if (!(whichclass %in% ayint)) stop(paste0(
    "Class number ", whichclass, " with label ",
    levels[whichclass], " has no objects to visualize."))
  if (is.null(main)) {
    main <- paste0(whichdata, " predictions of class ",
                   levels[whichclass])
  }
  clinds   <- which(yint == whichclass) # class indices
  predint  <- predint[clinds]
  casecols <- classCols[predint]
  coordinates <- coordinates[clinds, ]
  outliers <- which(outliers[clinds])
  #
  # Prepare for the plot:
  #
  PAClim <- PACsolid <- c(0, 1)
  mypch <- rep(19, length(predint)) # disk (16 is a bit smaller)
  mypch[outliers] <- 21 # filled circle
  mycex <- rep(cex, length(predint))
  mycex[outliers] <- 1.1*cex
  mycol <- mybg <- adjustcolor(casecols, alpha.f = opacity)
  mycol[outliers] <- "black"
  mycex.main <- 1
  if (!is.null(cex.main)) {mycex.main <- cex.main}
  mycex.lab <- 1
  if (!is.null(cex.lab)) {mycex.lab <- cex.lab}
  mycex.axis <- 1
  if (!is.null(cex.axis)) {mycex.axis <- cex.axis}
  #
  # Make actual plot
  #
  if (squareplot == TRUE) par(pty = "s")
  # makes the axes of the plot equally long
  farlim <- c(0, maxfar * maxfactor)
  plot(0, 0, type = "n", ann = FALSE, axes = FALSE,
       xlim = farlim, ylim = PAClim)
  # makes an empty plot
  title(main = main, line = 1, cex.main = mycex.main)
  title(ylab =  "P[alternative class]", line = 2.3,
        cex.lab = mycex.lab)
  title(xlab = "farness from given class",
        line = 2.3, cex.lab = mycex.lab)
  rect(farlim[1] - farlim[2], -2,farlim[2] * 2, 0.5, col = "lightgray",
       border = FALSE) # adds rectangle
  # col can be modified to other gray intensities
  par(new = TRUE) # overlay next plot on this one
  plot(coordinates, xlab = "", ylab = "", xlim = farlim,
       ylim = PAClim, col = mycol, bg = mybg, pch = mypch,
       cex = mycex, cex.axis = mycex.axis, xaxt = "n")
  probs <- c(0, 0.5, 0.75, 0.9, 0.99, 0.999, 1)
  axis(1, at = qfunc(probs), labels = probs,
       cex.axis = mycex.axis)
  if (plotcutoff) abline(v = fardashed, lty = 4)
  abline(h = PACsolid)
  attr(coordinates, "ids") <- integer(0)
  if (identify) {
    message("Press the escape key to stop identifying.")
    iout <- identify(coordinates, order = TRUE)
    ids <- iout$ind[order(iout$order)]
    if (length(ids) > 0) {
      wnq("Identified point(s): ")
      pnq(ids)
      attr(coordinates, "ids") <- ids
    } else {wnq("No points were identified.")}
  }
  invisible(coordinates)
}



silplot <- function(vcrout, classLabels = NULL, classCols = NULL,
                    showLegend = TRUE, showClassNumbers = FALSE,
                    showCases = FALSE, drawLineAtAverage = FALSE,
                    topdown = TRUE, main = NULL, summary = TRUE)
{
  nlab <- length(vcrout$levels)
  if (is.null(classLabels)) {
    lvls <- vcrout$levels
  }
  else {
    if (!is.vector(classLabels)) {
      stop("\n classLabels should be a vector")
    }
    lvls <- classLabels
    if (length(lvls) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ", nlab, "."))
    }
  }
  if (is.null(classCols)) {
    classCols <- rainbow(nlab)
  }
  else {
    if (!is.vector(classCols)) {
      stop("\n classCols should be a vector")
    }
    if (length(classCols) < nlab) {
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ", nlab,
                  "."))
    }
    classCols <- classCols[seq_len(nlab)]
  }
  if (is.null(vcrout$yintnew)) {
    if (is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew")
    }
    yintv <- vcrout$yint
    training <- TRUE
  }
  else {
    yintv <- vcrout$yintnew
    training <- FALSE
  }
  whichdata <- if (training) {
    "[training]"
  }
  else {
    "[newdata]"
  }
  whichyint <- if (training) {
    "yint"
  }
  else {
    "yintnew"
  }
  indsv <- which(!is.na(yintv))
  if (length(indsv) == 0) {
    stop(paste0(whichyint, " has no available values, ",
                " so no silhouette plot can be made."))
  }
  if (length(indsv) < 2) {
    stop(paste0("At least 2 cases with non-missing ",
                whichyint, " are required."))
  }
  yintv <- yintv[indsv]
  ayint <- sort(unique(yintv))
  if (sum(!(yintv %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ", whichyint, "[indsv] lie in 1:",
                nlab))
  }
  PAC <- vcrout$PAC[indsv]
  if (sum(is.na(PAC)) > 0)
    stop("PAC[indsv] has missing values.")
  si <- 1 - 2 * PAC
  if (is.null(main)) {
    main <- paste0(whichdata, " Silhouette plot of classification")
  }
  df <- as.data.frame(cbind(class = yintv, si, name = indsv))
  if (topdown) {
    df <- df[order(-df$class, df$si), ]
  }
  else {
    df <- df[order(df$class, -df$si), ]
  }
  avswidth <- mean(df$si)
  df$name  <- factor(df$name, levels = df$name)
  df$class <- as.factor(df$class)
  if (showClassNumbers) {
    mapping <- aes_string(x = "name", y = "si",
                          color = "class", fill = "class")
  }
  else {
    df$label <- factor(lvls[df$class], levels = lvls[seq_len(nlab)])
    mapping <- aes_string(x = "name", y = "si",
                          color = "label", fill = "label")
  }
  gg <- ggplot(df, mapping) + geom_bar(stat = "identity",
                                       show.legend = showLegend,
                                       size = 0.05, width = 0.75) +
    labs(y = paste0("Silhouette width s(i)"),
         caption = paste0("\nOverall average silhouette width: ",
                          round(avswidth, 2)), x = "",
         title = paste0(main))
  silwidths_perclass <- round(sapply(lvls,
                                     function(y) mean(df$si[which(df$label == y)])), 2)
  templabels <- c()
  # ?sprintf
  for (i in seq_len(length(unique(df$class)))) {
    sbar <- sprintf("%1.2f", silwidths_perclass[i])
    templabel <- paste('paste(bar(s), " = ", ',
                       '"', sprintf("%1.2f", silwidths_perclass[i]),'"',
                       ', " ", ', '"', lvls[i], '"', ")")
    templabels <- c(templabels, templabel)
  }
  gg <- gg + scale_fill_manual(values = classCols[ayint],
                               aesthetics = c("colour", "fill"),
                               labels = parse(text = templabels),
                               name = "Classes") +
    theme(legend.text.align = 0)
  if (topdown)
    gg <- gg + coord_flip()
  if (drawLineAtAverage) {
    gg <- gg + geom_hline(yintercept = avswidth, linetype = "dashed",
                          color = "red")
  }
  if (!showCases) {
    if (topdown) {
      gg <- gg + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
    }
    else {
      gg <- gg + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }
  }
  else {
    if (topdown) {
      gg <- gg + theme(axis.text.y = element_text(angle = 0))
    }
    else {
      gg <- gg + theme(axis.text.x = element_text(angle = 90))
    }
  }
  gg <- gg + theme(plot.title = element_text(hjust = 0.5))
  gg <- gg + theme(axis.line.x = element_line(size = 0.25))
  gg <- gg + scale_x_discrete(expand = c(0, 0.013 * length(indsv)))
  gg <- gg + scale_y_continuous(expand = c(0, 0), limits = c(-1, 1))
  gg <- gg + theme(panel.background = element_blank())
  gg <- gg + theme(plot.title.position = "plot")
  gg <- gg + theme(axis.title = element_text(size = 10))
  gg <- gg + theme(plot.caption = element_text(size = 10, hjust = 0,
                                               margin = ggplot2::margin(0, 0, 0, 0)))
  if (summary) {
    ave <- tapply(df$si, df$class, mean)
    n <- tapply(df$class, df$class, length)
    sil.sum <- data.frame(classNumber = names(ave), classLabel = lvls[as.numeric(names(ave))],
                          classSize = n, classAveSi = round(ave, 2), stringsAsFactors = TRUE)
    print(sil.sum, row.names = FALSE)
  }
  return(gg)
}



qresplot <- function(PAC, feat, xlab = NULL, xlim = NULL,
                     main = NULL, identify = FALSE, gray = TRUE,
                     opacity = 1, squareplot = FALSE, plotLoess = FALSE,
                     plotErrorBars = FALSE, plotQuantiles = FALSE,
                     grid = NULL, probs = c(0.5, 0.75),
                     cols = NULL, fac = 1, cex = 1,
                     cex.main = 1.2, cex.lab = 1,
                     cex.axis = 1, pch = 19){
  #
  # Draws a quasi residual plot of PAC versus a data feature
  #
  # Arguments:
  # PAC           : vector with the PAC values of a classification,
  #                 typically the $PAC of a function vcr.*.*()
  # feat          : the PAC will be plotted versus this data
  #                 feature. Note that feat does not have to be
  #                 one of the explanatory variables of the model.
  #                 It can be another variable, a combination of
  #                 variables (like a sum or a principal component
  #                 score), the row number of the cases if they
  #                 were recorded succesively, etc.
  # xlab          : label for the horizontal axis, i.e. the name
  #                 of variable feat.
  # xlim          : limits for the horizontal axis. If NULL, the
  #                 range of feat is used.
  # main          : title for the plot.
  # identify      : if TRUE, left-click on a point to get its number,
  #                 then ESC to exit.
  # gray          : logical, if TRUE (the default) the plot region
  #                 where PAC < 0.5 gets a light gray background.
  #                 Points in this region were classified into
  #                 their given class, and the points above this
  #                 region were misclassified.
  # opacity       : determines opacity of plotted dots.
  #                 Value between 0 and 1, where 0 is transparent
  #                 and 1 is opaque.
  # squareplot    : if TRUE, the horizontal and vertical axis will
  #                 get the same length.
  # plotLoess     : if TRUE, a standard loess curve is fitted and
  #                 superimposed on the plot. May not work well if
  #                 feat is discrete with few values.
  #                 At most one of the options plotLoess,
  #                 plotErrorbars, or plotQuantiles can be selected.
  # plotErrorBars : if TRUE, the average PAC and its standard
  #                 error are computed on the intervals of a grid
  #                 (see option grid). Then a red curve connecting
  #                 the averages is plotted, as well as two blue
  #                 curves corresponding to the average plus or
  #                 minus one standard error. At most one of the
  #                 options plotLoess, plotErrorbars, or
  #                 plotQuantiles can be selected.
  # plotQuantiles : if TRUE, one or more quantiles of the PAC
  #                 are computed on the intervals of a grid
  #                 (see option grid). The quantiles correspond
  #                 the probabilities in option probs.
  #                 Then the curves connecting the quantiles
  #                 are plotted. At most one of the options plotLoess,
  #                 plotErrorbars, or plotQuantiles can be selected.
  # grid          : only used when plotErrorBars or plotQuantiles
  #                 are selected. This is a vector with increasing
  #                 feat values, forming the grid. If NULL, the
  #                 grid consists of the minimum and the maximum
  #                 of feat, with 9 equispaced points between them.
  # probs         : only used when plotQuantiles is selected. This
  #                 is a vector with probabilities determining the
  #                 quantiles. If NULL, defaults to c(0.5, 0.75).
  # cols          : only used when plotquantiles is selected.
  #                 A vector with the colors of the quantile curves.
  #                 If NULL the cols are taken as 2,...
  # fac           : only used when ploTloess, plotErrorBars or
  #                 plotQuantiles are selected. A real number to
  #                 multiply the resulting curves. A value fac > 1
  #                 can be useful to better visualize the curves
  #                 when they would be too close to zero.
  #                 By default (fac=1) this is not done.
  # cex           : passed on to plot().
  # cex.main      : same, for title.
  # cex.lab       : same, for labels on horizontal and vertical axes.
  # cex.axis      : same, for axes.
  # pch           : plot character for the points, defaults to 19.

  # Returns:
  # coordinates  : a matrix with 2 columns containing the
  #                coordinates of the plotted points. This makes it
  #                easier to add text next to interesting points.
  #                If identify = T, the attribute ids of coordinates
  #                contains the row numbers of the identified points
  #                in the matrix coordinates.

  if (!is.vector(PAC)) stop(
    "The PAC on the vertical axis should be a vector.")
  n <- length(PAC)
  if (!is.vector(feat)) stop(
    "The feature on the horizontal axis should be a vector.")
  if (length(feat) != n) stop(paste0("The feature should have ",
                                     "length ", n, ", like PAC."))
  if (sum(c(plotLoess, plotErrorBars, plotQuantiles)) > 1) {
    stop(paste0("Only one of the options plotLoess, ",
                "plotErrorBars, or plotQuantiles can be selected.")) }
  if (is.null(main)) {
    mymain <- "quasi residual plot"
  } else {
    mymain <- main
  }
  if (is.null(xlab)) {
    myxlab <- "feature"
  } else {
    myxlab <- xlab
  }
  xran <- range(feat, na.rm = TRUE)
  myxlim <- xran
  if (!is.null(xlim)) {myxlim <- xlim}
  if (is.null(gray)) {gray <- TRUE}
  myopacity <- 1
  if (!is.null(opacity)) {myopacity <- opacity}
  if (squareplot) {par(pty = "s")} else {par(pty = "m")}
  # The "else" part is still missing in classmap,
  # but there the default is sqauareplot=T anyway.
  if (is.null(fac)) {fac <- 1}
  mycex <- 1
  if (!is.null(cex)) {mycex <- cex}
  mycex.main <- 1.2
  if (!is.null(cex.main)) {mycex.main <- cex.main}
  mycex.lab <- 1
  if (!is.null(cex.lab)) {mycex.lab <- cex.lab}
  mycex.axis <- 1
  if (!is.null(cex.axis)) {mycex.axis <- cex.axis}
  mypch <- 19 # disk (16 is a bit smaller)
  if (!is.null(pch)) {mypch <- pch}
  # Start the actual plot:
  plot(0, 0, type = "n", xlim = myxlim, ylim = c(0, 1),
       ann = FALSE, axes = FALSE)
  # makes an empty plot
  title(main = mymain, line = 1, cex.main = cex.main)
  title(ylab = "P[alternative class]", line = 2.3,
        cex.lab = 1)
  title(xlab = myxlab, line = 2.3, cex.lab = 1)
  if (gray) { # adds gray region
    rect(myxlim[1] - myxlim[2], -2, myxlim[2] * 2, 0.5,
         col = "lightgray", border = FALSE)
  }
  par(new = TRUE) # overlay next plot on this one
  mycol <- adjustcolor(1, alpha.f = myopacity)
  coordinates <- cbind(feat, PAC)
  plot(coordinates, xlab = "", ylab = "", xlim = myxlim,
       ylim = c(0, 1), pch = mypch, col = mycol,
       cex = mycex, cex.axis = mycex.axis)
  if (plotLoess) {# add a loess curve
    lofit <- loess(PAC ~ feat)
    grid <- seq(from = myxlim[1], to = myxlim[2], length.out = 201)
    lines(grid, fac * predict(lofit, grid), col = "red", lwd = 2)
  }
  if (plotErrorBars) { # add average plus/minus standard error
    if (is.null(grid)) grid <-
        seq(from = xran[1], to = xran[2], length.out = 11)
    nintervals <- length(grid) - 1
    midpoints <- avers <- sders <- rep(NA, nintervals)
    for (j in seq_len(nintervals)) {
      midpoints[j] <- (grid[j] + grid[j + 1]) / 2
      clinds <- which((grid[j] <= feat) & (feat < grid[j + 1]))
      avers[j] <- mean(PAC[clinds], na.rm = TRUE)
      sders[j] <- sd(PAC[clinds], na.rm = TRUE) / sqrt(length(clinds))
    }
    lines(midpoints, fac * (avers + sders), col = "blue", lwd = 2)
    lines(midpoints, fac * avers, col = "red", lwd = 2)
    lines(midpoints, fac * (avers - sders), col = "blue", lwd = 2)
  }
  if (plotQuantiles) { # add quantile-based curves
    if (is.null(grid)) grid <-
        seq(from = xran[1], to = xran[2], length.out = 11)
    nintervals <- length(grid) - 1
    midpoints <- rep(NA, nintervals)
    if (is.null(probs)) probs <- c(0.5, 0.75)
    nprobs <- length(probs)
    quants <- matrix(NA, ncol = nintervals, nrow = nprobs)
    for (j in seq_len(nintervals)) {
      midpoints[j] <- (grid[j] + grid[j + 1]) / 2
      clinds <- which((grid[j] <= feat) & (feat < grid[j + 1]))
      quants[, j] <- quantile(PAC[clinds], probs)
    }
    if (is.null(cols)) {
      cols <- seq_len(nprobs) + 1 # red, green, blue,...
    } else {
      if (!is.vector(cols)) {
        stop("\n cols should be a vector")
      }
      if (length(cols) < nprobs) {
        stop(paste0("\n The number of cols should be at",
                    " least length(probs) = ", nprobs, "."))
      }
      cols <- cols[seq_len(nprobs)]
    }
    for (k in seq_len(nprobs)) {
      lines(midpoints, fac * quants[k, ], col = cols[k], lwd = 2)
    }
  }
  attr(coordinates, "ids") <- integer(0)
  if (identify) {
    message("Press the escape key to stop identifying.")
    iout <- identify(coordinates, order = TRUE)
    ids <- iout$ind[order(iout$order)]
    if (length(ids) > 0) {
      write(noquote("Identified point(s): "))
      print(noquote(ids))
      attr(coordinates, "ids") <- ids
    } else {
      write(noquote("No points were identified."))
    }
  }
  invisible(coordinates)
}
