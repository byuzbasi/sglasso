######################################################################
## These functions are minor modifications or directly copied from the
## grpreg package:
## Breheny, P. and Huang, J. (2015) Group descent algorithms for nonconvex penalized linear 
## and logistic regression models with grouped predictors. Statistics and Computing, 25: 173-187.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.
newXG <- function(X, g, m, ncolY, bilevel, standardize = TRUE) {
  # Prepare the predictor matrix and keep the bookkeeping needed to map
  # solver coefficients back to the user's original column order and scale.
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg", call.=FALSE)
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)
  xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  
  # Normalize the supplied group labels before any column filtering/reordering.
  G <- setupG(g, m, bilevel)
  
  # Reconfigure for multiple outcomes, if necessary.
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    G <- multiG(G, ncolY)
  }
  
  #.Call("std_c",X)
  
  if (standardize) {
    # Feature-level standardization is performed before the C++ solver. Constant
    # columns are dropped here and later restored as zero coefficients.
    std <- std_c(X)
    XX <- std[[1]]
    center <- std[[2]]
    scale <- std[[3]]
    nz <- which(scale > 1e-6)                # non-constant columns
    if (length(nz) != ncol(X)) {
      XX <- XX[, nz, drop=FALSE]
      G <- subsetG(G, nz)
    }
  } else {
    XX <- X
    center <- rep(0, ncol(X))
    scale <- rep(1, ncol(X))
    nz <- seq_len(ncol(X))
  }
  
  # Reorder columns so each group is contiguous, matching the block-coordinate
  # updates used by the solver.
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder')) XX <- XX[, attr(G, 'ord')]
  
  # Group-level standardization uses an orthogonalized design for group lasso
  # fits; bilevel fits keep the feature-level representation.
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group")
  } else {
    g <- as.integer(G)
  }
  
  # Set group multiplier if missing.
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
  }
  
  # Return both transformed design data and inverse-transform metadata.
  return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'), names=xnames,
              center=center, scale=scale, nz=nz))
}

# Convert user group labels to compact integer labels and validate optional
# group multipliers.
setupG <- function(group, m, bilevel) {
  gf <- factor(group)
  if (any(levels(gf)=='0')) {
    g <- as.integer(gf) - 1
    lev <- levels(gf)[levels(gf)!='0']
  } else {
    g <- as.integer(gf)
    lev <- levels(gf)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    #if (all.equal(sort(names(m)), sort(group)))
    TRY <- try(as.integer(group)==g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call.=FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call.=FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  structure(g, levels=lev, m=m)
}

# Drop group labels/multipliers for columns removed during standardization and
# renumber the remaining groups.
subsetG <- function(g, nz) {
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz]
  dropped <- setdiff(g, new)
  if (length(dropped)) {
    lev <- lev[-dropped]
    m <- m[-dropped]
    gf <- factor(new)
    new <- as.integer(gf) - 1*any(levels(gf)=='0')
  }
  structure(new, levels=lev, m=m)
}

# Put columns belonging to the same group next to one another and store the
# inverse order so coefficients can be returned in the original order.
reorderG <- function(g, m, bilevel) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g==0)) {
    g <- as.integer(relevel(factor(g), "0"))-1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels=lev, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}

# Validate and center the response. The Gaussian solver works with centered y
# and stores the mean so the intercept can be reconstructed later.
newY <- function(y, family) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }
  
  # Convert fuzzy binomial data
  if (family=="binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y=", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
      if (d[2] > 1) attr(y, "dim") <- d
    }
  }
  
  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning=function(w) {stop("y must be numeric or able to be coerced to numeric", call.=FALSE)})
  }
  if (any(is.na(y))) stop("Missing data (NA's) detected in outcome y.  You must eliminate missing data (e.g., by removing cases or imputation) before passing y to grpreg", call.=FALSE)
  
  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y))) paste("Y", 1:ncol(y), sep="")
    }
    attributes(y) <- NULL
  }
  
  if (family=="gaussian") {
    meanY <- mean(y)
    y <- y - meanY
    attr(y, "mean") <- meanY
  }
  attr(y, "m") <- d[2]
  y
}

# Orthogonalize each penalized group with an SVD transform. This improves the
# group update geometry and records the transform for coefficient recovery.
orthogonalize <- function(X, group) {
  out <- orthogonalize_c(X, as.integer(group))
  XX <- out$X
  attr(XX, "T") <- out$T
  attr(XX, "group") <- out$group
  XX
}

# Map coefficients from the orthogonalized basis back to the standardized
# feature basis.
unorthogonalize <- function(b, XX, group, intercept=TRUE, Tmat = NULL) {
  ind <- !sapply(attr(XX, "T"), is.null)
  if (!is.null(Tmat)) {
    T <- Tmat
    if (intercept) {
      ind0 <- c(1, 1+which(group==0))
      val <- Matrix::as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
    } else if (sum(group==0)) {
      ind0 <- which(group==0)
      val <- Matrix::as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
    } else {
      val <- as.matrix(T %*% b)
    }
    return(val)
  }

  T_list <- attr(XX, "T")
  group0 <- which(group == 0)
  offset <- if (intercept) 1L else 0L
  if (intercept) {
    ind0 <- c(1, 1+which(group==0))
    penalized_b <- b[-ind0, , drop=FALSE]
    val <- matrix(0, nrow = length(ind0) + sum(vapply(T_list[ind], nrow, integer(1))), ncol = ncol(b))
    val[seq_along(ind0), ] <- b[ind0, , drop=FALSE]
  } else if (sum(group==0)) {
    ind0 <- which(group==0)
    penalized_b <- b[-ind0, , drop=FALSE]
    val <- matrix(0, nrow = length(ind0) + sum(vapply(T_list[ind], nrow, integer(1))), ncol = ncol(b))
    val[seq_along(ind0), ] <- b[ind0, , drop=FALSE]
  } else {
    ind0 <- integer(0)
    penalized_b <- b
    val <- matrix(0, nrow = sum(vapply(T_list[ind], nrow, integer(1))), ncol = ncol(b))
  }

  src_start <- 1L
  dst_start <- length(ind0) + 1L
  for (Tg in T_list[ind]) {
    src_end <- src_start + ncol(Tg) - 1L
    dst_end <- dst_start + nrow(Tg) - 1L
    val[dst_start:dst_end, ] <- Tg %*% penalized_b[src_start:src_end, , drop=FALSE]
    src_start <- src_end + 1L
    dst_start <- dst_end + 1L
  }
  val
}

# Undo feature-level standardization and reconstruct the intercept on the
# original data scale.
unstandardize <- function(b, XG) {
  beta <- matrix(0, nrow=1+length(XG$scale), ncol=ncol(b))
  beta[1 + XG$nz,] <- b[-1,] / XG$scale[XG$nz]
  beta[1,] <- b[1,] - crossprod(XG$center, beta[-1, , drop=FALSE])
  beta
}
