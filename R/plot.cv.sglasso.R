#' @export
plot.cv.sglasso <- function(
    x,
    log.l = TRUE,
    type.tun = c("lambda", "d"),
    type = c("cve"),
    selected = TRUE,
    vertical.line = TRUE,
    col = "red",
    ...
) {
  
  type <- match.arg(type)
  type.tun <- match.arg(type.tun, c("lambda", "d"))
  
  new.args <- list(...)
  
  if (type.tun == "lambda") {
    
    l <- x$lambda
    
    if (log.l) {
      l <- log(l)
      xlab <- expression(log(lambda))
    } else {
      xlab <- expression(lambda)
    }
    
    cve_use  <- x$cve[, x$min_ind[2]]
    cvse_use <- x$cvse[, x$min_ind[2]]
    
    L <- cve_use - cvse_use
    U <- cve_use + cvse_use
    y <- cve_use
    
    ind <- is.finite(l[seq_along(cve_use)])
    ylim <- range(c(L[ind], U[ind]), na.rm = TRUE)
    aind <- ((U - L) / diff(ylim) > 1e-3) & ind
    
    plot.args <- list(
      x = l[ind],
      y = y[ind],
      ylim = ylim,
      xlab = xlab,
      ylab = "Cross-validation error",
      type = "n",
      xlim = rev(range(l[ind], na.rm = TRUE)),
      las = 1,
      bty = "n"
    )
    
    if (length(new.args)) {
      plot.args[names(new.args)] <- new.args
    }
    
    do.call("plot", plot.args)
    
    if (vertical.line) {
      abline(v = l[x$min_ind[1]], lty = 2, lwd = 0.5)
    }
    
    suppressWarnings(
      arrows(
        x0 = l[aind],
        x1 = l[aind],
        y0 = L[aind],
        y1 = U[aind],
        code = 3,
        angle = 90,
        col = "gray80",
        length = 0.05
      )
    )
    
    points(l[ind], y[ind], col = col, pch = 19, cex = 0.5)
    
    if (selected) {
      pred_ns <- predict(
        x$fit,
        lambda = x$lambda,
        d = x$d.min,
        type = "groups"
      )
      
      if (class(pred_ns)[1] != "list") {
        n.s <- rep(dim(pred_ns)[1], length(x$lambda))
      } else {
        n.s <- sapply(pred_ns, length)
      }
      
      axis(
        3,
        at = l,
        labels = n.s,
        tick = FALSE,
        line = -0.5,
        cex.axis = 0.8
      )
      
      mtext("Groups selected", cex = 0.9, line = 1.5)
    }
  } else {
    
    d.val <- x$d
    xlab <- expression(d)
    
    cve_use  <- x$cve[x$min_ind[1], ]
    cvse_use <- x$cvse[x$min_ind[1], ]
    
    L <- cve_use - cvse_use
    U <- cve_use + cvse_use
    y <- cve_use
    
    ind <- is.finite(d.val[seq_along(cve_use)])
    ylim <- range(c(L[ind], U[ind]), na.rm = TRUE)
    aind <- ((U - L) / diff(ylim) > 1e-3) & ind
    
    plot.args <- list(
      x = d.val[ind],
      y = y[ind],
      ylim = ylim,
      xlab = xlab,
      ylab = "Cross-validation error",
      type = "l",
      xlim = range(d.val[ind], na.rm = TRUE),
      las = 1,
      bty = "n"
    )
    
    if (length(new.args)) {
      plot.args[names(new.args)] <- new.args
    }
    
    do.call("plot", plot.args)
    
    if (vertical.line) {
      abline(v = d.val[x$min_ind[2]], lty = 2, lwd = 0.5)
    }
    
    suppressWarnings(
      arrows(
        x0 = d.val[aind],
        x1 = d.val[aind],
        y0 = L[aind],
        y1 = U[aind],
        code = 3,
        angle = 90,
        col = "gray80",
        length = 0.05
      )
    )
    
    points(d.val[ind], y[ind], col = col, pch = 19, cex = 0.5)
    
    if (selected) {
      pred_ns <- predict(
        x$fit,
        lambda = x$lambda.min,
        d = x$d,
        type = "groups"
      )
      
      if (class(pred_ns)[1] != "list") {
        n.s <- rep(dim(pred_ns)[1], length(x$d))
      } else {
        n.s <- sapply(pred_ns, length)
      }
      
      axis(
        3,
        at = d.val,
        labels = n.s,
        tick = FALSE,
        line = -0.5,
        cex.axis = 0.8
      )
      
      mtext("Groups selected", cex = 0.9, line = 1.5)
    }
  }
}