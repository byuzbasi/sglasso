#' @export
plot.cv.sglasso <- function(x,
                            type.tun = c("lambda", "d"),
                            log.l = TRUE,
                            col = "red",
                            selected = TRUE,
                            ...) {

  type.tun <- match.arg(type.tun)

  if (type.tun == "d") {
    d <- x$d
    cve <- x$cve[x$min_ind[1], ]
    cvse <- x$cvse[x$min_ind[1], ]

    L <- cve - cvse
    U <- cve + cvse

    ind <- is.finite(d) & seq_along(cve) <= length(d)
    ylim <- range(c(L[ind], U[ind]), na.rm = TRUE)

    plot(
      d[ind],
      cve[ind],
      type = "n",
      xlab = "d",
      ylab = "CV Error",
      ylim = ylim,
      las = 1,
      bty = "n"
    )

    abline(v = d[x$min_ind[2]], lty = 2, col = "black")

    arrows(
      d[ind],
      L[ind],
      d[ind],
      U[ind],
      code = 3,
      angle = 90,
      col = "gray80",
      length = 0.05
    )

    points(d[ind], cve[ind], col = col, pch = 19, cex = 0.6)
    return(invisible())
  }
  
  # Plot the lambda path on the same scale used for model selection.
  l <- x$lambda
  
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
  }
  
  # The d dimension is fixed at the cross-validated optimum for this curve.
  cve <- x$cve[, x$min_ind[2]]
  cvse <- x$cvse[, x$min_ind[2]]
  
  L <- cve - cvse
  U <- cve + cvse
  
  ind <- is.finite(l) & seq_along(cve) <= length(l)
  
  ylim <- range(c(L[ind], U[ind]), na.rm = TRUE)
  
  plot(
    l[ind],
    cve[ind],
    type = "n",
    xlab = xlab,
    ylab = "CV Error",
    xlim = rev(range(l[ind])),
    ylim = ylim,
    las = 1,
    bty = "n"
  )
  
  abline(v = l[x$min_ind[1]], lty = 2, col = "black")
  
  arrows(
    l[ind],
    L[ind],
    l[ind],
    U[ind],
    code = 3,
    angle = 90,
    col = "gray80",
    length = 0.05
  )
  
  points(l[ind], cve[ind], col = col, pch = 19, cex = 0.6)
  
  # Read selected-group counts directly from the coefficient path to avoid
  # requiring new prediction data.
  
  if (selected) {
    
    beta <- x$fit$betas
    
    if (length(dim(beta)) == 3) {
      beta <- beta[, , x$min_ind[2]]
    }
    
    beta <- beta[-1, , drop = FALSE]
    support <- beta != 0

    if (!is.null(x$fit$group) && length(x$fit$group) == nrow(support)) {
      group <- x$fit$group
      n.s <- apply(support, 2, function(active) {
        length(unique(group[active]))
      })
    } else {
      n.s <- apply(support, 2, sum)
    }
    
    axis(
      3,
      at = l[ind],
      labels = n.s[ind],
      tick = FALSE,
      line = -0.5,
      cex.axis = 0.8
    )
    
    mtext("Groups selected", side = 3, line = 1.5)
  }
}
