#' @export
plot.cv.sglasso <- function(x, log.l = TRUE, col = "red", selected = TRUE, ...) {
  
  l <- x$lambda
  
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
  }
  
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
  
  # ============================================================
  # GROUPS (NO predict!)
  # ============================================================
  
  if (selected) {
    
    beta <- x$fit$betas
    
    if (length(dim(beta)) == 3) {
      beta <- beta[, , x$min_ind[2]]
    }
    
    beta <- beta[-1, , drop = FALSE]
    
    n.s <- apply(beta != 0, 2, sum)
    
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