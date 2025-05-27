#' @export
plot.cv.sglasso <- function(x,log.l=T, type.tun=c("lambda","d"),
                            type=c("cve"),selected=TRUE,vertical.line=TRUE,col="red"){
  type     = match.arg(type)
  type.tun = match.arg(type.tun, c("lambda","d"))
  if(type.tun == "lambda"){
    l <- x$lambda
    if (log.l) {
      l <- log(l)
      xlab <- expression(log(lambda))
    } else xlab <- expression(lambda)
    ## Calculate y
    
    x$cve <- x$cve[,x$min_ind[2]]
    x$cvse <- x$cvse[,x$min_ind[2]]
    L.cve <- x$cve - x$cvse
    U.cve <- x$cve + x$cvse
    y <- x$cve
    L <- L.cve
    U <- U.cve
    ylab <- "Cross-validation error"
    #ind <- if (type=="pred") is.finite(l[1:length(x$pe)]) else is.finite(l[1:length(x$cve)])
    ind <- is.finite(l[1:length(x$cve)])
    ylim <- range(c(L[ind], U[ind]))
    aind <- ((U-L)/diff(ylim) > 1e-3) & ind
    plot.args = list(x=l[ind], y=y[ind], ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l[ind])), las=1, bty="n")
    #new.args = list(...)
    #if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    if (vertical.line) abline(v=l[x$min_ind[1]], lty=2, lwd=.5)
    suppressWarnings(arrows(x0=l[aind], x1=l[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
    points(l[ind], y[ind], col=col, pch=19, cex=.5)
    if (selected) {
      pred_ns <- predict(x$fit, lambda=x$lambda, d=x$d.min, type="groups")
      if(class(pred_ns)[1] != "list"){
        n.s <- rep(dim(pred_ns)[1],length(x$lambda))
      }else{
        n.s <- sapply(pred_ns, length)  
      }
      axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
      mtext("Groups selected", cex=0.8, line=1.5)
    }
  }
  else {
    d.val <- x$d
    xlab <- expression(d)
    x$cve <- x$cve[x$min_ind[1],]
    x$cvse <- x$cvse[x$min_ind[1],]
    L.cve <- x$cve - x$cvse
    U.cve <- x$cve + x$cvse
    y <- x$cve
    L <- L.cve
    U <- U.cve
    ylab <- "Cross-validation error"
    ind <- is.finite(d.val[1:length(x$cve)])
    ylim <- range(c(L[ind], U[ind]))
    aind <- ((U-L)/diff(ylim) > 1e-3) & ind
    plot.args = list(x=d.val[ind], y=y[ind], ylim=ylim, xlab=xlab, ylab=ylab, type="l", xlim=range(d.val[ind]), las=1, bty="n")
    #new.args = list(...)
    #if (length(new.args)) plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    if (vertical.line) abline(v=d.val[x$min_ind[2]], lty=2, lwd=.5)
    suppressWarnings(arrows(x0=d.val[aind], x1=d.val[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
    points(d.val[ind], y[ind], col=col, pch=19, cex=.5)
    if (selected) {
      pred_ns <- predict(x$fit, lambda=x$lambda.min, d=x$d, type="groups")
      if(class(pred_ns)[1] != "list"){
        n.s <- rep(dim(pred_ns)[1],length(x$d))
      }else{
        n.s <- sapply(pred_ns, length)  
      }
      axis(3, at=d.val, labels=n.s, tick=FALSE, line=-0.5)
      mtext("Groups selected", cex=0.8, line=1.5)
    }
  }
}
