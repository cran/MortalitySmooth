print.summary.Mort1Dsmooth <-
function(x, digits=max(3, getOption("digits")-3), ...){
    if(!is.null(cl <- x$call)){
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    cat("\nNumber of Observations                  :", x$n, "\n")
    cat("Effective dimension                     :", format(signif(x$df, 3)), "\n")
    cat("(Selected) smoothing parameter          :", format(signif(x$lambda, 3)), "\n")
    cat("Bayesian Information Criterion (BIC)    :", signif(x$bic,3), "\n")
    cat("Akaike's Information Criterion (AIC)    :", signif(x$aic,3), "\n")
    cat("(Estimated) dispersion parameter (psi^2):", signif(x$psi2,3), "\n")
    cat("\nResiduals:\n", sep="")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(x$residuals), names = nam)
    print(rq, digits = digits, ...)
    cat("\nSettings and control:\n")
    cat("  number of B-splines    :", x$ndx + x$deg, "\n")
    cat("  degree of the B-splines:", x$deg, "\n")
    cat("  order of differences   :", x$pord, "\n")
    cat("  convergence tolerance  :", x$tolerance)
    cat("\n")
    invisible(x)
}

