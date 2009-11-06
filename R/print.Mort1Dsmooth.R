print.Mort1Dsmooth <-
function(x, digits=max(3, getOption("digits")-3), ...){
    if(!is.null(cl <- x$call)){
        cat("Call:\n")
       dput(cl, control=NULL)
    }
    cat("\nNumber of Observations              :", x$n, "\n")
    cat("Effective dimension                 :", format(round(x$df, 2)), "\n")
    cat("(Selected) smoothing parameter      :", round(x$lambda,3), "\n")
    cat("Bayesian Information Criterion (BIC):", round(x$bic,3))
    cat("\n")
    invisible(x)
}

