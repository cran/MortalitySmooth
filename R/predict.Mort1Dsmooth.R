predict.Mort1Dsmooth <-
function(object, newdata = NULL,
                                 type = c("link", "response"),
                                 se.fit = FALSE, ...){
    # Input:
        # object: a Mort1Dsmooth object
        # newdata: optionally, a vector in which to look for x with
        # which to predict. If omitted, the fitted linear
        # predictors are used  
        # type: the type of prediction required. The default ("link") is on the scale of the linear predictors;
        #       the alternative "response" is on the scale of the response variable.
        # se.fit: logical switch indicating if standard errors are required

    # Output: a list with components 
        # fit: Predictions
        # se.fit: Estimated standard errors
    type <- match.arg(type)
    
    if (!se.fit) {
    ## No standard errors
        if(missing(newdata)) {
            pred <- switch(type,
                           link = object$linear.predictors - object$offset,
                           response = object$fitted.values)
        }else{
            ran.newdata <- range(newdata) 
            ran.x <- range(object$x)
            if(ran.newdata[1] < ran.x[1] | ran.newdata[2] > ran.x[2])
                stop("newdata must be within the range of x")
            xl <- min(object$x)
            xr <- max(object$x)
            xmax <- xr + 0.01 * (xr - xl)
            xmin <- xl - 0.01 * (xr - xl)
            Bnew <- MortSmooth.bbase(newdata, xmin, xmax, ndx=object$ndx, deg=object$deg)
            etanew <- as.vector(Bnew %*% as.vector(object$coef))
            pred <- switch(type,
                           link = etanew,
                           response = exp(etanew))
        }
    }else{
    ## Standard errors
        xl <- min(object$x)
        xr <- max(object$x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        B <- MortSmooth.bbase(object$x, xmin, xmax, ndx=object$ndx, deg=object$deg)
        nb <- ncol(B)
        BtWB <- t(B) %*% (as.vector(object$fitted.values) * B)
        D. <- diff(diag(nb), diff=object$pord)
        PtP <- object$lambda * t(D.) %*% D.
        solBtWB <- solve(BtWB + PtP)
        if(missing(newdata)){
            se <- sqrt(diag(B %*% solBtWB %*% BtWB %*% solBtWB %*% t(B)))
            fit <- switch(type,
                          link = object$linear.predictors - object$offset,
                          response = object$fitted.values)
            se.fit <- switch(type,
                             link = se,
                             response = object$fitted.values * (exp(se) -1))
            pred <- list(fit=fit, se.fit=se.fit)
        }else{
            ran.newdata <- range(newdata) 
            ran.x <- range(object$x)
            if(ran.newdata[1] < ran.x[1] | ran.newdata[2] > ran.x[2])
                stop("newdata must be within the range of x")
            Bnew <- MortSmooth.bbase(newdata, xmin, xmax, ndx=object$ndx, deg=object$deg)
            etanew <- as.vector(Bnew %*% as.vector(object$coef))
            se <- sqrt(diag(Bnew %*% solBtWB %*% BtWB %*% solBtWB %*% t(Bnew)))
            se.fit <- switch(type,
                             link = se,
                             response = (exp(se)-1) * exp(etanew))
            fit <- switch(type,
                          link = etanew,
                          response = exp(etanew))
            pred <- list(fit=fit, se.fit=se.fit)
            }
        }
    pred
}

