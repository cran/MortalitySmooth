predict.Mort2Dsmooth <-
function(object, newdata = NULL,
                                 type = c("link", "response"),
                                 se.fit = FALSE, ...){
    # Input:
        # object: a Mort2Dsmooth object
        # newdata: an optional list in which to look for x and/or y variables with which to predict. If omitted, the fitted values are used.
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
            # check newdata as a list
            if(!is.list(newdata))
                stop("newdata must be a data frame")
            # check newdata names and variables
            NEWdata <- list(x=object$x, y=object$y)
            namesNEW <- names(NEWdata)
            NEWdata[(namesnew <- names(newdata))] <- newdata
            if (length(noNms <- namesnew[!namesnew %in% namesNEW]) > 0) 
                stop("unknown names in newdata: ", paste(noNms, collapse = ", "))
            # check newdata ranges
            NEWx <- NEWdata$x
            NEWy <- NEWdata$y
            ran.NEWx <- range(NEWx)
            ran.NEWy <- range(NEWy)            
            ran.x <- range(object$x)
            ran.y <- range(object$y)
            if(ran.NEWx[1] < ran.x[1] | ran.NEWx[2] > ran.x[2] | ran.NEWy[1] < ran.y[1] | ran.NEWy[2] > ran.y[2])
                stop("newdata must be within the range of the original axes")
            # B-splines basis for the new abscissa (x)
            xl <- min(object$x)
            xr <- max(object$x)
            xmax <- xr + 0.01 * (xr - xl)
            xmin <- xl - 0.01 * (xr - xl)
            Bx <- MortSmooth.bbase(NEWx, xmin, xmax, object$ndx[1], object$deg[1])
            # B-splines basis for the new ordinate (y)
            yl <- min(object$y)
            yr <- max(object$y)
            ymax <- yr + 0.01 * (yr - yl)
            ymin <- yl - 0.01 * (yr - yl)
            By <- MortSmooth.bbase(NEWy, ymin, ymax, object$ndx[2], object$deg[2])
            # new linear predictor
            etanew <- matrix(MortSmooth.BcoefB(Bx, By, object$coef), length(NEWx), length(NEWy), dimnames=list(NEWx, NEWy))
            # switching according to "type"
            pred <- switch(type,
                           link = etanew,
                           response = exp(etanew))
        }
    }else{
    ## Standard errors
        # B-splines basis for the abscissa (x)
        xl <- min(object$x)
        xr <- max(object$x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        Bx <- MortSmooth.bbase(object$x, xmin, xmax, object$ndx[1], object$deg[1])
        nbx <- ncol(Bx)
        # B-splines basis for the ordinate (y)
        yl <- min(object$y)
        yr <- max(object$y)
        ymax <- yr + 0.01 * (yr - yl)
        ymin <- yl - 0.01 * (yr - yl)
        By <- MortSmooth.bbase(object$y, ymin, ymax, object$ndx[2], object$deg[2])
        nby <- ncol(By)
        # Row tensors of B-splines basis
        Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), Bx)
        Bx2 <- kronecker(Bx, matrix(1, ncol=nbx,nrow=1))
        RTBx <- Bx1*Bx2
        By1 <- kronecker(matrix(1, ncol=nby, nrow=1), By)
        By2 <- kronecker(By, matrix(1, ncol=nby, nrow=1))
        RTBy <- By1*By2
        # penalty stuff
        Dx <- diff(diag(nbx), diff=object$pord[1])
        Dy <- diff(diag(nby), diff=object$pord[2])
        Px <- kronecker(diag(nby), t(Dx)%*%Dx)
        Py <- kronecker(t(Dy)%*%Dy, diag(nbx))
        P <- (object$lambdas[1] * Px) + (object$lambdas[2] * Py)
        # 
        eta <- MortSmooth.BcoefB(Bx, By, object$coef)
        mu <- exp(object$offset + eta)
        W <- mu
        WW <- object$w*W
        BWB <- MortSmooth.BWB(RTBx, RTBy, nbx, nby, WW)
        # 
        BWB.P1 <- solve(BWB + P)
        if(missing(newdata)){
            se <- matrix(Mort2Dsmooth.se(RTBx=RTBx, RTBy=RTBy, nbx=nbx, nby=nby, BWB.P1=BWB.P1), object$m, object$n, dimnames=list(object$x, object$y))
            fit <- switch(type,
                          link = object$linear.predictors - object$offset,
                          response = object$fitted.values)
            se.fit <- switch(type,
                             link = se,
                             response = object$fitted.values * (exp(se) -1))
            pred <- list(fit=fit, se.fit=se.fit)
        }else{
            # check newdata as a list
            if(!is.list(newdata))
                stop("newdata must be a data frame")
            # check newdata names and variables
            NEWdata <- list(x=object$x, y=object$y)
            namesNEW <- names(NEWdata)
            NEWdata[(namesnew <- names(newdata))] <- newdata
            if (length(noNms <- namesnew[!namesnew %in% namesNEW]) > 0) 
                stop("unknown names in newdata: ", paste(noNms, collapse = ", "))
            # check newdata ranges
            NEWx <- NEWdata$x
            NEWy <- NEWdata$y
            ran.NEWx <- range(NEWx)
            ran.NEWy <- range(NEWy)            
            ran.x <- range(object$x)
            ran.y <- range(object$y)
            if(ran.NEWx[1] < ran.x[1] | ran.NEWx[2] > ran.x[2] | ran.NEWy[1] < ran.y[1] | ran.NEWy[2] > ran.y[2])
                stop("newdata must be within the range of the original axes")
                        # B-splines basis for the new abscissa (x)
            xl <- min(object$x)
            xr <- max(object$x)
            xmax <- xr + 0.01 * (xr - xl)
            xmin <- xl - 0.01 * (xr - xl)
            Bx <- MortSmooth.bbase(NEWx, xmin, xmax, object$ndx[1], object$deg[1])
            # B-splines basis for the new ordinate (y)
            yl <- min(object$y)
            yr <- max(object$y)
            ymax <- yr + 0.01 * (yr - yl)
            ymin <- yl - 0.01 * (yr - yl)
            By <- MortSmooth.bbase(NEWy, ymin, ymax, object$ndx[2], object$deg[2])
            # Row tensors of B-splines basis
            Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), Bx)
            Bx2 <- kronecker(Bx, matrix(1, ncol=nbx,nrow=1))
            RTBx <- Bx1*Bx2
            By1 <- kronecker(matrix(1, ncol=nby, nrow=1), By)
            By2 <- kronecker(By, matrix(1, ncol=nby, nrow=1))
            RTBy <- By1*By2
            # new linear predictor
            etanew <- matrix(MortSmooth.BcoefB(Bx, By, object$coef), length(NEWx), length(NEWy), dimnames=list(NEWx, NEWy))
            # standard errors ver the new data
            se <- matrix(Mort2Dsmooth.se(RTBx=RTBx, RTBy=RTBy, nbx=nbx, nby=nby, BWB.P1=BWB.P1), length(NEWx), length(NEWy), dimnames=list(NEWx, NEWy))
            # 
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

