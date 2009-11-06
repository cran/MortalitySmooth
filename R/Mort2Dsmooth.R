Mort2Dsmooth <-
function(x, y, Z, offset, W, overdispersion=FALSE,
                         ndx=c(floor(length(x)/5), floor(length(y)/5)),
                         deg=c(3,3), 
                         pord=c(2,2),
                         lambdas=NULL, df=NULL, method=1, 
                         control=list()){
    # Input:
        # x: abcissa of data (commonly ages)
        # y: ordinate of data (commonly years)
        # Z: matrix of count response (commonly death counts)
        # offset: an a priori known component (optional) (commonly log-exposures)
        # W: a matrix of weights to be used in the fitting process (optional)
        # overdispersion: logical on the presence of an overdispersion parameter. Default:FALSE

        # ndx: a vector for the numbers of internal knots -1 for both axes. Default: floor(length(x)/5) and floor(length(y)/5)
        # deg: a vector for the degrees of the B-splines for the x-axis and y-axis. Default: c(3,3)
        # pord: a vector for the order of differences for both axes. Default: c(2,2)
        
        # lambdas: a vector of smoothing parameters for both axes (optional)
        # df: degree of freedom for both axes (optional)
        # method: the method for controlling the amount of smoothing. Default: 1
        
        # control: a list of control parameters
            # MON: logical on monitoring
            # TOL1: convergence of the IWLS algorithm. Default: 1e-06.
            # TOL2: difference between two adjacent smoothing parameter in the grid search, log-scale. Default: 0.5.
            # MAX.IT: the maximum number of iterations. Default: 50
            
    # Output: a Mort1Dsmooth object containing
        # aic: Akaike Information Criterion
        # bic: Bayesian Information Criterion
        # dev: Poisson-deviance
        # lev: diag of the hat-matrix
        # df: effective dimension
        # lambda: the selected (given) lambda
        # ndx: number of internal knots -1
        # deg: degree of the B-splines
        # pord: order of differences
        # x: abcissae of data
        # y: count responses
        # offset: an a priori known component
        # fitted.values: fitted counts
        # eta.hat: fitted linear predictor
        # coef: fitted (penalized) B-splines coefficients
        # psi2: (estimated) overdispersion parameter
        # W: (only for weighted fits) the specified weights
        # call: the matched call
        # n: number of observations
        # tolerance: convergence tolerance
        # residuals: the deviance residuals


    # checkings:
    if(missing(W)){
        W <- matrix(1, length(x), length(y))
    }
    check <- Mort2Dsmooth.checker(x=x, y=y, Z=Z, offset=offset, W=W,
                                  overdispersion=overdispersion,
                                  ndx=ndx, deg=deg, pord=pord,
                                  lambdas=lambdas, df=df, method=method, 
                                  control=control)
    x <- check$x
    y <- check$y
    Z <- check$Z
    m <- check$m
    n <- check$n
    offset <- check$offset
    wei <- check$W
    over <- check$overdispersion
    ndx <- check$ndx
    deg <- check$deg
    lambdas <- check$lambdas
    df <- check$df
    MET <- check$method
    MON <- check$control$MON
    TOL1 <- check$control$TOL1
    TOL2 <- check$control$TOL2
    MAX.IT <- check$control$MAX.IT
    call <- match.call()
    # B-splines basis for the abscissa (x)
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    Bx <- MortSmooth.bbase(x, xmin, xmax, ndx[1], deg[1])
    nbx <- ncol(Bx)
    # B-splines basis for the ordinate (y)
    yl <- min(y)
    yr <- max(y)
    ymax <- yr + 0.01 * (yr - yl)
    ymin <- yl - 0.01 * (yr - yl)
    By <- MortSmooth.bbase(y, ymin, ymax, ndx[2], deg[2])
    nby <- ncol(By)
    # Row tensors of B-splines basis
    Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), Bx)
    Bx2 <- kronecker(Bx, matrix(1, ncol=nbx,nrow=1))
    RTBx <- Bx1*Bx2
    By1 <- kronecker(matrix(1, ncol=nby, nrow=1), By)
    By2 <- kronecker(By, matrix(1, ncol=nby, nrow=1))
    RTBy <- By1*By2
    # penalty stuff
    Dx <- diff(diag(nbx), diff=pord[1])
    Dy <- diff(diag(nby), diff=pord[2])
    Px <- kronecker(diag(nby), t(Dx)%*%Dx)
    Py <- kronecker(t(Dy)%*%Dy, diag(nbx))
    # General initialize:
    Z[is.na(Z)] <- 0
    eta0 <- log((Z + 1)) - offset
    # simple ridge penalized regression
    BBx <- solve(t(Bx)%*%Bx + diag(nbx) * 1e-6, t(Bx))
    BBy <- solve(t(By)%*%By + diag(nby) * 1e-6, t(By))
    a.init <- MortSmooth.BcoefB(BBx, BBy, eta0)
    psi2 <- 1
    # optimize AIC or BIC
    if(MET==1|MET==2){
        by.lambda <- length(seq(-8,8,by=TOL2))
        Mort2Dsmooth.opt.ic <- function(X){
          FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                       Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                       lambdas=c(X[1], X[2]), Px=Px, Py=Py, a.init=a.init, 
                                       MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
          return(ifelse(MET==2, FIT$aic, FIT$bic))
        }
        # if overdisperion is true
        if(over){
          tol.over <- 10
          i.over <- 0
          while(tol.over > 1e-03 && i.over < 5){
            i.over <- i.over+1
            lambdas.hat <- cleversearch(fn=Mort2Dsmooth.opt.ic, lower=c(-8, -8), upper=c(8,8), ngrid=by.lambda, 
                                        logscale=TRUE, verbose=FALSE)[[1]]
            FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                         Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                         lambdas=lambdas.hat, Px=Px, Py=Py, a.init=a.init, 
                                         MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
            # recalculating overdispersion parameter
            psi2.old <- psi2
            psi2 <- FIT$dev / (m*n - FIT$df)
            tol.over <- abs(psi2 - psi2.old)/abs(psi2)
          }
        }else{# if psi2==1
          lambdas.hat <- cleversearch(fn=Mort2Dsmooth.opt.ic, lower=c(-8, -8), upper=c(8,8), ngrid=by.lambda, 
                                      logscale=TRUE, verbose=FALSE)[[1]]
          FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                       Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                       lambdas=lambdas.hat, Px=Px, Py=Py, a.init=a.init, 
                                       MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
          # a posteriori calculation of overdispersion parameter
          psi2 <- FIT$dev / (m*n - FIT$df)
        }
        if(log10(lambdas.hat)[1]==8 | log10(lambdas.hat)[1]==-8) {
          warning(paste("optimal lambda for x at the edge of the grid."))
        }
        if(log10(lambdas.hat)[2]==8 | log10(lambdas.hat)[2]==-8) {
          warning(paste("optimal lambda for y at the edge of the grid."))
        }
      }
    # given lambda
    if(MET==3){
      lambdas.hat <- lambdas
      FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                   Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                   lambdas=lambdas.hat, Px=Px, Py=Py, a.init=a.init, 
                                   MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
      # a posteriori calculation of overdispersion parameter
      psi2 <- FIT$dev / (m*n - FIT$df)
    }
    # optimize given df
    if(MET==4){
        by.lambda <- length(seq(-8,8,by=TOL2))
        Mort2Dsmooth.opt.df <- function(X){
            FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                         Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                         lambdas=c(X[1], X[2]), Px=Px, Py=Py, a.init=a.init, 
                                         MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
            return(abs(FIT$df - df))
        }
        lambdas.hat <- cleversearch(fn=Mort2Dsmooth.opt.df, lower=c(-8, -8), upper=c(8,8), ngrid=by.lambda,
                                    logscale=TRUE, verbose=FALSE)[[1]]
        if(log10(lambdas.hat)[1]==8 | log10(lambdas.hat)[1]==-8) {
            warning(paste("optimal lambda for x at the edge of the grid."))
        }
        if(log10(lambdas.hat)[2]==8 | log10(lambdas.hat)[2]==-8) {
            warning(paste("optimal lambda for y at the edge of the grid.."))
        }
        FIT <- Mort2Dsmooth.estimate(x=x, y=y, Z=Z, offset=offset, wei=wei, psi2=psi2,
                                     Bx=Bx, By=By, nbx=nbx, nby=nby, RTBx=RTBx, RTBy=RTBy, 
                                     lambdas=lambdas.hat, Px=Px, Py=Py, a.init=a.init, 
                                     MON=MON, TOL=TOL1, MAX.IT=MAX.IT)
        # a posteriori calculation of overdispersion parameter
        psi2 <- FIT$dev / (m*n - FIT$df)
    }
    aic <- FIT$aic
    bic <- FIT$bic
    df <- FIT$df
    dev <- FIT$dev
    coef <- FIT$a
    psi2 <- psi2
    h <- FIT$h
    tolerance <- FIT$tol
    eta.hat <- matrix(MortSmooth.BcoefB(Bx, By, coef),m,n,dimnames = list(x,y))
    fitted.values <- matrix(exp(offset + eta.hat),m,n,dimnames = list(x,y))
    lambdas.hat <- lambdas.hat
    res <- sign(Z - fitted.values) * sqrt(2 * (Z * log(ifelse(Z == 0, 1, Z/fitted.values)) - (Z - fitted.values)))
    res <- matrix(res,m,n,dimnames = list(x,y))
    wei <- matrix(wei,m,n,dimnames = list(x,y))
    # output
    object <- list(call=call, m=m, n=n, tolerance=tolerance, residuals=res,
                   # diagnostic outcomes
                   aic=aic, bic=bic, lev=h, df=df, dev=dev, psi2=psi2,
                   # smoothing specifications
                   lambdas=lambdas.hat, ndx=ndx, deg=deg, pord=pord,
                   # observed values
                   x=x, y=y, Z=Z, offset=offset, w=wei,
                   # fitted values
                   fitted.values=fitted.values, linear.predictors=eta.hat + offset, coefficients=coef
                   )
    class(object) <- "Mort2Dsmooth"
    object
}

