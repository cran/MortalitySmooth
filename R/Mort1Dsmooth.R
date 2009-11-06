Mort1Dsmooth <-
function(x, y, offset, w, overdispersion=FALSE,
                         ndx=floor(length(x)/5), deg=3, pord=2, 
                         lambda=NULL, df=NULL, method=1, 
                         control=list()){
    # Input:
        # x: abcissae of data
        # y: count response
        # offset: an a priori known component (optional)
        # w: a vector of weights to be used in the fitting process (optional)
        # overdispersion: logical on the presence of an overdispersion parameter. Default: FALSE
        # ndx: number of internal knots -1. Default: floor(length(x)/5)
        # deg: degree of the B-splines. Default: 3
        # pord: order of differences. Default: 2
        # lambda: smoothing parameter (optional)
        # df: a number which specifies the degrees of freedom (optional)
        # method: the method for controlling the amount of smoothing. Default: 1
        # control: a list of control parameters
            # MON: logical on monitoring
            # TOL1: convergence of the IWLS algorithm. Default: 1e-06.
            # TOL2: difference between two adjacent smoothing parameters in the grid search, log-scale. Default: 0.1.
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
        # w: (only for weighted fits) the specified weights
        # call: the matched call
        # n: number of observations
        # tolerance: convergence tolerance
        # residuals: the deviance residuals
     
    # checkings:
    if(missing(w)){
        w <- rep(1, length(y))
    }
    check <- Mort1Dsmooth.checker(x=x, y=y,
                                  offset=offset, w=w,
                                  overdispersion=overdispersion,
                                  ndx=ndx, deg=deg, pord=pord, 
                                  lambda=lambda, df=df,
                                  method=method, 
                                  control=control)
    x <- check$x
    y <- check$y
    m <- check$m
    offset <- check$offset
    wei <- check$w
    over <- check$overdispersion
    ndx <- check$ndx
    deg <- check$deg
    pord <- check$pord
    lambda <- check$lambda
    df <- check$df
    MET <- check$method
    MON <- check$control$MON
    TOL1 <- check$control$TOL1
    TOL2 <- check$control$TOL2
    MAX.IT <- check$control$MAX.IT
    call <- match.call()
    # B-splines basis
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    B <- MortSmooth.bbase(x, xmin, xmax, ndx, deg)
    # penalty stuff
    nb <- ncol(B)
    D. <- diff(diag(nb), diff=pord)
    DtD <- t(D.)%*%D.
    # General initialize:
    y[is.na(y)] <- 0
    eta0 <- log(y+1)
    mu0 <- exp(eta0 + offset)
    w0 <- wei*mu0
    z0 <- wei*((y - mu0)/mu0 + eta0)
    # simple penalized-poisson-GLM
    P <- 1e+08 * DtD
    W0 <- diag(w0)
    BtWB <- t(B) %*% W0 %*% B
    BtWz <- t(B) %*% (W0 %*% z0)
    a.init <- solve(BtWB + P, BtWz)
    psi2 <- 1
    # optimize AIC or BIC
    if(MET==1|MET==2){    
        by.lambda <- length(seq(-8,8,by=TOL2))
        Mort1Dsmooth.opt.ic <- function(X){
              FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                               wei=wei,
                               psi2=psi2,
                               B=B,lambda=X, DtD=DtD,
                               a.init=a.init,
                               MON=MON, TOL=TOL1,
                               MAX.IT=MAX.IT)
              return(ifelse(MET==2, FIT$aic, FIT$bic))
            }
        # if overdisperion is true
        if(over){
          tol.over <- 10
          i.over <- 0
          while(tol.over > 1e-03 && i.over < 5){
            i.over <- i.over+1
            lambda.hat <- cleversearch(fn=Mort1Dsmooth.opt.ic,
                                       lower=-8, upper=8,
                                       ngrid=by.lambda,
                                       logscale=TRUE,
                                       verbose=FALSE)[[1]]
            FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                                         wei=wei,
                                         psi2=psi2,
                                         B=B, 
                                         lambda=lambda.hat,
                                         DtD=DtD,
                                         a.init=a.init, 
                                         MON=MON, TOL=TOL1,
                                         MAX.IT=MAX.IT)
            # recalculating overdispersion parameter
            psi2.old <- psi2
            psi2 <- FIT$dev / (m - FIT$df)
            tol.over <- abs(psi2 - psi2.old)/abs(psi2)
          }
        }else{# if psi2==1
          lambda.hat <- cleversearch(fn=Mort1Dsmooth.opt.ic,
                                     lower=-8, upper=8,
                                     ngrid=by.lambda,
                                     logscale=TRUE,
                                     verbose=FALSE)[[1]]
          FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                                       wei=wei,
                                       psi2=psi2,
                                       B=B, 
                                       lambda=lambda.hat,
                                       DtD=DtD,
                                       a.init=a.init, 
                                       MON=MON, TOL=TOL1,
                                       MAX.IT=MAX.IT)
          psi2 <- FIT$dev / (m - FIT$df)
        }
    if(log10(lambda.hat)==8 | log10(lambda.hat)==-8) {
            warning(paste("optimal lambda at the edge of the grid."))
          }
      }
    # given lambda
    if(MET==3){
        lambda.hat <- lambda
        FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                                     wei=wei,
                                     psi2=psi2,
                                     B=B, 
                                     lambda=lambda.hat, DtD=DtD,
                                     a.init=a.init, 
                                     MON=MON, TOL=TOL1,
                                     MAX.IT=MAX.IT)
        psi2 <- FIT$dev / (m - FIT$df)
      }
    # optimize given df
    if(MET==4){
      Mort1Dsmooth.opt.df <- function(X){
        FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                                     wei=wei,
                                     psi2=psi2,
                                     B=B, 
                                     lambda=X, DtD=DtD,
                                     a.init=a.init, 
                                     MON=MON, TOL=TOL1,
                                     MAX.IT=MAX.IT)
        return(abs(FIT$df - df))
      }
      by.lambda <- length(seq(-8,8,by=TOL2))
      lambda.hat <- cleversearch(fn=Mort1Dsmooth.opt.df,
                                 lower=-8,
                                 upper=8,
                                 ngrid=by.lambda,
                                 logscale=TRUE,
                                 verbose=FALSE)[[1]]
      if(log10(lambda.hat)==8 | log10(lambda.hat)==-8){
        warning(paste("optimal lambda at the edge of the grid."))
      }
      FIT <- Mort1Dsmooth.estimate(x=x, y=y, offset=offset,
                                   wei=wei,
                                   psi2=psi2,
                                   B=B, 
                                   lambda=lambda.hat,
                                   DtD=DtD,
                                   a.init=a.init, 
                                   MON=MON, TOL=TOL1,
                                   MAX.IT=MAX.IT)
      psi2 <- FIT$dev / (m - FIT$df)
    }
    aic <- FIT$aic
    bic <- FIT$bic
    df <- FIT$df
    dev <- FIT$dev
    coef <- FIT$a
    psi2 <- psi2
    h <- FIT$h
    tolerance <- FIT$tol
    eta.hat <- B%*%coef
    fitted.values <- exp(eta.hat + offset)
    res <- sign(y - fitted.values) * sqrt(2 * (y * log(ifelse(y == 0, 1, y/fitted.values)) - (y - fitted.values)))
    # output
    object <- list(# fitted values
                   coefficients=as.vector(coef),
                   residuals=as.vector(res),
                   fitted.values=as.vector(fitted.values),
                   linear.predictors=as.vector(eta.hat) + as.vector(offset),
                   # diagnostics
                   lev=h, df=df, deviance=dev, aic=aic, bic=bic,
                   psi2=psi2, lambda=lambda.hat,
                   # general
                   call=call, n=m, tolerance=tolerance,  
                   # smoothing specifications
                   ndx=ndx, deg=deg, pord=pord,
                   # observed values
                   x=x, y=as.vector(y),
                   offset=as.vector(offset), w=as.vector(wei)
                   )
    class(object) <- "Mort1Dsmooth"
    object
}

