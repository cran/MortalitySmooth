Mort1Dsmooth.checker <-
function(x, y, offset, w,
                                 overdispersion,
                                 ndx, deg, pord, 
                                 lambda, df, method, 
                                 control){
    # Input:
        # x: abcissae of data
        # y: count response
        # offset: an a priori known component (optional)
        # w: weights
        # overdispersion: logical on the presence of overdispersion
        # ndx: number of internal knots -1. Default: floor(length(x)/5)
        # deg: degree of the B-splines. Default: 3
        # pord: order of differences. Default: 2
        # lambda: smoothing parameter. Default: NULL (optional)
        # df: a number which specifies the degrees of freedom. Default: NULL (optional)
        # method: the method for controlling the amount of smoothing. Default: 4
        # control: a list of control parameters
            
    # Output: a list containing CHECKED arguments for the Mort1Dsmooth function
   
    # about x and y:
    if(missing(y)){
        if(is.list(x)){
            if(any(is.na(match(c("x", "y"), names(x)))))
               stop("cannot find x and y in list")
            x <- x$x
            y <- x$y
        }
        else if(is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
        }
        else if(is.matrix(x) && ncol(x) == 2) {
            y <- x[, 2]
            x <- x[, 1]
        }
        else {
            y <- x
            x <- time(x)
        }
    }
    # about the offset
    m <- length(y)
    if (missing(offset)) offset <- rep(0, m)
    if (length(offset) == 1) offset <- rep(offset, m)
    # about lengths and wrong values
    if (length(x)!=length(y)) 
        stop("Arguments must have same length")
    if (length(y) != m | length(offset) != m) 
        stop("Argument arrays of wrong length")
    if (deg <= 1 | deg >= 10) 
        stop("Wrong value for deg")
    if (pord <= 0 | pord >= 5) 
        stop("Wrong value for pord")
    if (ndx <= 3 | ndx >= floor(m*.8)) 
        stop("Wrong value for ndx")
    # about method
    if (method != 1 & method != 2 & method != 3 & method != 4) 
        stop("Wrong value for method")
    # method = 1 adjusts lambda so that the BIC is minimized.
    # method = 2 adjusts lambda so that the AIC is minimized.
    # method = 3 uses the value supplied for lambda. 
    # method = 4 adjusts lambda so that the degrees of freedom is equal to df.
    # check-point methods
    lambda.check <- is.null(lambda)
    df.check <- is.null(df)
    MET <- NULL
    if(lambda.check & df.check & method==3){MET=3}
    if(lambda.check & df.check & method==4){MET=4}
    if(!lambda.check & df.check & method==1){MET=1}
    if(lambda.check & !df.check & method==2){MET=2}
    # impossible values for lambda and df
    if(!lambda.check && lambda<0) stop("lambda must be positive")
    if(!df.check && df<pord) stop("df must be larger than pord")
    if(!df.check && df>ndx) stop("df must be smaller than ndx")
    # about the overdispersion parameter
    if(overdispersion & method==3) warning("given method 3, overdispersion is computed a posteriori")
    if(overdispersion & method==4) warning("given method 4, overdispersion is computed a posteriori")
    # error messages for methods
    if(is.null(MET)){
        if(lambda.check & df.check & method==3) stop("provide lambda")
        if(lambda.check & df.check & method==4) stop("provide df")
        if(!lambda.check & df.check & method==1) stop("lambda and method 1 cannot be chosen together")
        if(!lambda.check & df.check & method==2) stop("lambda and method 2 cannot be chosen together")
        if(!lambda.check & df.check & method==4) stop("lambda and method 4 cannot be chosen together")
        if(lambda.check & !df.check & method==1) stop("df and method 1 cannot be chosen together")
        if(lambda.check & !df.check & method==2) stop("df and method 2 cannot be chosen together")
        if(lambda.check & !df.check & method==3) stop("df and method 3 cannot be chosen together")
        if(!lambda.check & !df.check & method==1) stop("lambda, df and method 1 cannot be chosen together")
        if(!lambda.check & !df.check & method==2) stop("lambda, df and method 2 cannot be chosen together")
        if(!lambda.check & !df.check & method==3) stop("lambda, df and method 3 cannot be chosen together")
        if(!lambda.check & !df.check & method==4) stop("lambda, df and method 4 cannot be chosen together")
    }
    if (!df.check & length(df)!=1)
        stop("df must be length 1")
    if (!lambda.check & length(lambda)!=1)
        stop("lambda must be length 1")
    # setting control-parameters
    con <- list(MON=FALSE, TOL1=1e-06, TOL2=0.1, MAX.IT=50)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
    # about weights
    if(min(w) < 0) {
        warning(paste("At least one weight entry is negative"))
    }
    if(length(w)!=m){ 
        stop("length of w and y must be equal")
    }
    # returning
    llist <- list(x=x, y=y, offset=offset, w=w,
                  overdispersion=overdispersion, m=m,
                  ndx=ndx, deg=deg, pord=pord, 
                  lambda=lambda, df=df, method=method, 
                  control=con)
    llist
}

