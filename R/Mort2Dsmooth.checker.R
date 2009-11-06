Mort2Dsmooth.checker <-
function(x, y, Z, offset, W,
                                 overdispersion,
                                 ndx, deg, pord,
                                 lambdas, df, method, 
                                 control){
    # Input:
        # x: abcissa of data
        # y: ordinate of data
        # Z: matrix of count response
        # offset: an a priori known component (optional)
        # overdispersion: logical on the presence of overdispersion
        # W: a matrix of weights to be used in the fitting process 

        # ndx: a vector for the numbers of internal knots -1 for both axes. 
        # deg: a vector for the degrees of the B-splines for the x-axis and y-axis. 
        # pord: a vector for the order of differences for both axes. 
        
        # lambdas: a vector of smoothing parameters for both axes (optional)
        # df: degree of freedom for both axes (optional)
        # method: the method for controlling the amount of smoothing. 
        # control: a list of control parameters
        
    # Output: a list containing CHECKED arguments for the Mort2Dsmooth function
   
    # about x and y:
    if(missing(x)){
        x <- time(nrow(Z))
    }
    if(missing(y)){
        y <- time(ncol(Z))
    }
    # about the offset
    m <- length(x)
    n <- length(y)
    if (missing(offset)) offset <- matrix(0, m, n)
    if (length(offset) == 1) offset <- matrix(offset, m, n)
    # about lengths and wrong values
    if (length(x)!=nrow(Z)) 
        stop("length of x must be equal to number of rows in Z")
    if (length(y)!=ncol(Z)) 
        stop("length of y must be equal to number of columns in Z")
    if (dim(Z)[1] != m | dim(Z)[2] != n)
        stop("Argument arrays of wrong length")
    if (dim(offset)[1] != m | dim(offset)[2] != n)
        stop("Argument arrays of wrong length")
    if (deg[1] <= 1 | deg[1] >= 10 | deg[2] <= 1 | deg[2] >= 10) 
        stop("Wrong values for deg")    
    if (pord[1] <= 0 | pord[1] >= 5 | pord[2] <= 0 | pord[2] >= 5) 
        stop("Wrong value for pord")
    if (ndx[1] <= 3 | ndx[1] >= floor(m*.8) | ndx[2] <= 3 | ndx[2] >= floor(n*.8)) 
        stop("Wrong value for ndx")
    # about method
    if (method != 1 & method != 2 & method != 3 & method != 4) 
        stop("Wrong value for method")
    # method = 1 adjusts lambda so that the BIC is minimized.
    # method = 2 adjusts lambda so that the AIC is minimized.
    # method = 3 uses the value supplied for lambdas. 
    # method = 4 adjusts lambda so that the degrees of freedom is equal to df.
    # check-point methods
    lambdas.check <- is.null(lambdas)
    df.check <- is.null(df)
    MET <- NULL
    if(lambdas.check & df.check & method==1){MET=1}
    if(lambdas.check & df.check & method==2){MET=2}
    if(!lambdas.check & df.check & method==3){MET=3}
    if(lambdas.check & !df.check & method==4){MET=4}
    # impossible values for lambdas and df
    if(!lambdas.check && lambdas<0) stop("lambdas must be positive")
    if(!df.check && df < (pord[1] + pord[2])) stop("df must be larger than the sum of pord")
    if(!df.check && df>(ndx[1]*ndx[2])) stop("df must be smaller than the product of ndx values")
    # about the overdispersion parameter
    if(overdispersion & method==3) warning("given method 3, overdispersion is computed a posteriori")
    if(overdispersion & method==4) warning("given method 4, overdispersion is computed a posteriori")
    # error messages for methods
    if(is.null(MET)){
        if(lambdas.check & df.check & method==3) stop("provide lambdas")
        if(lambdas.check & df.check & method==4) stop("provide df")
        if(!lambdas.check & df.check & method==1) stop("lambdas and method 1 cannot be chosen together")
        if(!lambdas.check & df.check & method==2) stop("lambdas and method 2 cannot be chosen together")
        if(!lambdas.check & df.check & method==4) stop("lambdas and method 4 cannot be chosen together")
        if(lambdas.check & !df.check & method==1) stop("df and method 1 cannot be chosen together")
        if(lambdas.check & !df.check & method==2) stop("df and method 2 cannot be chosen together")
        if(lambdas.check & !df.check & method==3) stop("df and method 3 cannot be chosen together")
        if(!lambdas.check & !df.check & method==1) stop("lambdas, df and method 1 cannot be chosen together")
        if(!lambdas.check & !df.check & method==2) stop("lambdas, df and method 2 cannot be chosen together")
        if(!lambdas.check & !df.check & method==3) stop("lambdas, df and method 3 cannot be chosen together")
        if(!lambdas.check & !df.check & method==4) stop("lambdas, df and method 4 cannot be chosen together")
    }
    if (!df.check & length(df)!=1)
        stop("df must be length 1")
    if (!lambdas.check & length(lambdas)!=1 & length(lambdas)!=2)
        stop("lambda must be length 1 or 2")
    if (!lambdas.check & length(lambdas)==1){
        lambdas <- rep(lambdas, 2)
        warning("Isotropic smoothing is applied", call.=FALSE)
    }
    # setting control-parameters
    con <- list(MON=FALSE, TOL1=1e-06, TOL2=0.5, MAX.IT=50)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
    # about weights
    if(min(W) < 0) {
        warning(paste("At least one weight entry is negative"))
    }
    if(nrow(W)!=m | ncol(W)!=n){
        stop("dimensions of W and Z must be equal")
    }
    # returning
    llist <- list(x=x, y=y, Z=Z, offset=offset, W=W, m=m, n=n,
                  overdispersion=overdispersion,
                  ndx=ndx, deg=deg, pord=pord,
                  lambdas=lambdas, df=df, method=method, 
                  control=con)
    llist
}

