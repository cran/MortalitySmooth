Mort1Dsmooth.update <-
function(x, y, offset, wei, psi2, B,
                                lambdaP, a){
    # Input:
        # x: abcissae of data
        # y: count response
        # offset: an a priori known component
        # wei: weights
        # psi2: overdispersion parameter
        # B: B-splines basis
        # P: penalty term/matrix
        # a: coefficients
    
    # Output:
        # a: updated coefficients
    
    # linear predictor
    eta <- B%*%a
    # expected values
    mu <- exp(eta + offset)
    # weights
    w <- wei*mu
    W <- diag(c(w))
    # working response
    z <- wei*((y - mu)/mu + eta)
    # regression
    BtWB <- t(B) %*% W %*% B
    BtWz <- t(B) %*% (W %*% z)
    a <- solve(BtWB + psi2*lambdaP, BtWz)
    # coefficients
    a <- matrix(a, ncol = 1)
    a
}

