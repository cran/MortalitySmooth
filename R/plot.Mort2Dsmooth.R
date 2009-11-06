plot.Mort2Dsmooth <-
function(x, 
                              type=c("logrates", "deaths"), ...){
    type <- match.arg(type)
    xx <- x$x
    yy <- x$y
    ZZ <- x$Z*x$w
    E <- exp(x$offset)*x$w
    fitted.values <- x$fitted.values * x$w
    list. <- list(xx=xx, yy=yy, type=c("Actual", "Fitted"))
    grid. <- expand.grid(list.)
    Plot <- switch(type,
                   logrates = 1,
                   deaths = 2)
    if(Plot==1){
        grid.$ZZZ <- c(as.vector(log(ZZ/E)), as.vector(log(fitted.values/E)))
        my.breaks <- quantile(grid.$ZZZ, prob=seq(0,1,0.1))
        my.col  <- rainbow(length(my.breaks)-1)
    }
    if(Plot==2){
        grid.$ZZZ <- c(as.vector(ZZ), as.vector(fitted.values))
        my.breaks <- quantile(grid.$ZZZ, prob=seq(0,1,0.1))
        my.col  <- rainbow(length(my.breaks)-1)
    }
    levelplot(ZZZ ~ yy * xx | type, grid., layout=c(2,1), at=my.breaks, col.regions=my.col, 
              colorkey=list(col=my.col), ...)
}

