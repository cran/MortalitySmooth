plot.Mort1Dsmooth <-
function(x, 
                              type=c("logrates", "deaths"),
                              ...){
    type <- match.arg(type)
    xx <- x$x
    yy <- x$y
    fitted.values <- x$fitted.values*x$w
    Plot <- switch(type,
                   deaths = 1,
                   logrates = 2)
    if(Plot==1){
        ran <- range(yy, fitted.values)
        plot(xx, yy, ylim=c(ran[1], ran[2]), ...)
        lines(xx, fitted.values, col=2, lwd=2)
    }
    if(Plot==2){
        whi <- which(x$y!=0)
        offset <- x$offset
        lograte <- log(yy) - offset
        ran <- range(lograte[whi], log(x$fitted) - offset)
        plot(xx, lograte, ylim=c(ran[1], ran[2]), ...)
        lines(xx, log(x$fitted) - offset, col=2, lwd=2)
    }
}

