plot.HMDdata <-
function(x, ...){
    at <- attributes(x)
    type <- strsplit(at$`country-data-sex`, split="-")[[1]][2]
    if(is.null(dim(x))){
        # one-dimensional graph
        xx <- as.numeric(at$names)
        if(type=="Rates"){
            logrates <- log(as.vector(x))
            plot(xx, logrates, ...)
        }
        if(type=="Population"){
            population <- as.vector(x)
            plot(xx, population, ...)
        }
        if(type=="Deaths"){
            deaths <- as.vector(x)
            plot(xx, deaths, ...)
        }
        if(type=="Exposures"){
            exposures <- as.vector(x)
            plot(xx, exposures, ...)
        }
    }else{
        x1 <- as.numeric(at$dimnames[[1]])
        x2 <- as.numeric(at$dimnames[[2]])
        listMSF <- list(A=x1, Y=x2)
        gridMSF <- expand.grid(listMSF)
        if(type=="Rates"){
            gridMSF$logrates <- log(as.vector(x))
        }
        if(type=="Population"){
            gridMSF$population <- as.vector(x)
        }
        if(type=="Deaths"){
            gridMSF$deaths <- as.vector(x)
        }
        if(type=="Exposures"){
            gridMSF$exposures <- as.vector(x)
        }
    levelplot(gridMSF[,3] ~ Y * A, gridMSF, ...)
    }
}

