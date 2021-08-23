## Florian Gerber, gerber@mines.edu, July 17, 2019
## -----------------------------------------------

## load packages and set defaults ----------------
rm(list=ls())

library(LatticeKrig)
library(fields)
library(optimParallel)
options(optimParallel.loginfo=TRUE)

## a plot function for optimParallel
plot.optimParallel <- function(oo, legendPosition="topleft", legendText=NULL, lwd=1.5, main="", ...){
    opar <- par("xpd", "mai")
    on.exit(par(opar))

    par(mai=c(.8,1,.8,1))
    indexPar <- 1:(ncol(oo$loginfo)/2-1)+1
    nPar <- length(indexPar)
    indexFn <- (ncol(oo$loginfo)/2)+1
    colPar <- colorRampPalette(c(c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", 
"#44AA99", "#999933", "#AA4499")))(length(indexPar))
    colFn <- "black"

    if(is.null(legendText))
        legendText <- colnames(oo$loginfo)[c(indexPar, indexFn)]
    else if(length(legendText) != nPar+1)
        stop("argument 'legendText' has to be of length #parameters + 1.")
    matplot(oo$loginfo[,indexPar], type="l", lty=1, yaxt="n", col=colPar,
            xlab="", ylab="pars", panel.first={grid(col="lightblue", lty=1, lwd=.3)}, lwd=lwd, main=main, ...)
    mtext("step", side=1, line=2)
    axis(2, las=2)
    par(new=TRUE)
    plot(oo$loginfo[,1], oo$loginfo[,indexFn], xlab="", xaxt="n", yaxt="n", ylab="", type="l", col=colFn, lty=2, lwd=lwd, ...)
    axis(4, las=2, col=colFn)
    corners  <-  par("usr") # Gets the four corners of plot area (x1, x2, y1, y2)
    par(xpd = TRUE) # Draw outside plot area
    text(x = corners[2]*1.15, y = mean(corners[3:4]), "fn()", srt = 270)
    par(xpd = FALSE)
    legend(legendPosition, legend=legendText, col=c(colPar, colFn), lty=c(rep(1, nPar), 2) , bg="white", lwd=lwd, ...)
}
## check for global variables
codetools::findGlobals(plot.optimParallel, merge=FALSE)$variables

set.seed(144)

formals(quilt.plot)$nx <- 150
formals(quilt.plot)$ny <- 150

LKrigCV <- function(xTrain=NULL, xVali=NULL,
                    zTrain=NULL, zVali=NULL,
                    yTrain=NULL, yVali=NULL,
                    beta=NULL,
                    a.wght=NULL,
                    lambda=NULL,
                    LKinfo=NULL,
                    ...,
                    cache=NULL,
                    verbose = FALSE,
                    optimOut=FALSE){

    ## setup LKinfo ---------------------------------------------------
    if(is.null(cache)){
        xTrain <- as.matrix(xTrain)
        if (any(duplicated(cat.matrix(xTrain)))) 
            warning("Not all xTrain locations are unique: see the results of\n                           duplicated(cat.matrix(xTrain)) ")
        if (any(is.na(xTrain))) {
            stop("Missing values in xTrain not allowed ")
        }
        
        if(is.null(LKinfo)){
            LKinfo <- do.call("LKrigSetup",
                              c(list(x = apply(rbind(xTrain, xVali), 2,"range")),
                                a.wght=a.wght, lambda=lambda, 
                                list(...), 
                                list(verbose = verbose)))
        } else {
            LKinfo <- do.call("LKinfoUpdate", c(list(LKinfo = LKinfo), 
                                                a.wght=a.wght, lambda=lambda, list(...)))
        }
        LKinfo$fixedFunction <- LKinfo$fixedFunctionArgs <- NULL
        
        if (is.na(LKinfo$lambda) || is.null(LKinfo$lambda) || length(LKinfo$lambda) > 1) 
            stop("Must specify a scalar lambda in call to LKrig or in LKinfo")
        
        if (LKinfo$dense) 
            stop("'LKinfo$dense=TRUE' not supported")
    } else {
        LKinfo <- cache$LKinfo
        LKarg <- c(list(LKinfo = LKinfo), list(...))
        if(!is.null(lambda))
            LKarg <- c(LKarg, lambda=lambda)
        if(!is.null(a.wght))
            LKarg <- c(LKarg, a.wght=a.wght)
        LKinfo <- do.call("LKinfoUpdate", LKarg)

        if(!is.null(xTrain) || !is.null(xVali) || !is.null(zTrain) || !is.null(zVali) || !is.null(yTrain) || !is.null(yVali))
            stop("if 'cache' is provided 'xTrain', 'xVali', 'zTrain', 'zVali', 'yTrain', and 'yVali' have to be NULL")
        
        ## we assume that LKinfo has not changed!!!
        ## if(!identical(LKinfo[!(names(LKinfo) %in% c("a.wght", "lambda"))],
        ##               cache$LKinfo[!(names(cache$LKinfo) %in% c("a.wght", "lambda"))]))
        ##     stop("LKinfo other than 'a.wght' or 'lambda' changed: cache cannot be used!")

        xTrain <- cache$xTrain
        xVali <- cache$xVali
        zTrain <- cache$zTrain
        zVali <- cache$zVali
        yTrain <- cache$yTrain
        yVali <- cache$yVali
    }

    ## decide what to update -------------------------------
    updateInfo <- c(basis=TRUE, Q=TRUE, G=TRUE, beta=TRUE)
    if(!is.null(cache)){
        updateInfo["basis"] <- FALSE
        if(identical(LKinfo$a.wght, cache$LKinfo$a.wght)){
            updateInfo["Q"] <- FALSE
            if(identical(LKinfo$lambda, cache$LKinfo$lambda))
                updateInfo["G"] <- FALSE
        }
        if(is.null(beta) || identical(beta, cache$beta)){
            beta <- cache$beta
            updateInfo["beta"] <- FALSE
        }
    }

    ## basis ----------------------    
    timeBasis <- system.time({
        if(updateInfo["basis"]){
            wX <- LKrig.basis(xTrain, LKinfo); 
            wXcrossprod <- crossprod(wX);
            PHIg <- LKrig.basis(xVali, LKinfo);
        } else {
            wX <- cache$wX;
            wXcrossprod <- cache$wXcrossprod;
            PHIg <- cache$PHIg;
        }
    })
    
    ## Q and G and Chol(G)-----------
    timeQ <- system.time({
        if(updateInfo["Q"]){
            Q <- LKrig.precision(LKinfo, verbose = verbose)
        } else {
            Q <- cache$Q
        }
    })
    timeG <- system.time({
        if(updateInfo["G"]){
            G <- wXcrossprod + LKinfo$lambda * Q
        } else {
            G <- cache$G
        }
    })

    timeChol <- system.time({
        if(updateInfo["G"]){
            if(is.null(cache)){
                GCholesky <- chol.spam(G, memory = LKinfo$choleskyMemory)
            } else {
                GCholesky <- update.spam.chol.NgPeyton(cache$GCholesky, G)
            } 
        } else {
            GCholesky <- cache$GCholesky
        }
    })

    ## combine with linear term -------
    timePred <- system.time({
        if(updateInfo["beta"]){
            y <- as.matrix(yTrain - zTrain %*% beta);
            if (!is.null(y) && any(is.na(y))) 
                stop("Missing values in y not allowed ");
        } else {
             y <- cache$y
        }
        c.coef <- solve.spam(GCholesky, crossprod(wX, y));
        pred <- c(PHIg %*% c.coef) + zVali %*% beta;
        SSE <- sum((yVali - pred)^2)
    })
    if(optimOut)
        return(SSE)
    time <- rbind(basis=timeBasis, Q=timeQ, G=timeG, Chol=timeChol, pred=timePred)[,1:3]
    time <- rbind(time, all=colSums(time))
    list(pred=pred, SSE=SSE, time=time, updateInfo=updateInfo,
         cache=list(LKinfo=LKinfo, wX=wX, wXcrossprod=wXcrossprod,
                    PHIg=PHIg, Q=Q, G=G, GCholesky=GCholesky, y=y, beta=beta,
                    xTrain=xTrain, xVali=xVali, zTrain=zTrain, zVali=zVali, yTrain=yTrain, yVali=yVali))
}

## check for global variables
codetools::findGlobals(LKrigCV, merge=FALSE)$variables


## simulate a spatial field ------------------------------------
n <- 30000
x <- cbind(runif(n), runif(n)) # spatial locations
beta <- c(1,2,3)               # linear predictor
z <- cbind(1, x)               # covariates
LKinfo <- LKrigSetup(cbind(c(0,1), c(0,1)),
                     NC=8, nlevel=3, alpha=(.5)^(0:2), lambda=.2, a.wght=15)
resid <- LKrig.sim(x1=x, LKinfo=LKinfo)       # spatial residuals without nugget
y <- z %*% beta + resid + rnorm(n, sd=.2)     # measured data
quilt.plot(x[,1], x[,2], y, main="simulated data")

## divide into training/validation data
index <- sample(c(TRUE,FALSE), size=n, replace=TRUE)   
xTrain <- x[index, ]
xValid <- x[!index, ]
yTrain <- y[index, ]
yValid <- y[!index, ]
zTrain <- z[index, ]
zValid <- z[!index, ]

## predict validation data with LatticeKrig, compute sum of squared errors SSE
## 1. use parameters from simulation
pp <- LKrigCV(x=xTrain, xVali=xValid, zTrain=zTrain, zVali=zValid,
              yTrain=yTrain, yVali=yValid, beta=c(1,2,3), LKinf=LKinfo)
## quilt.plot(xTrain[,1], xTrain[,2], z=yTrain, zlim=range(y, pp$pred))
## quilt.plot(xValid[,1], xValid[,2], z=pp$pred, zlim=range(y, pp$pred), add=TRUE)
## points(xValid, pch=".", cex=1)
pp$SSE
pp$updateInfo
pp$time

## 2. update SSE for a different 'beta'
pp2 <- LKrigCV(beta=c(1,2,3), a.wght=15, lambda=3, cache=pp$cache)
## quilt.plot(xTrain[,1], xTrain[,2], z=yTrain, zlim=range(y, pp2$pred))
## quilt.plot(xValid[,1], xValid[,2], z=pp2$pred, zlim=range(y, pp2$pred), add=TRUE)
## points(xValid, pch=".")
pp2$SSE
pp2$updateInfo # only a small part has to be computed...
pp2$time 

## update SSE for a different 'lambda'
pp3 <- LKrigCV(lambda=0.01, cache=pp$cache)
## quilt.plot(xTrain[,1], xTrain[,2], z=yTrain, zlim=range(y, pp3$pred))
## quilt.plot(xValid[,1], xValid[,2], z=pp3$pred, zlim=range(y, pp3$pred), add=TRUE)
## points(xValid, pch=".")
pp3$SSE
pp3$updateInfo # Q can be reused...
pp3$time       
1-pp3$time["Chol", 3]/pp$time["Chol", 3] # savings due to re-using symbolic factorization


## update SSE for a different 'a.wght'
pp4 <- LKrigCV(a.wght=13, cache=pp$cache)
quilt.plot(xTrain[,1], xTrain[,2], z=yTrain, zlim=range(yTrain, pp4$pred))
quilt.plot(xValid[,1], xValid[,2], z=pp4$pred, zlim=range(yTrain, pp4$pred), add=TRUE)
points(xValid, pch=".")
pp4$SSE
pp4$updateInfo
pp4$time       
pp4$time["Chol", 3]/pp$time["Chol", 3] # savings due to re-using symbolic factorization


## find 'beta' that lead to optimal smales SSE
cl <- makeForkCluster(detectCores(), outfile="")
setDefaultCluster(cl)
oo <- optimParallel(coef(lm(y ~ z-1)), 
                    fn=function(x) {
    print(unname(x), digits=3)
    LKrigCV(beta=x, cache=pp$cache, optimOut=TRUE)})

oo$par
lm(y ~ z-1) # quite different from start values ... 

plot.optimParallel(oo, legendText=c(expression(beta[1]), expression(beta[2]),
                                    expression(beta[3]), "fn"))# 

## find 'beta', 'a.wght', 'lambda' that lead to optimal (small) SSE
## takes 1.5 minutes...
oo2 <- optimParallel(c(coef(lm(y ~ z-1)), 4.5, .1), fn=function(x) {
    print(unname(x), digits=3)
    LKrigCV(beta=x[1:3], a.wght=x[4], lambda=x[5],
              cache=pp$cache, optimOut=TRUE)},
    lower=c(-Inf, -Inf, -Inf, 4.001, 0.001))

oo2$par
plot.optimParallel(oo2, legendText=c(expression(beta[1]), expression(beta[2]),
                                     expression(beta[3]), "a.wght", expression(lambda), "fn"))
stopCluster(cl)


