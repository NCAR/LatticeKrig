print.LKinfo <- function(x, ...) {
    LKinfo <- x
    L <- LKinfo$nlevel
    cat("Number of levels:", L, fill = TRUE)
    cat(" ", fill = TRUE)
    cat("alpha:", unlist(x$alpha), fill = TRUE)
    if (!is.null(x$nu)) {
        cat("based on nu: ", x$nu, fill = TRUE)
    }
    cat("a.wght", unlist(x$a.wght), fill = TRUE)
    cat(" ", fill = TRUE)
    cat("grid sizes and number of basis functions", fill = TRUE)
    temp <- cbind(LKinfo$mx, LKinfo$my, diff(LKinfo$offset))
    dimnames(temp) <- list(paste("level", 1:LKinfo$nlevel), c("mx", 
        "my", "total"))
    print(temp)
    cat("total number of basis functions:", LKinfo$m, fill = TRUE)
    cat(" ", fill = TRUE)
    #
    if (LKinfo$NC.buffer > 0) {
        cat("A buffer of ", LKinfo$NC.buffer, " is included.", 
            fill = TRUE)
        cat("grid sizes and number of basis functions within spatial domain", 
            fill = TRUE)
        temp <- cbind(LKinfo$mx - 2 * LKinfo$NC.buffer.x, LKinfo$my - 
            2 * LKinfo$NC.buffer.y, LKinfo$m.domain)
        dimnames(temp) <- list(paste("level", 1:LKinfo$nlevel), 
            c("mx", "my", "total"))
        print(temp)
        cat("total number of basis functions inside domain: ", 
            sum(temp[, 3]), fill = TRUE)
    }
    #
    cat(" ", fill = TRUE)
    cat("grid info: ranges of spatial domain (excluding the buffer regions)", 
        fill = TRUE)
    temp <- unlist(LKinfo$grid.info)
    names(temp) <- names(LKinfo$grid.info)
    print(temp)
    cat(" ", fill = TRUE)
    #
    if (LKinfo$NC.buffer > 0) {
        cat("Spatial ranges of the grid including buffer regions", 
            fill = TRUE)
        temp <- matrix(NA, nrow = LKinfo$nlevel, ncol = 4)
        dimnames(temp) <- list(paste("level", 1:LKinfo$nlevel), 
            c("xmin", "xmax", "ymin", "ymax"))
        for (k in 1:LKinfo$nlevel) {
            grid.temp <- (LKinfo$grid[[k]])
            temp[k, ] <- c(range(grid.temp$x), range(grid.temp$y))
        }
        print(temp)
    }
    cat(" ", fill = TRUE)
    cat("Distance type: ", LKinfo$distance.type, fill = TRUE)
}



