#' @export
plot.extDiagram <- function(x,
                            barcode = FALSE,
                            rotated = FALSE,
                            col = NULL,
                            showLegend = FALSE,
                            legendPos = NULL,
                            lab.line = 2.2, ...) {
    # Add checks for diagram list
    # Add checks for diagram Names
    X = vector()
    for (i in 2:5) {
        if (length(x[[i]]) > 2) {
            sortedx <- x[[i]][order(x[[i]][,1], decreasing = FALSE), ]
            sortedx <- cbind(rep(i-1, length(sortedx)/2), sortedx)
            X <- rbind(X, sortedx)
        } else if (length(x[[i]]) == 2) {
            sortedx <- cbind(i-1, x[[i]])
            X <- rbind(X, sortedx)
        }
    }
    
    if (is.null(col)) {
        col <- c(rep("#000000", length(which(X[,1] == 1))),
                 rep("#e69f00", length(which(X[,1] == 2))),
                 rep("#009e73", length(which(X[,1] == 3))),
                 rep("#d55e00", length(which(X[,1] == 4))))
    }
    
    if (barcode) {
        birth <- X[,2]
        death <- X[,3]
        
        first <- min(birth)
        last <- max(death)
        
        n <- length(birth)
        
        graphics::plot(c(first,last),
                       c(1,n+1),
                       type = "n",
                       xlab = "",
                       ylab = "",
                       xlim = c(first, last),
                       ylim = c(0, n+1),
                       xaxt = "n",
                       yaxt = "n", ...)
        
        graphics::axis(1)
        
        graphics::title(xlab = "height",
                        line = lab.line)
        
        graphics::segments(birth, 1:n,
                           death, 1:n,
                           lwd = rep(2,n),
                           lty = rep(1,n),
                           col = col)
        
    }
}