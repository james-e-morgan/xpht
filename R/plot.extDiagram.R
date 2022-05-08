#' @export
plot.extDiagram <- function(x,
                            barcode = FALSE,
                            rotated = FALSE,
                            histogram = FALSE,
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
        
    } else {
        
        if (rotated == TRUE) {
            extrema <- c(X[,2] + X[,3], X[,3] - X[,2])
            
            first <- min(extrema)
            last <- max(extrema)
            
            symbs <- 15:18
            col <- unique(col)
            
            graphics::plot(0,0,
                           type = "n",
                           xlab = "",
                           ylab = "",
                           xlim = c(first, last),
                           ylim = c(first, last),
                           ...)
            
            for (i in 1:4) {
                subDgm <- X[X[,1] == i,]
                birth <- (subDgm[,2] + subDgm[,3])
                death <- (subDgm[,3] - subDgm[,2])
                
                graphics::points(birth,
                                death,
                                pch = symbs[i],
                                lwd = 2,
                                cex = 1,
                                col = col[i])
            }
            
            graphics::title(xlab = "Birth + Death",
                            ylab = "Death - Birth")
            
            if (showLegend) {
                graphics::legend(legendPos,
                                legend = c("Ord0", "Rel1", "Ext0", "Ext1"),
                                pch = symbs,
                                col = c("#000000","#e69f00","#009e73","#d55e00"))
            }
            
        } else {
            first <- min(X[,2])
            last <- max(X[,3])
            symbs <- 15:18
            col <- unique(col)
            #initialise empty plot
            graphics::plot(0, 0,
                        type = "n",
                        xlab = "",
                        ylab = "",
                        xlim = c(first, last),
                        ylim = c(first, last),
                        ...)
            
            for (i in 1:4) {
                subDgm <- X[X[,1] == i,]
                birth <- subDgm[,2]
                death <- subDgm[,3]
                
                graphics::points(birth,
                                death,
                                pch = symbs[i],
                                lwd = 2,
                                cex = 1,
                                col = col[i])
            }
            
            graphics::abline(0,1,
                            lty = 2,
                            col = "#A9A9A9")
            
            graphics::title(xlab = "birth",
                            ylab = "death")
            
            if (showLegend) {
                graphics::legend(legendPos,
                                legend = c("Ord0", "Rel1", "Ext0", "Ext1"),
                                pch = symbs,
                                col = c("#000000","#e69f00","#009e73","#d55e00"))
            } 
        }
    }
}