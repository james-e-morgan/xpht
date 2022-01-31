#' @export
plot.extDiagram <- function(x, maxDirection,
                            directions = NULL,
                            barcode = FALSE,
                            rotated = FALSE,
                            col = NULL,
                            legendPos = NULL,
                            lab.line = 2.2, ...) {
    # Add checks for diagram list
    # Add checks for diagram Names
    
    if (is.null(directions)) {
        directions <- 1:maxDirection
    } else if (length(directions) == 1) {
        directions <- c(directions)
    }
 
    for (d in directions) {
        X = vector()
        
        imgName <- head(strsplit(x[[1]][[d]], "-")[[1]], -1)
        pName <- imgName[1]
        if (length(imgName) > 1) {
            for (i in 2:length(imgName)) {
                pName <- paste(pName, imgName[i], sep = " ")
            }
        }
        #print(pName)
        denom <- as.character(maxDirection/2)
        numer <- as.character(d)
        directionName <- bquote(.(pName) ~ (.(numer) * pi / .(denom)))
        
        for (i in 2:5) {

            if (length(x[[i]][[d]]) > 2) {
                if (i %% 2 == 1) {
                    sortedx <- x[[i]][[d]][order(x[[i]][[d]][,1], decreasing = TRUE),]
                } else {
                    sortedx <- x[[i]][[d]][order(x[[i]][[d]][,1], decreasing = FALSE),]
                }
                sortedx <- cbind(rep(i-1, length(sortedx)/2),
                                 sortedx)
                X <- rbind(X, sortedx)
            } else if (length(x[[i]][[d]]) == 2) {
                sortedx <- cbind(i-1, x[[i]][[d]])
                X <- rbind(X, sortedx)
            }
        }
        
        if (length(X) > 0) {
            plotLim <- c(min(X[,2:3]), max(X[,2:3]))
        } else {
            plotLim <- c(0,0)
        }
        
        if (is.null(col)) {
            col <- vector()
            for (i in 1:length(X[,1])) {
                if (X[i,1] == 1) {
                    col <- append(col, "#000000")
                } else if (X[i,1] == 2) {
                    col <- append(col, "#e69f00")
                } else if (X[i,1] == 3) {
                    col <- append(col, "#009e73")
                } else {
                    col <- append(col, "#d55e00")
                }
            }
        }
        
        if (barcode) {
            birth <- X[,2]
            death <- X[,3]
            
            first <- min(birth)
            last <- max(death)
            
            n <- length(birth)
            
            if (is.null(legendPos)) {
                legendPos <- "topright"
            }
            
            graphics::plot(c(first, last),
                           c(1,n+1),
                           type = "n",
                           xlab = "",
                           ylab = "",
                           xlim = c(first, last),
                           ylim = c(0, n+1),
                           xaxt = "n",
                           yaxt = "n", ...)
            
            graphics::axis(1)
            graphics::title(main = directionName,
                            xlab = "height",
                            line = lab.line)
            graphics::segments(birth, 1:n,
                               death, 1:n,
                               lwd = rep(2,n),
                               lty = rep(1,n),
                               col = col)
            graphics::legend(legendPos,
                             legend = c(expression('Ord'[0]),
                                        expression('Rel'[1]),
                                        expression('Ext'[0]),
                                        expression('Ext'[1])),
                             col = col,
                             lwd = 2,
                             lty = 1)
        }
    }
}