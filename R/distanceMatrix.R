#' @export
distanceMatrix <- function(diagrams,
                           nObjects,
                           q = 1,
                           aligned = FALSE) {
    # The extended diagram for an object in n directions
    nDirs <- length(diagrams)/nObjects
    
    distMatrix <- matrix(0, nrow = nObjects, ncol = nObjects)
    
    # Distance matrix is symmetric, so only compute upper triangular part.
    # Diagonal entries are zero.
    for (i in 1:(nObjects-1)) {
        idx1 <- ((i-1)*nDirs + 1) : (i*nDirs)
        for (j in (i+1):(nObjects)){
            idx2 <- ((j-1)*nDirs + 1) : (j*nDirs)
            
            obj1 <- lapply(idx1, function(k) diagrams[[k]])
            obj2 <- lapply(idx2, function(k) diagrams[[k]])
            
            if (aligned) {
                d_ij <- alignedDistance(object1 = obj1,
                                        object2 = obj2,
                                        q = q,
                                        K = nDirs)
            } else {
                d_ij <- unalignedDistance(object1 = obj1,
                                          object2 = obj2,
                                          q = q,
                                          K = nDirs)
            }
            
            d_ij <- (d_ij/nDirs)^(1/q)
            
            distMatrix[i,j] <- d_ij
            distMatrix[j,i] <- d_ij
            
            print(paste("Computed distance (", i, ",", j, ") and (",
                        j, ",", i, ")", sep=""))
        }
    }
    
    return(distMatrix)
}

unalignedDistance <- function(object1,
                              object2,
                              q,
                              K) {
    ds <- c("Ext0", "Ext1", "Ord0", "Rel1")
    # Compute aligned distance to give a starting point
    minDist <- alignedDistance(object1 = object1,
                               object2 = object2,
                               q = q,
                               K = K)
    if (minDist == 0) {
        return(minDist)
    }
    # Only compute distances until min_dist is exceeded
    for (m in 0:(K-1)) { # m controls offset
        total <- 0
        for (x in ds) {
            d <- sum(sapply(1:K,
                            function(k) pointDistance(X = object1[[k]][[x]],
                                                      Y = object2[[((k+m)%%K)+1]][[x]],
                                                      q = q)))

            total <- total + d
            if (total >= minDist) {
                break
            }                                                                       
        }
        if (total < minDist) {
            minDist <- total
        }
    }
    return(minDist)
}

alignedDistance <- function(object1,
                            object2,
                            q,
                            K) {
    ds <- c("Ord0", "Rel1", "Ext0", "Ext1")
    total <- 0
    
    for (x in ds) {
        d <- sum(sapply(1:K,
                        function(k) pointDistance(X = object1[[k]][[x]],
                                                  Y = object2[[k]][[x]],
                                                  q = q)))
        
        total <- total + d
    }
    
    return(total)
}

pointDistance <- function(X, Y, q) {
    nX <- length(X)/2
    nY <- length(Y)/2
    
    if (nX + nY == 0) {
        # No points
        return(0)
    } else if (nX > 0 && nY == 0) {
        # Only X has points
        d <- sum(sapply(1:nX,
                        function(i) diagonalDist(X[i,],q)))
       
        return(d)
    } else if (nX == 0 && nY > 0) {
        # Only Y has points
        d <- sum(sapply(1:nY,
                        function(i) diagonalDist(Y[i,],q)))
        
        return(d)
    } else {
        # Hungarian Algorithm
        costMatrix <- vector()
        
        for (i in 1:nX) {
            r1 <- apply(Y, 1, function(x) distanceLp(X[i,],x,q))
            
            distxy <- diagonalDist(X[i,],q)
            r2 <- rep(distxy,nX)
            
            costMatrix <- rbind(costMatrix, c(r1,r2))
        }
        
        r1 <- apply(Y, 1, function(x) diagonalDist(x,q))
        r2 <- rep(0,nX)
        
        for (j in 1:nY) {
            costMatrix <- rbind(costMatrix, c(r1,r2))
        }
        
        pairing <- RcppHungarian::HungarianSolver(costMatrix)
        idxs <- pairing[["pairs"]]
        
        pairVals <- apply(idxs, 1,
                          function(x) costMatrix[x[1],x[2]])
        
        return(sum(pairVals))
    }
}

diagonalDist <- function(p, q) {
    d <- 2 * ((abs(p[2]-p[1])/2)^q)
    return(d)
}

distanceLp <- function(p1, p2, q) {
    Lp <- abs(p1[1] - p2[1])^q + abs(p1[2] - p2[2])^q
}