extendedPersistence <- function(bdryCurves,
                                imgName,
                                nDirections,
                                tolerance = 1/sqrt(2),
                                saveOutput = FALSE,
                                outputDir = NULL,
                                fName = NULL,
                                verbose = TRUE) {
    if (nDirections < 2 || nDirections %% 2 == 1) {
        stop("Number of directions must be non-zero and even")
    } else {
        directions <- t(sapply(1:n, function(r) c(cos(2 * pi * r / n),
                                                  sin(2 * pi * r / n))))
    }
    
    xDiagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
    xDiagram <- sapply(xDiagram.names, function(x) NULL)
    
    oneSkeleton <- parseSkeleton(bdryCurves)
    
    acCurves <- bdryCurves[[1]][1]
    cCurves <- bdryCurves[[1]][2]
    
    for (d in 1:(directions/2)) {
        diagramName <- paste(imgName, "-", tostring(d), sep="")
        xDiagram[["Name"]][[d]] <- diagramName
        
        dirVector <- directions[d,]
        
        nComponents <- length(oneSkeleton)
        zerothDiagram <- vector(mode = "list", length = nComponents)
        
        count <- 1
        
        for (i in 1:nComponents) {
            skeletonComponent <- one_skeleton[i][[1]]
            if (i <= acCurves) {
                orientation = 0
            } else {
                orientation = 1
            }
            heightFiltration <- computeHeightFiltration(curve = skeletonComponent,
                                                        direction = dirVector,
                                                        orientation = orientation)
        }
    }
}

parseSkeleton <- function(bdryCurves) {
    
    complex <- vector(mode = "list")
    
    skeleton.names <- c("vertex", "edge", "coords")
    
    acCurves <- bdryCurves[[1]][1]
    cCurves <- bdryCurves[[2]][2]
    nCurves <- acCurves + cCurves
    
    for (i in 2:(nCurves + 1)) {
        curve <- bdryCurves[[i]]
        skeleton <- sapply(skeleton.names, function(x) NULL)
        
        if (all(curve[1,1:2] == curve[nrow(curve),1:2])) {
            np <- nrow(curve) - 1
        } else {
            np <- nrow(curve)
        }
        
        skeleton[["vertex"]] <- 1:np
        skeleton[["edge"]] <- vector()
        skeleton[["coords"]] <- vector()
        
        for (j in 1:np) {
            skeleton[["coords"]] <- rbind(skeleton[["coords"]],
                                          curve[j,1:2])
            skeleton[["edge"]] <- rbind(skeleton[["edge"]],
                                        c(j, (j+1)%%np))
        }
        
        complex[[i-1]] <- skeleton
    }
}

computeHeightFiltration <- function(curve,
                                    direction,
                                    orientation) {
    filtration.names <- c("height", "lowerNbrs", "coords", "minimal")
    filtration <- sapply(filtration.names, function(x) NULL)
    
    filtration[["coords"]] <- curve[["coords"]]
    
    filtration[["height"]] <- apply(curve[["coords"]], 1, function(v) dotProduct(v, direction))
    
    N <- max(curve[["vertex"]])
    
}

dotProuct <- function(v1, v2) {
    return(sum(v1 * v2))
}

testMinimality <- function(vertices,
                           heights,
                           direction,
                           colinearCond = 1e-8) {
    vk <- vertices[2,]
    vPrev <- vertices[1,]
    vNext <- vertices[3,]
    
    hk <- heights[2]
    hPrev <- heights[1]
    hNext <- heights[3]
    
    if (hk <= hPrev && hk <= hNext) {
        if (abs(hk - hPrev) <= colinearCond && abs(hk - hNext) <= colinearCond) {
            return(testNormalVector(vk, vNext, direction))
        } else {
            P <= c(0.5 * (hNext[1] + hPrev[1]),
                   0.5 * (hNext[2] + hPrev[2]))
            delta <- det(matrix(c(vk, 1, vNext, 1, P, 1),
                                nrow = 3,
                                ncol = 3,
                                byrow = TRUE))
            if (delta > 0) {
                return(TRUE)
            } else if (delta < 0) {
                return(FALSE)
            } else {
                stop("Cannot determine minimality.")
            }
        }
    }
}

testNormalVector <- function(v1,
                             v2,
                             direction) {
    normalVect <- c(v2[2] - v1[2], v2[1] - v1[1]) * c(-1,1)
    
    return(dotProduct(direction, normalVect) > 0)
}