#' Compute the Extended PHT of a given image
#'
#' Given a disjoint collection of simple closed curves representing the boundaries of
#' a two dimensional shape, compute the zeroth order extended persistent homology transform
#' of the shape in a given number of directions.
#' 
#' Computing the extended persistent homology in a given direction requires information
#' about the orientation of a curve. This information is used in detemining whether a vertex
#' is a local minimum with respect to the foreground of the original image.
#'
#' @param bdryCurves A list containing points along the boundary curves of the image. The list
#' must have length number_of_curves + 1 with the first entry a vector whose entries are the numbers
#' of positively oriented curves and negatively oriented curves, respectively. All of the positively 
#' oriented curves must be listed before the negatively oriented curves.
#' @param imgName A string giving the name of the image.
#' @param nDirections An integer giving the number of directions to use. This must be positive and even.
#' @param tolerance Parameter to control the noise from the image. Any vertices closer than this in a particular
#' direction will be classified as noise and ignored. (Default: 1/sqrt(2))
#' @param saveOutput If TRUE, will save output to directory specified by outputDir. (Default: FALSE)
#' @param outputDir The directory to save the output. If saveOutput is TRUE and no directory is specified, saves to working directory. (Default: NULL)
#' @param fName The name of the output file saved. If saveOutput is TRUE and no filename specified, prompts user for filename. (Default: NULL)
#' @param verbose If TRUE, prints indictors of progress throughout. (Default: TRUE)
#' @return A list of length 5 with names "Name", "Ord0", "Rel1", "Ext0", "Ext1". Each entry is a list
#' containing the relevant part of the extended persistence diagram. The i-th entry corresponds to the i-th direction.
#' @export
extendedPersistence <- function(bdryCurves,
                                imgName,
                                nDirections,
                                tolerance = 1/sqrt(2),
                                saveOutput = FALSE,
                                outputDir = NULL,
                                fName = NULL,
                                verbose = TRUE) {
    if (saveOutput) {
        if (!dir.exists(outputDir)) {
            outputDir <- getwd()
            print("Output directory doesn't exist. Saving to working directory.")
        }
        if (is.null(fName)) {
            fName <- readLine(prompt="Please provide a filename for save: ")
        }                          
        outFile <- paste(outputDir, "/", fName, ".RDS", sep = "")
    }
    if (nDirections < 2 || nDirections %% 2 == 1) {
        stop("Number of directions must be non-zero and even")
    } else {
        directions <- t(sapply(1:nDirections, function(r) c(cos(2 * pi * r / nDirections),
                                                            sin(2 * pi * r / nDirections))))
        midPoint <- (nDirections/2) + 1
    }
    
    #xDiagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
    #xDiagram <- sapply(xDiagram.names, function(x) NULL)
    xDiagram <- vector(mode = "list")
    
    oneSkeleton <- parseSkeleton(bdryCurves)
    
    acCurves <- bdryCurves[[1]][1]
    cCurves <- bdryCurves[[1]][2]
    
    for (d in 1:(nDirections/2)) {
        diagramName <- paste(imgName, "-", toString(d), sep="")
        #xDiagram[["Name"]][[d]] <- diagramName
        
        dirVector <- directions[d,]
        
        nComponents <- length(oneSkeleton)
        zerothDiagram <- vector(mode = "list", length = nComponents)
        
        count <- 1
        
        for (i in 1:nComponents) {
            skeletonComponent <- oneSkeleton[i][[1]]
            if (i <= acCurves) {
                orientation = 0
            } else {
                orientation = 1
            }
            heightFiltration <- computeHeightFiltration(curve = skeletonComponent,
                                                        direction = dirVector,
                                                        orientation = orientation)
            
            zerothDiagram[[i]] <- computeDiagram(filtration = heightFiltration,
                                                 tolerance = tolerance)
        }
        extendedDiagram <- computeExtendedDiagram(zerothDiagram = zerothDiagram,
                                                  nComponents = nComponents,
                                                  diagramName = diagramName)
        
        #xDiagram[["Ord0"]][[d]] <- extendedDiagram[["Ord0"]]
        #xDiagram[["Rel1"]][[d]] <- extendedDiagram[["Rel1"]]
        #xDiagram[["Ext0"]][[d]] <- extendedDiagram[["Ext0"]]
        #xDiagram[["Ext1"]][[d]] <- extendedDiagram[["Ext1"]]
        class(extendedDiagram) <- "extDiagram"
        xDiagram[[d]] <- extendedDiagram
        
    }
    
    for (d in midPoint:nDirections) {
        diagramName <- paste(imgName, "-", toString(d), sep="")
        #xDiagram[["Name"]][[d]] <- diagramName
        negDiagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
        negDiagram <- sapply(negDiagram.names, function(x) NULL)
        
        negDiagram[["Name"]] <- diagramName
        k <- d - midPoint + 1
        
        #xDiagram[["Ord0"]][[d]] <- -xDiagram[["Rel1"]][[k]]
        #xDiagram[["Rel1"]][[d]] <- -xDiagram[["Ord0"]][[k]]
        negDiagram[["Ord0"]] <- -xDiagram[[k]][["Rel1"]]
        negDiagram[["Rel1"]] <- -xDiagram[[k]][["Ord0"]]
        negDiagram[["Ext0"]] <- -matrix(xDiagram[[k]][["Ext0"]][,c(2,1)], ncol = 2)
        
        #xDiagram[["Ext0"]][[d]] <- -matrix(xDiagram[["Ext0"]][[k]][,c(2,1)], ncol=2)
        
        if (length(xDiagram[[k]][["Ext1"]]) > 0) {
            #xDiagram[["Ext1"]][[d]] <- -matrix(xDiagram[["Ext1"]][[k]][,c(2,1)], ncol=2)
            negDiagram[["Ext1"]] <- -matrix(xDiagram[[k]][["Ext1"]][,c(2,1)], ncol = 2)
        } else {
            #xDiagram[["Ext1"]][[d]] <- vector()
            negDiagram[["Ext1"]] <- vector()
        }
        class(negDiagram) <- "extDiagram"
        xDiagram[[d]] <- negDiagram
    }
    
    class(xDiagram) <- "extDiagram"
    
    if (saveOutput) {
        saveRDS(xDiagram, file = outFile)
        if (verbose) {
            cat("Successfully saved ", outFile, "\n", sep = "")
        }
    } else {
        if (verbose) {
            cat("Extended persistence diagrams successfully computed for",
                nDirections, "directions.\n", sep = " ")
        }
        return(xDiagram)
    }
}

#' Compute the Extended PHT for multiple images
#'
#' Runs \code{\link{extractBoundary}} on all .RDS files in a given directory.
#' 
#' Each persistence diagram is \emph{saved} in in the specified output directory.
#'
#' For more information, see \code{\link{extractBoundary}}.
#'
#' @param inputDir The directory containing the extracted boundary curves.
#' @param outputDir The directory to save the output. If saveOutput is TRUE and no directory is specified, saves to working directory. (Default: NULL)
#' @param nDirections An integer giving the number of directions to use. This must be positive and even.
#' @param tolerance Parameter to control the noise from the image. Any vertices closer than this in a particular
#' direction will be classified as noise and ignored. (Default: 1/sqrt(2))
#' @param verbose If TRUE, prints indictors of progress throughout. (Default: TRUE)
#' @export
multiExtendedPersistence <- function(inputDir,
                                     outputDir,
                                     nDirections,
                                     tolerance = 1/sqrt(2),
                                     verbose = TRUE) {
    if (!dir.exists(outputDir)) {
        outputDir <- getwd()
        print("Output directory doesn't exist. Saving to working directory.")
    }
        
    files <- list.files(path = inputDir,
                        pattern = "*.RDS",
                        full.names = TRUE,
                        recursive = FALSE)
    
    for (i in seq_along(files)) {
        f <- files[[i]]
        if (verbose) {
            cat("Commencing,", f, "\n", sep = " ")
        }
        
        fName <- tail(strsplit(strsplit(f, ".", fixed = TRUE)[[1]][1],
                               "/", fixed = TRUE)[[1]], n = 1)
        
        bdryCurves <- readRDS(f)
        
        extendedPersistence(bdryCurves = bdryCurves,
                            imgName = fName,
                            nDirections = nDirections,
                            tolerance = tolerance,
                            saveOutput = TRUE,
                            outputDir = outputDir,
                            fName = fName,
                            verbose = verbose)
    }
}
                                     
parseSkeleton <- function(bdryCurves) {
    
    complex <- vector(mode = "list")
    
    skeleton.names <- c("vertex", "edge", "coords")
    
    acCurves <- bdryCurves[[1]][1]
    cCurves <- bdryCurves[[1]][2]
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
        
        for (j in 1:(np-1)) {
            skeleton[["coords"]] <- rbind(skeleton[["coords"]],
                                          curve[j,1:2])
            skeleton[["edge"]] <- rbind(skeleton[["edge"]],
                                        c(j, j+1))
        }
        skeleton[["edge"]] <- rbind(skeleton[["edge"]],
                                    c(np, 1))
        skeleton[["coords"]] <- rbind(skeleton[["coords"]],
                                      curve[np,1:2])
        
        complex[[i-1]] <- skeleton
    }
    
    return(complex)
}

computeHeightFiltration <- function(curve,
                                    direction,
                                    orientation) {
    filtration.names <- c("height", "lowerNbrs", "coords", "minimal")
    filtration <- sapply(filtration.names, function(x) NULL)
    
    filtration[["coords"]] <- curve[["coords"]]
    
    filtration[["height"]] <- apply(curve[["coords"]], 1, function(v) dotProduct(v, direction))
    
    N <- max(curve[["vertex"]])
    filtration[["minimal"]] <- vector(length = N)
    
    vtxs <- rbind(curve[["coords"]][N,],
              curve[["coords"]][1,],
              curve[["coords"]][2,])
    hs <- c(filtration[["height"]][N],
            filtration[["height"]][1],
            filtration[["height"]][2])
    
    filtration[["minimal"]][1] <- testMinimality(vtxs, hs, direction, orientation)
    
    vtxs <- rbind(curve[["coords"]][N-1,],
              curve[["coords"]][N,],
              curve[["coords"]][1,])
    hs <- c(filtration[["height"]][N-1],
            filtration[["height"]][N],
            filtration[["height"]][1])
    
    filtration[["minimal"]][N] <- testMinimality(vtxs, hs, direction, orientation)
    
    for (i in 2:(N-1)) {
        vtxs <- rbind(curve[["coords"]][i-1,],
                curve[["coords"]][i,],
                curve[["coords"]][i+1,])
        hs <- c(filtration[["height"]][i-1],
                filtration[["height"]][i],
                filtration[["height"]][i+1])
    
        filtration[["minimal"]][i] <- testMinimality(vtxs, hs, direction, orientation)
    }
    
    filtration[["lowerNbrs"]] <- vector(mode = "list", length = N)
    
    for (i in 1:N) {
        e <- curve[["edge"]][i,]
        
        h1 <- filtration[["height"]][e[1]]
        h2 <- filtration[["height"]][e[2]]
        
        if (h1 < h2) {
            filtration[["lowerNbrs"]][[e[2]]] <- append(filtration[["lowerNbrs"]][[e[2]]],
                                                        e[1])
        } else if (h1 == h2) {
            # Vertex with lower index is lower neighbour.
            if (e[1] < e[2]) {
                filtration[["lowerNbrs"]][[e[2]]] <- append(filtration[["lowerNbrs"]][[e[2]]],
                                                            e[1])
            } else {
                filtration[["lowerNbrs"]][[e[1]]] <- append(filtration[["lowerNbrs"]][[e[1]]],
                                                            e[2])
            }
        } else {
            filtration[["lowerNbrs"]][[e[1]]] <- append(filtration[["lowerNbrs"]][[e[1]]],
                                                        e[2])
        }
    }
    
    return(filtration)
}

dotProduct <- function(v1, v2) {
    return(sum(v1 * v2))
}

testMinimality <- function(vertices,
                           heights,
                           direction,
                           orientation,
                           colinearCond = 1e-8) {
    vk <- vertices[2,]
    vPrev <- vertices[1,]
    vNext <- vertices[3,]
    
    hk <- heights[2]
    hPrev <- heights[1]
    hNext <- heights[3]
    
    if (hk > hPrev || hk > hNext) {
        return(FALSE)
    }
    
    if (abs(hk - hPrev) <= colinearCond && abs(hk - hNext) <= colinearCond) {
        return(testNormalVector(vk, vNext, direction))
    } else {
        P <- c(0.5 * (vNext[1] + vPrev[1]),
               0.5 * (vNext[2] + vPrev[2]))

        delta <- (vNext[1] - vk[1])*(P[2] - vk[2]) - (vNext[2] - vk[2])*(P[1] - vk[1])
            
        if (delta != 0) {
            return(delta > 0)
        } else {
            stop("Numerical error in minimality test. Could not determine minimality.")
        }
    }
}

testNormalVector <- function(v1,
                             v2,
                             direction) {
    normalVect <- c(v2[2] - v1[2], v2[1] - v1[1]) * c(-1,1)
    
    return(dotProduct(direction, normalVect) > 0)
}

computeDiagram <- function(filtration,
                           tolerance) {
    sortedHeights <- unique(sort(filtration[["height"]]))
    
    diagram.names <- c("finite", "extended", "minimal", "exMinimal")
    diagram <- sapply(diagram.names, function(x) NULL)
    
    parents <- vector(mode = "list",
                      length = length(filtration[["height"]]))
    
    diagram[["extended"]] <- c(head(sortedHeights,1),
                               tail(sortedHeights,1))
    
    for (hv in sortedHeights) {
        
        hVertex <- which(filtration[["height"]] == hv) #all vertices at height hv
        
        for (v in hVertex) {
            
            if (is.null(filtration[["lowerNbrs"]][[v]])) {
                # The vertex v has no lower neighbours. 
                parents[[v]] <- v
            } else {
                
                components <- sapply(filtration[["lowerNbrs"]][[v]],
                                     function(x) findParent(x, parents))
                components <- unique(components) #ignore multiple paths to same parent
                
                if (length(components) == 1) {
                    parents[[v]] <- components
                } else {
                    birthTimes <- sapply(components,
                                         function(x) filtration[["height"]][x])
                    minBirthTime <- min(birthTimes)
                    
                    components <- sort(components)
                    
                    # If there are two candidates for the birth of a component,
                    # we take the one with the lower index.
                    # This is for consistency.
                    count = 0
                    
                    for (x in components) {
                        hx <- filtration[["height"]][x]
                        if (hx > minBirthTime) {
                            if (hx < hv) {
                                # Component born at height of x dies at 
                                # the current height.
                                if (abs(hx - hv) > tolerance) {
                                    diagram[["finite"]] <- rbind(diagram[["finite"]],
                                                                 c(hx, hv))
                                    diagram[["minimal"]] <- append(diagram[["minimal"]],
                                                                   filtration[["minimal"]][x])
                                }
                            }
                        } else if (hx == minBirthTime) {
                            if (count == 0) {
                                newComponent <- x
                                parents[[v]] <- newComponent
                                count <- 1
                            } else {
                                if (hx < hv) {
                                    if (abs(hx - hv) > tolerance) {
                                        diagram[["finite"]] <- rbind(diagram[["finite"]],
                                                                     c(hx, hv))
                                        diagram[["minimal"]] <- append(diagram[["minimal"]],
                                                                       filtration[["minimal"]][x])
                                    }
                                }
                            }
                        }
                    }
                    # All components found are part of same connected component.
                    # Update this information
                    for (x in components) {
                        parents[[x]] <- newComponent
                    }
                }
            }
        }
    }
    
    birthPoint <- unique(sapply(unique(parents),
                                function(x) findParent(x, parents)))
    
    if (length(birthPoint) > 1) {
        stop("Simple closed curve has more than one essential class.")
    }
    
    diagram[["exMinimal"]] <- filtration[["minimal"]][[birthPoint[1]]]
    
    return(diagram)
}

findParent <- function(x,
                       parents) {
    if (parents[[x]] == x) {
        return(x)
    } else {
        px <- findParent(parents[[x]], parents)
        return(px)
    }
}

computeExtendedDiagram <- function(zerothDiagram,
                                   nComponents,
                                   diagramName) {
    exDiagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
    exDiagram <- sapply(exDiagram.names, function(x) NULL)
    
    Ord0 <- vector()
    Rel1 <- vector()
    Ext0 <- vector()
    Ext1 <- vector()
    
    exDiagram[["Name"]] <- diagramName
    for (i in 1:nComponents) {
        diagram <- zerothDiagram[[i]]
        
        if (length(diagram[["finite"]]) > 0) {
            nFinite <- nrow(diagram[["finite"]])
            
            for (j in 1:nFinite) {
                if (diagram[["minimal"]][j]) {
                    Ord0 <- rbind(Ord0,
                                  diagram[["finite"]][j,])
                } else {
                    Rel1 <- rbind(Rel1,
                                  rev(diagram[["finite"]][j,]))
                }
            }
        }
        
        if (diagram[["exMinimal"]]) {
            Ext0 <- rbind(Ext0,
                          diagram[["extended"]])
        } else {
            Ext1 <- rbind(Ext1,
                          rev(diagram[["extended"]]))
        }
    }
    
    #class(Ord0) <- "extDiagram"
    #class(Rel1) <- "extDiagram"
    #class(Ext0) <- "extDiagram"
    #class(Ext1) <- "extDiagram"
    
    exDiagram[["Ord0"]] <- Ord0
    exDiagram[["Rel1"]] <- Rel1
    exDiagram[["Ext0"]] <- Ext0
    exDiagram[["Ext1"]] <- Ext1
    
    #attributes(exDiagram)[["extended"]] <- TRUE        
    
    return(exDiagram)                        
}