#' Extract the boundary curves of a binary image
#' 
#' Given a binary image (pixel values 0 and 1) extractBoundary will trace all boundary curves along the
#' midpoints of pixel edges between foreground and background regions. Curves which bound the outside 
#' of a foreground region are traced in the anticlockwise direction while curves which bound the inside
#' of a foreground region are traced in the clockwise direction.
#' 
#' The image is padded with a 1-pixel wide border of background. Each connected component of the image
#' is labelled as follows:
#'      - Foreground pixels are labelled with positive integers (1,2,3,...)
#'      - Background pixels are labelled with negative integers (-1,-2,...)
#'      - The outermost background region is labelled with -1.
#' 
#' OUTPUT
#' ------
#' The boundary information is stored in the list boundary, structured as follows
#'      - boundary[[1]] is a vector (m,n) containing the number of anticlockwise and clockwise curves.
#'      - boundary[[2]] to boundary[[m+1]] contain the anticlockwise curves.
#'      - boundary[[m+2]] to boundary[[n+1]] contain the clockwise curves.
#' If set to save, the list is saved as an RDS file.
#' 
#' @param image A binary image read in from the imager package.
#' @param background The value of the background, 0 for black and 1 for white. (Default: 0)
#' @param saveOutput If TRUE, will save output to directory specified by outputDir. (Default: FALSE)
#' @param outputDir The directory to save the output. If saveOutput is TRUE and no directory is specified, saves to working directory. (Default: NULL)
#' @param fName The name of the output file saved. If saveOutput is TRUE and no filename specified, prompts user for filename. (Default: NULL)
#' @param verbose If TRUE, prints indictors of progress throughout. (Default: TRUE)
#' @return A list containing the the points around each boundary curve of the image as matrices.
#' @export
extractBoundary <- function(image,
                           background = 0,
                           saveOutput = FALSE,
                           outputDir  = NULL,
                           fName = NULL,
                           verbose = TRUE) {
    
    if (saveOutput) {
        if (!dir.exists(outputDir)) {
            outputDir <- getwd()
            if (!is.character(fName)) {
                fName <- readline(prompt="Please provide a filename for save: ")
            }
            if (nchar(fname) == 0) {
                stop("No filename provided. Aborting extractBoundary.")
            }
        }
    }
    # Pad image
    image <- imager::pad(image, nPix = 1, val = background)

    imgMatrix <- componentLabelling(image, background, verbose)
    if (any(imgMatrix > 1) || any(imgMatrix < 0)) {
        stop("Image not binary. Pixel values must be either 0 or 1.")
    }
    
    boundary <- boundaryTrace(imgMatrix, verbose)

    if (saveOutput) {
        if (substr(outputDir, nchar(outputDir), nchar(outputDir)) == "/") {
            outFile <- paste(outputDir, fname, ".RDS", sep = "")
        } else {
            outFile <- paste(outputDir, "/", fname, ".RDS", sep = "")
        }
        saveRDS(boundary, file = outFile)
    } else {
        return(boundary)
    }
}

#' Extract the boundary curves of binary images in a given directory
#' 
#' Runs extractBoundary on all images of a specified type in a given directory.
#' Images must be binary (pixel values 0 and 1). 
#'
#' Output is either returned as a list of lists or saved as individual RDS
#' files in the specified output directory.
#'
#' For more information, run help(extractBoundary)
#'
#' @param inputDir The directory containing the images.
#' @param imgType The image file type for input files, specified without a "." (Default: "png")
#' @param background The value of the background, 0 for black and 1 for white. (Default: 0)
#' @param saveOutput If TRUE, will save output to directory specified by outputDir. (Default: FALSE)
#' @param outputDir The directory to save the output. If saveOutput is TRUE and no directory is specified, saves to working directory. (Default: NULL)
#' @param verbose If TRUE, prints indictors of progress throughout. (Default: TRUE)
#' @return A list for each image containing the points of the boundary curves. Stored in individual files (if saveOutput = TRUE) or a list.
#' @export
multiExtractBoundary <- function(inputDir,
                                 imgType = 'png',
                                 background = 0,
                                 saveOutput = FALSE,
                                 outputDir = NULL,
                                 verbose = TRUE) {
                                     
    fType <- paste("*.", imgType, sep = "")
    files <- list.files(path = inputDir,
                        pattern = fType,
                        full.names = TRUE,
                        recursive = FALSE)
    
    if (saveOutput) {
        if (!dir.exists(outputDir)) {
            outputDir <- getwd()
        }
    } else {
        boundaries <- vector(mode = "list", length = length(files))
    }
    
    for (i in seq_along(files)) {
        if (verbose) {
            cat("Commencing", files[[i]], sep = " ")
        }
        
        image <- imager::load.image(files[[i]])
        
        if (saveOutput) {
            fString <- strsplit(files[[i]], "/", fixed = TRUE)
            fName <- fString[[1]][length(fString[[1]])]
            extractBoundary(image = image,
                            background = background,
                            saveOutput = TRUE,
                            outputDir = outputDir,
                            fName = fName,
                            verbose = verbose)
        } else {
            boundaries[[i]] <- extractBoundary(image = image,
                                               background = background,
                                               verbose = verbose)
        }
    }
    
    if(saveOutput) {
        cat("All boundaries successfully saved in:\n", outputDir, sep="")
    } else {
        return(boundaries)   
    }
}
componentLabelling <- function(image,
                               background,
                               verbose) {

    w <- imager::width(image)
    h <- imager::height(image)

    imgMatrix <- matrix(image, nrow = h, ncol = w, byrow = TRUE)
    imgMatrix <- apply(imgMatrix, 2, rev)
    
    equivLabelsFG <- vector(mode = "list")
    equivLabelsBG <- vector(mode = "list")
    
    currentFG <- 1
    currentBG <- -1
    
    foreground <- 1-background

    # image is scanned left-to-right, top-to-bottom
    for (i in 1:h) {
        for (j in 1:w) {
            if (imgMatrix[i,j] == foreground) {
                # Hit a foreground pixel. Assess the mask (X's):
                # X X X
                # X * o
                # o o o
                mask <- c(imgMatrix[i,j-1],
                          imgMatrix[i-1,j-1],
                          imgMatrix[i-1,j],
                          imgMatrix[i,j+1])
                
                if (mask[3] > 0) {
                    imgMatrix[i,j] <- mask[3]
                } else if (mask[1] > 0) {
                    imgMatrix[i,j] <- mask[1]
                    if (mask[4] > 0) {
                        equivLabelsFG <- resolveLabels(mask[4], mask[1], equivLabelsFG, bg = FALSE)
                    }
                } else if (mask[2] > 0) {
                    imgMatrix[i,j] <- mask[2]
                    if (mask[4] > 0) {
                        equivLabelsFG <- resolveLabels(mask[2], mask[4], equivLabelsFG, bg = FALSE)
                    }
                } else if (mask[4] > 0) {
                    imgMatrix[i,j] <- mask[4]
                } else {
                    imgMatrix[i,j] <- currentFG
                    equivLabelsFG[[currentFG]] <- currentFG
                    currentFG <- currentFG + 1
                }
            } else {
                # Found a background pixel
                mask <- c(0,0)
                if (i > 1) {
                    mask[2] <- imgMatrix[i-1,j]
                }
                if (j > 1) {
                    mask[1] <- imgMatrix[i,j-1]
                }
                
                if (mask[2] < 0) {
                    imgMatrix[i,j] <- mask[2]
                    if (mask[1] < 0) {
                        equivLabelsBG <- resolveLabels(-mask[1], -mask[2], equivLabelsBG, bg = TRUE)
                    }
                } else if (mask[1] < 0) {
                    imgMatrix[i,j] <- mask[1]
                } else {
                    imgMatrix[i,j] <- currentBG
                    x <- -currentBG
                    equivLabelsBG[[x]] <- currentBG
                    currentBG <- currentBG - 1
                }
            } 
        }
    }
    if (verbose) {
        print("Completed first raster for component labelling. Now constructing representative tables.")
    }
    
    for (i in 1:length(equivLabelsFG)) {
        for (lab in equivLabelsFG[[i]]) {
            equivLabelsFG[[lab]] <- unique(c(equivLabelsFG[[lab]], equivLabelsFG[[i]]))
        }
    }
    
    for (i in 1:length(equivLabelsBG)) {
        for (lab in equivLabelsBG[[i]]) {
            x <- -lab
            equivLabelsBG[[x]] <- unique(c(equivLabelsBG[[x]], equivLabelsBG[[i]]))
        }
    }
    
    repTableFG <- vector()
    repTableBG <- vector()
    
    for (i in 1:length(equivLabelsFG)) {
        repTableFG <- append(repTableFG, min(equivLabelsFG[[i]]))
    }
    for (i in 1:length(equivLabelsBG)) {
        repTableBG <- append(repTableBG, max(equivLabelsBG[[i]]))
    }
    
    lab <- 1
    for (i in 1:(currentFG - 1)) {
        if (repTableFG[i] == i) {
            repTableFG[i] <- lab
            j <- j + 1
        } else {
            repTableFG[i] <- repTableFG[repTableFG[i]]
        }
    }
    
    lab <- -1
    for (i in 1:(abs(currentBG)-1)) {
        if (repTableBG[i] == -i) {
            repTablBG <- lab
            lab <- lab - 1
        } else {
            x <- -repTableBG[i]
            repTableBG[i] <- repTableBG[x]
        }
    }
    
    if (verbose) {
        print("Representative tables constructed. Now relabelling components.")
    }
    
    for (i in 1:h) {
        for (j in 1:w) {
            if (imgMatrix[i,j] > 0) {
                imgMatrix[i,j] <- repTableFG[imgMatrix[i,j]]
            } else  {
                x <- -imgMatrix[i,j]
                imgMatrix[i,j] <- repTableBG[x]
            }
        }
    }
    
    return(imgMatrix)
}
        
resolveLabels <- function(lab1,
                          lab2,
                          labList,
                          bg) {
    if (bg) {
        labList[[lab1]] <- unique(c(labList[[lab1]], labList[[lab2]], -lab2))
        labList[[lab2]] <- unique(c(labList[[lab1]], labList[[lab2]], -lab1))
    }  else {
        labList[[lab1]] <- unique(c(labList[[lab1]], labList[[lab2]], lab1, lab2))
        lablist[[lab2]] <- unique(c(labList[[lab1]], labList[[lab2]], lab1, lab2))
    }
    return(labList)
}

boundaryTrace <- function(imgMatrix) {
    
    w <- ncol(imgMatrix) - 1
    h <- hrow(imgMatrix) - 1
    
    startPx <- vector()
    regions <- vector()
    
    maxRegionLabels <- c(min(imgMatrix), max(imgMatrix))
    final <- FALSE
    
    cCurves <- vector(mode = "list") # clockwise curve, i.e. inner
    acCurves <- vector(most = "list") # anticlockwise curve, i.e. outer.
    
    for (i in 2:h) {
        for (i in 2:w) {
            if (imgMatrix[i,j] != -1) {
                if (!any(regions == imgMatrix[i,j])) {
                    lab <- imgMatrix[i,j]
                    if (lab > 0) {
                        acCurves[[lab]] <- traceCurve(imgMatrix, i, j, lab)
                    } else {
                        cCurves[[-lab]] <- traceCurve(imgMatrix, i, j, lab)
                    }
                    regions <- append(regions,
                                      imgMatrix[i,j])
                }
                if (any(regions == maxRegionLabels[1]) && any(regions == maxRegionLabels[2])) {
                    final <- TRUE
                    break
                }
            }
        }
        if (final) {
            break
        }
    }
    
    totalCurves <- length(cCurves) + length(acCurves)
    
    if (verbose) {
        cat("Successfully extracted", totalCurves, ":\n", length(acCurves), "anticlcokwise curves\n", length(cCurves), "clockwise curves", sep = " ")
    }
    
    boundaryCurves <- vector(mode = "list")
    boundaryCurves[[1]] <- c(length(acCurves), length(cCurves))
    boundaryCurves <- append(boundaryCurves, acCurves, cCurves)
    
    return(boundaryCurves)
}

traceCurve <- function(imgMatrix,
                       i,
                       j,
                       lab) {
    p0 <- c(i-0.5, j) # start with North egde
    
    p1 <- c(i,j)
    pc <- c(0,0)
    pxLoc <- 1 # N = 1, W = 2, S = 3, E = 4
    
    curvePts <- matrix(data = p0, nrow = 1, ncol = 2)
    
    dir <- 7
    prevDir <- 0
    
    while (any(pc != p1)) {
        nbhd <- matrix(c(i, i-1, i-1, i-1, i, i+1, i+1, i+1,
                         j+1, j+1, j, j-1, j-1, j-1, j, j+1),
                       nrow = 8,
                       ncol = 2)
        
        if (dir %% 2 == 1) {
            dir <- dir - 1
        }
        
        for (k in ((dir+7):(dir+14) %% 8)) {
            nbr <- FALSE
            x <- nbhd[k+1,][1]
            y <- nbhd[k+1,][2]
            
            if (imgMatrix[x,y] == lab && any(c(x,y) != pc)) {
                nbr <- TRUE
                prev <- switch(pxLoc,
                               switch(k+1,
                                      rbind(c(i,j-0.5),
                                            c(i+0.5,j),
                                            c(x+0.5,y)),
                                      rbind(c(i,j-0.5),
                                            c(i+0.5,j),
                                            c(i,j+0.5),
                                            c(x+0.5,y)),
                                      c(),
                                      c(x,y+0.5),
                                      c(x-0.5,y),
                                      rbind(c(i,j-0.5),
                                            c(x-0.5,y)),
                                      rbind(c(i,j-0.5),
                                            c(x,y-0.5)),
                                      rbind(c(i,j-0.5),
                                            c(i+0.5,j),
                                            c(x,y-0.5))),
                               switch(k+1,
                                      rbind(c(i+0.5,j),
                                            c(x+0.5,y)),
                                      rbind(c(i+0.5,j),
                                            c(i,j+0.5),
                                            c(x+0.5,y)),
                                      c(x+0.5,y),
                                      rbind(c(i+0.5,j),
                                            c(i,j+0.5),
                                            c(i-0.5,j),
                                            c(x,y+0.5)),
                                      c(),
                                      c(x-0.5,y),
                                      c(x,y-0.5),
                                      rbind(c(i+0.5,j),
                                            c(x,y-0.5))),
                               switch(k+1,
                                      c(x+0.5,y),
                                      rbind(c(i,j+0.5),
                                            c(x+0.5,y)),
                                      rbind(c(i,j+0.5),
                                            c(x,y+0.5)),
                                      rbind(c(i,j+0.5),
                                            c(i-0.5,j),
                                            c(x,y+0.5)),
                                      rbind(c(i,j+0.5),
                                            c(i-0.5,j),
                                            c(x-0.5,y)),
                                      rbind(c(i,j+0.5),
                                            c(i-0.5,j),
                                            c(i,j-0.5),
                                            c(x-0.5,y)),
                                      c(),
                                      c(x,y-0.5)),
                               switch(k+1,
                                      c(),
                                      c(x+0.5,y),
                                      c(x,y+0.5),
                                      rbind(c(i-0.5,j),
                                            c(x,y+0.5)),
                                      rbind(c(i-0.5,j),
                                            c(x-0.5,y)),
                                      rbind(c(i-0.5,j),
                                            c(i,j-0.5),
                                            c(x-0.5,y)),
                                      rbind(c(i-0.5,j),
                                            c(i,j-0.5),
                                            c(x,y-0.5)),
                                      rbind(c(i-0.5,j),
                                            c(i,j-0.5),
                                            c(i+0.5,j),
                                            c(x,y-0.5))))
                
                curvePts <- unname(rbind(curvePts, prev))
                pxLoc <- switch(pxLoc,
                                swich(k+1,
                                      3,
                                      3,
                                      NULL,
                                      4,
                                      1,
                                      1,
                                      2,
                                      2),
                                switch(k+1,
                                       3,
                                       3,
                                       4,
                                       4,
                                       NULL,
                                       1,
                                       2,
                                       2),
                                switch(k+1,
                                       3,
                                       3,
                                       4,
                                       4,
                                       1,
                                       1,
                                       NULL,
                                       2),
                                switch(k+1,
                                       NULL,
                                       3,
                                       4,
                                       4,
                                       1,
                                       1,
                                       2,
                                       2))
                dir <- k1
                i <- x
                j <- y
                pc <- c(i,j)
                break                
            }
        }
        if (!nbr) {
            # The region is a single pixel
            curvePts <- rbind(curvePts,
                              c(i,j-0.5),
                              c(i+0.5,j),
                              c(i,j+0.5))
            pc <- p1
        }
    }
    
    if (lab > 0) {
        curvePts[,1] <- rev(curvePts[,1])
        curvePts[,2] <- rev(curvePts[,2])
    }
    
    return(curvePts)
}