#' Centre and Scale the Extended Persistence Diagrams of an image.
#'
#' Given a collection of persistence diagrams in N evenly spaced directions, shift the points
#' in each diagram so that the correspond to the points of the persistence diagrams of the
#' centred image. Then scale the birth and death times of each point by some factor.
#'
#' Centring and scaling the persistence diagrams allows for a more accurate comparison of their shapes.
#'
#' @param diagrams A list of extended persistence diagrams on a single image computed in N evenly spaced directions.
#' @param scale_constant A constant to control the scaling of the diagram. Default is 1. Note this is not the
#' percentage of the original scale of the image.
#' @return A list of extended persistence diagrams of the same image, centred at the origin and scaled by an amount proportional
#' to the sum of the minimum birth times.
#' @export
centreScaleDiagrams <- function(diagrams,
                                scale = TRUE,
                                scale_constant = 1) {
    # Diagram has fields [['name']], [['Ord0']], [['Rel1']], [['Ext0']], [['Ext1']]
    nDirs <- length(diagrams)
    lambda <- minBirthTimes(nDirs = nDirs,
                            diagrams = diagrams)
    
    csDiagrams <- centreDiagrams(diagrams = diagrams,
                                 nDirs = nDirs,
                                 lambda = lambda)
    print("Diagrams successfully centred.")
    if (scale) {
        csDiagrams <- scaleDiagrams(diagrams = diagrams,
                                    nDirs = nDirs,
                                    lambda = lambda,
                                    scale_constant = scale_constant)
        print("Diagrams successfully scaled.")
    }
    
    return(csDiagrams)
      
}

minBirthTimes <- function(nDirs,
                          diagrams) {
    lambda <- vector()
    # The first component born in any direction belongs to Ext0.
    for (i in 1:nDirs) {
        birthTimes <- c(diagrams[[i]][["Ext0"]][,1])
        lambda <- append(lambda, min(birthTimes))
    }                              
    
    return(lambda)
}

centreDiagrams <- function(diagrams,
                           nDirs,
                           lambda) {
                           
    cp <- findCentre(nDirs = nDirs,
                     lambda = lambda)
    
    directions <- t(sapply(1:nDirs, function(r) c(cos(2 * pi * r / nDirs),
                                                  sin(2 * pi * r / nDirs))))
    
    for (i in 1:nDirs) {
        shift <- sum(cp, directions[i,])
        
        diagrams[[i]][["Ext0"]] <- diagrams[[i]][["Ext0"]] - shift
        
        if (length(diagrams[[i]][["Ord0"]]) > 0) {
            diagrams[[i]][["Ord0"]] <- diagrams[[i]][["Ord0"]] - shift
        }
        
        if (length(diagrams[[i]][["Rel1"]]) > 0) {
            diagrams[[i]][["Rel1"]] <- diagrams[[i]][["Rel1"]] - shift
        }
        
        if (length(diagrams[[i]][["Ext1"]]) > 0) {
            diagrams[[i]][["Ext1"]] <- diagrams[[i]][["Ext1"]] - shift
        }
    }
    
    return(diagrams)             
}

findCentre <- function(nDirs,
                       lambda) {
    directions <- t(sapply(1:nDirs, function(r) c(cos(2 * pi * r / nDirs),
                                                        sin(2 * pi * r / nDirs))))
    
    M <- lambda * directions
    
    cp <- colSums(M) / (nDirs * 0.5)
    
    return(cp)                       
}

scaleDiagrams <- function(diagrams,
                          nDirs,
                          lambda,
                          scale_constant) {
    L <- (-1)*sum(lambda)
    
    if (L <= 0) {
        stop("Incorrect lambda's computed. Cannot scale by a negative value.")
    }
    
    scale_value <- scale_constant / L
    
    for (i in 1:nDirs) {
        diagrams[[i]][["Ext0"]] <- diagrams[[i]][["Ext0"]] * scale_value
        
        if (length(diagrams[[i]][["Ord0"]]) > 0) {
            diagrams[[i]][["Ord0"]] <- diagrams[[i]][["Ord0"]] * scale_value
        }
        
        if (length(diagrams[[i]][["Rel1"]]) > 0) {
            diagrams[[i]][["Rel1"]] <- diagrams[[i]][["Rel1"]] * scale_value
        }
        
        if (length(diagrams[[i]][["Ext1"]]) > 0) {
            diagrams[[i]][["Ext1"]] <- diagrams[[i]][["Ext1"]] * scale_value
        }
    }
    
    return(diagrams)
}