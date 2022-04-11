centreScaleDiagrams <- function(diagrams,
                                scaleConstant = 1) {
    
    nDirs <- length(diagrams)
    lambda <- minBirthTimes(nDirs = nDirs,
                            diagrams = diagrams)
    
    csDiagrams <- centreDiagrams(diagrams = diagrams,
                                 nDirs = nDirs,
                                 lambda = lambda)       
}

minBirthTimes <- function(nDirs,
                          diagrams) {
    lambda <- vector()
    
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
    
                                 
}

findCentre <- function(nDirs,
                       lambda) {
    directions <- t(sapply(1:nDirs, function(r) c(cos(2 * pi * r / nDirs),
                                                        sin(2 * pi * r / nDirs))))
    
    M <- lambda * directions
    
    cp <- colSums(M) / (nDirs * 0.5)
    
    return(cp)                       
}