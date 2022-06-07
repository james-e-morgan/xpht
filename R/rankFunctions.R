###################################
######### RANK FUNCTIONS ##########
###################################
# Author: Yossi Bokor
# Version: 1
# Date: May 30, 2022
#' discretisedRankAbove
#'
#' Evaluate the a discretised rank function, evaluated on the cartesian product of xcoord and ycoord.
#'
#' @param diagram A set of n points above the diagonal, stored as an n x 2 array.
#' @param xcoords A list of x-coordinates
#' @param ycoords A list of y-coordinates
#' @return An array containing the values of the rank function, where the (i,j)-entry is the value at the point (x_i, y_j).
#'
#' @export

discretisedRankAbove <- function(diagram, xcoords,ycoords){
  nx = length(xcoords)
  ny = length(ycoords)
  np = length(diagram[,1])
  rk = matrix(0,nrow=ny, ncol=nx)

  for (i in 1:np){
    xi = Position(function(x) x >= diagram[i,1], xcoords)
    yi = Position(function(y) y >= diagram[i,2], ycoords)
    xinc = c(xi:nx)
    yinc = c(yi:nx)
    rk[xinc,yinc] <- rk[xinc,yinc]+matrix(1,nrow=length(xinc),ncol=length(yinc))
    print(rk[xinc,yinc])
  }

  for (i in 1:length(xcoords)){
    for (j in 1:length(ycoords)){
      if (xcoords[i]>ycoords[j]){
        rk[i,j] <- 0
      }
    }
  }
  return(rk)
}

#' discretisedRankBelow
#'
#' Evaluate the a discretised rank function, evaluated on the cartesian product of xcoord and ycoord.
#'
#' @param diagram A set of n points below the diagonal, stored as an n x 2 array.
#' @param xcoords A list of x-coordinates
#' @param ycoords A list of y-coordinates
#' @return An array containing the values of the rank function, where the (i,j)-entry is the value at the point (x_i, y_j).
#'
#' @export

discretisedRankBelow <- function(diagram, xcoords,ycoords){
  nx = length(xcoords)
  ny = length(ycoords)
  np = length(diagram[,1])
  rk = matrix(0,nrow=ny, ncol=nx)

  for (i in 1:np){
    xi = Position(function(x) x <= diagram[i,1], xcoords)
    yi = Position(function(y) y <= diagram[i,2], ycoords)
    xinc = c(xi:nx)
    yinc = c(yi:nx)
    rk[xinc,yinc] <- rk[xinc,yinc]+matrix(1,nrow=length(xinc),ncol=length(yinc))
    print(rk[xinc,yinc])
  }

  for (i in 1:length(xcoords)){
    for (j in 1:length(ycoords)){
      if (xcoords[i]>ycoords[j]){
        rk[i,j] <- 0
      }
    }
  }
  return(rk)
}

#' discretisedRank
#'
#' Evaluate the a discretised rank function, evaluated on the cartesian product of xcoord and ycoord.
#'
#' @param diagram A persistence diagam
#' @param xcoords A list of x-coordinates
#' @param ycoords A list of y-coordinates
#' @param name Name of the diagram
#' @return An array containing the values of the rank function, where the (i,j)-entry is the value at the point (x_i, y_j).
#'
#' @export

discretisedRank <- function(diagram, xcoords,ycoords, name){
  if (name in c('Ord0', 'Ext0')) {
    rk <- discretisedRankAbove(diagram, xcoords, ycoords)
  } else if (name in c('Rel1', 'Ext1')) {
    rk <- discretisedRankAbove(diagram, xcoords, ycoords)
  } else {
    # THROW ERROR NAME NOT A CORRECT DIAGRAM
  }
  return(rk)
}

#' discretisedRankAll
#'
#' Evaluate the a discretised rank function, evaluated on the cartesian product of xcoord and ycoord.
#'
#' @param diagram A persistence diagam
#' @param xcoords A list of x-coordinates
#' @param ycoords A list of y-coordinates
#' @param name Name of the diagram
#' @return An array containing the values of the rank function, where the (i,j)-entry is the value at the point (x_i, y_j).
#'
#' @export

discretisedRankAll <- function(diagram, xcoords,ycoords) {
  rks <- c()
  for (name in c('Ord0', 'Ext0', 'Rel1', 'Ext1')) {
    rk <- discretisedRank(diagram, xcoords, ycoords, name)
    rks <- append(rks, rk)
  }
  return(rks)
}

#' displayRanks
#' 
#' Use a heatmap to display a the rank functions of a diagram.
#' 
#' @param diagram A persistence diagam
#' @param xcoords A list of x-coordinates
#' @param ycoords A list of y-coordinates
#' @return displays a heat map of the rank functions
#' 
#' @export 
#' 

displayRanks <- function(ranks, xcoords, ycoords) {
  rk_Ord0 <- ranks[[1]]
  rk_Ext0 <- ranks[[2]]
  rk_Rel1 <- ranks[[3]]
  rk_Ext1 <- ranks[[4]]

  heatmap_Ord0 <- ComplexHeatmap::Heatmap(rk_Ord0, name = "Ord0")
  heatmap_Ext0 <- ComplexHeatmap::Heatmap(rk_Ext0, name = "Ext0")
  heatmap_Rel1 <- ComplexHeatmap::Heatmap(rk_REl1, name = "Rel1")
  heatmap_Ext1 <- ComplexHeatmap::Heatmap(rk_Ext1, name = "Ext1")

  heatmap_Ord0 + heatmap_Ext0 + heatmap_Rel1 + heatmap_Ext1
}
