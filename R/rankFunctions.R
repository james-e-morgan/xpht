###################################
######### RANK FUNCTIONS ##########
###################################
# Author: Yossi Bokor
# Version: 1
# Date: February 17, 2022

#' Obtain the discretised rank function of a persistence diagram
#' 
#' The rank function of a diagram at a point (x,y) is the number of points above
#' and to the left of this point. We discretise it by evaluating on a grid of points.
#' The grid is obtained by taking the cartesian product of a list of x-coordinates and 
#' y-coordinates. 
#' 
#' OUTPUT
#' ------
#' The values of 
#'      - boundary[[1]] is a vector (m,n) containing the number of anticlockwise and clockwise curves.
#'      - boundary[[2]] to boundary[[m+1]] contain the anticlockwise curves.
#'      - boundary[[m+2]] to boundary[[n+1]] contain the clockwise curves.
#' If set to save, the list is saved as an RDS file.
#' 
#' @param diagram A n x 2 array of points in the persistence diagram
#' @param xcoords The x-coordinates of the grid you want to evaluate the rank function on, in increasing order.
#' @param ycoords The y-coordinates of the grid you want to evaluate the rank function on, in increasing order.
 
#'
#' @param diagram A single persistence diagram.
#' @return

#' @export

discretisedRank <- function(diagram, xcoords,ycoords){
  nx = length(xcoords)
  ny = length(ycoords)
  np = length(diagram[,1])
  rk = matrix(0,nrow=ny, ncol=nx)

  for (i in 1:np){
    xi = Position(function(x) x >= diagram[i,1], xcoords)
    yi = Position(function(y) y >= diagram[i,2], ycoords)
    xinc = c(xi:nx)
    yinc = c(1:yi)
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