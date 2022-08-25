#' Plot an Extended Persistence Diagram or Barcode
#'
#' The function [plot.extendedPHT()] plots the extended persistence diagram
#'  or extended persistence barcode of an `extendedPHT` object.
#'
#' This is a generic function. Plot the the extended persistence diagram or the
#'  extended persistence barcode of an object. If the diagram is plotted, then
#'  the diagonal \eqn{birth = death} is also plotted. The classes in Ord0 and
#'  Ess0 have \eqn{birth < death}  and are plotted above the diagonal, whereas
#'  classes in Rel1 and Ess1 have \eqn{birth > death} and are plotted below the
#'  diagonal. If the barcode is plotted then classes will be stacked in the
#'  order Ord0, Rel1, Ess0, Ess1.
#'
#' @param x `extendedPHT` object to plot.
#' @param barcode If `TRUE`, then barcode is plotted. Default value is `FALSE`.
#' @param col `vector` containing the four colours to use in plotting in the
#'  order Ord0, Rel1, Ess0, Ess1. The default value is
#'    * "#000000" for Ord0
#'    * "#e69f00" for Rel1
#'    * "#009e73" for Ess0
#'    * "#d55e00" for Ess1
#' @param pch Symbols to use if plotting extended persistence diagram in order
#'  Ord0, Rel1, Ess0, Ess1. The default value is `15:18`.
#' @param lab.line Control distance of axis labels to axes. The default value is
#'  2.2.
#' @inheritParams graphics::plot
#' @export
plot.extendedPHT <- function(x,
                            barcode = FALSE,
                            col = c( "#000000",
                                     "#e69f00",
                                     "#009e73",
                                     "#d55e00"),
                            pch = 15:18,
                            lwd = 2,
                            cex = 1,
                            lab.line = 2.2, ...) {
  # Add checks for diagram list
  # Add checks for diagram Names
  X <- vector()

  if (length(col) != 4) {
    stop("Must provide four colours for plotting.")
  }

  for (i in 2:5) {
    if (length(x[[i]]) > 2) {
      sortedx <- x[[i]][order(x[[i]][, 1], decreasing = FALSE), ]
      sortedx <- cbind(rep(i - 1, length(sortedx) / 2), sortedx)
      X <- rbind(X, sortedx)
    } else if (length(x[[i]]) == 2) {
      sortedx <- cbind(i - 1, x[[i]])
      X <- rbind(X, sortedx)
    }
  }

  if (barcode) {
    temp_col <- col
    col <- c(
      rep(temp_col[1], length(which(X[, 1] == 1))),
      rep(temp_col[2], length(which(X[, 1] == 2))),
      rep(temp_col[3], length(which(X[, 1] == 3))),
      rep(temp_col[4], length(which(X[, 1] == 4)))
    )

    birth <- X[, 2]
    death <- X[, 3]

    first <- min(birth)
    last <- max(death)

    n <- length(birth)

    graphics::plot(c(first, last),
      c(1, n + 1),
      type = "n",
      xlab = "",
      ylab = "",
      xlim = c(first, last),
      ylim = c(0, n + 1),
      xaxt = "n",
      yaxt = "n", ...
    )

    graphics::axis(1)

    graphics::title(
      xlab = "height",
      line = lab.line
    )

    graphics::segments(birth, 1:n,
      death, 1:n,
      lwd = rep(lwd, n),
      lty = rep(1, n),
      col = col
    )

  } else {

    first <- min(X[, 2])
    last <- max(X[, 3])

    # initialise empty plot
    graphics::plot(0, 0,
      type = "n",
      xlab = "",
      ylab = "",
      xlim = c(first, last),
      ylim = c(first, last),
      ...
    )

    for (i in 1:4) {
      subDgm <- X[X[, 1] == i, ]
      if (length(subDgm) == 3) {
        subDgm <- matrix(subDgm, nrow = 1, ncol = 3)
      }

      if (length(subDgm) > 0) {
        birth <- subDgm[, 2]
        death <- subDgm[, 3]

        graphics::points(birth,
          death,
          pch = pch[i],
          lwd = 2,
          cex = cex,
          col = col[i]
        )
      }
    }

    graphics::abline(0, 1,
      lty = 2,
      col = "#A9A9A9"
    )

    graphics::title(
      xlab = "birth",
      ylab = "death",
      line = lab.line
    )
  }
}

#' Get the Default Colours used in Plotting Extended Persistence Diagrams
#' 
#' @return Vector containing the default colours used for plotting extended
#' persistence diagrams in the order (Ord0, Rel1, Ess0, Ess1).
#' @export
getDefaultColours <- function() {
  return(c( "#000000", "#e69f00", "#009e73", "#d55e00"))
}