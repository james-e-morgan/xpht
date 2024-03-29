#' Centre and Scale Extended Persistence Diagrams
#'
#' The function [centreScaleDiagrams()] translates the points in an extended
#'  persistence diagram to be centred about the origin and scales the points
#'  to make different diagrams comparable.
#'
#' Given a collection of persistence diagrams in \eqn{N} evenly spaced
#'  directions, shift the points in each diagram so that they correspond to the
#'  points of the persistence diagrams of the centred image. This is equivalent
#'  to centring the convex hull of the original image. Then scale the
#'  birth and death times of each point by some factor.
#'
#' @param diagrams A list of extended persistence diagrams on a single image
#' computed in \eqn{N} evenly spaced directions.
#' @param scale Flag to indicate whether diagrams are to be scaled or not. The
#'  default value is `TRUE`.
#' @param scaleConstant `numeric` constant to control the scaling of the
#'  diagram. This is *not* the percentage of the original scale of the image.
#'  The default value is 1.
#' @return A list of extended persistence diagrams centred at the origin and
#' with birth and death times scaled by an amount proportional to the sum of the
#' minimum birth times.
#' @export
centreScaleDiagrams <- function(diagrams,
                                scale = TRUE,
                                scaleConstant = 1) {
  n_dirs <- length(diagrams)
  lambda <- findMinBirthTimes(
    n_dirs = n_dirs,
    diagrams = diagrams
  )

  adjusted_diagrams <- centreDiagrams(
    diagrams = diagrams,
    n_dirs = n_dirs,
    lambda = lambda
  )
  cat("Diagrams successfully centred.")
  if (scale) {
    adjusted_diagrams <- scaleDiagrams(
      diagrams = adjusted_diagrams,
      n_dirs = n_dirs,
      lambda = lambda,
      scaleConstant = scaleConstant
    )
    cat("Diagrams successfully scaled.")
  }

  return(adjusted_diagrams)
}

findMinBirthTimes <- function(n_dirs, diagrams) {
  lambda <- vector()
  # The first component born in any direction belongs to Ess0
  for (i in 1:n_dirs) {
    birth_times <- c(diagrams[[i]][["Ess0"]][, 1])
    lambda <- append(lambda, min(birth_times))
  }
  return(lambda)
}

centreDiagrams <- function(diagrams, n_dirs, lambda) {
  cp <- findCentre(
    n_dirs = n_dirs,
    lambda = lambda
  )

  directions <- t(sapply(1:n_dirs, function(r) {
    c(
      cos(2 * pi * r / n_dirs),
      sin(2 * pi * r / n_dirs)
    )
  }))

  for (i in 1:n_dirs) {
    shift <- sum(cp * directions[i, ])

    diagrams[[i]][["Ess0"]] <- diagrams[[i]][["Ess0"]] - shift

    if (length(diagrams[[i]][["Ord0"]]) > 0) {
      diagrams[[i]][["Ord0"]] <- diagrams[[i]][["Ord0"]] - shift
    }

    if (length(diagrams[[i]][["Rel1"]]) > 0) {
      diagrams[[i]][["Rel1"]] <- diagrams[[i]][["Rel1"]] - shift
    }

    if (length(diagrams[[i]][["Ess1"]]) > 0) {
      diagrams[[i]][["Ess1"]] <- diagrams[[i]][["Ess1"]] - shift
    }
  }
  return(diagrams)
}

findCentre <- function(n_dirs, lambda) {
  directions <- t(sapply(1:n_dirs, function(r) {
    c(
      cos(2 * pi * r / n_dirs),
      sin(2 * pi * r / n_dirs)
    )
  }))
  K <- sum(directions[,1]^2)
  directions <- lambda * directions
  cp <- colSums(directions) / K

  return(cp)
}

scaleDiagrams <- function(diagrams, n_dirs, lambda, scaleConstant) {
  scale_denom <- (-1) * sum(lambda)

  if (scale_denom <= 0) {
    stop("Incorrect lambda's computed. Cannot scale by a negative value.")
  }

  scale_value <- scaleConstant / scale_denom

  for (i in 1:n_dirs) {
    diagrams[[i]][["Ess0"]] <- diagrams[[i]][["Ess0"]] * scale_value

    if (length(diagrams[[i]][["Ord0"]]) > 0) {
      diagrams[[i]][["Ord0"]] <- diagrams[[i]][["Ord0"]] * scale_value
    }

    if (length(diagrams[[i]][["Rel1"]]) > 0) {
      diagrams[[i]][["Rel1"]] <- diagrams[[i]][["Rel1"]] * scale_value
    }

    if (length(diagrams[[i]][["Ess1"]]) > 0) {
      diagrams[[i]][["Ess1"]] <- diagrams[[i]][["Ess1"]] * scale_value
    }
  }
  return(diagrams)
}
