#' Centre and Scale the Extended Persistence Diagrams of an image.
#'
#' Given a collection of persistence diagrams in $N$ evenly spaced directions,
#' shift the points in each diagram so that the correspond to the points of the
#' persistence diagrams of the centred image. Then scale the birth and death
#' times of each point by some factor.
#'
#' Centring and scaling the persistence diagrams allows for a more accurate
#' comparison of their shapes.
#'
#' @param diagrams A list of extended persistence diagrams on a single image
#' computed in N evenly spaced directions.
#' @param scaleConstant A constant to control the scaling of the diagram.
#' Default is 1. Note this is not the percentage of the original scale of the
#' image.
#' @return A list of extended persistence diagrams of the same image, centred
#' at the origin and scaled by an amount proportional to the sum of the
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
  print("Diagrams successfully centred.")
  if (scale) {
    adjusted_diagrams <- scaleDiagrams(
      diagrams = adjusted_diagrams,
      n_dirs = n_dirs,
      lambda = lambda,
      scaleConstant = scaleConstant
    )
    print("Diagrams successfully scaled.")
  }

  return(adjusted_diagrams)
}

findMinBirthTimes <- function(n_dirs, diagrams) {
  lambda <- vector()
  # The first component born in any direction belongs to Ext0
  for (i in 1:n_dirs) {
    birth_times <- c(diagrams[[i]][["Ext0"]][, 1])
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
    shift <- sum(cp, directions[i, ])

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

findCentre <- function(n_dirs, lambda) {
  directions <- t(sapply(1:n_dirs, function(r) {
    c(
      cos(2 * pi * r / n_dirs),
      sin(2 * pi * r / n_dirs)
    )
  }))
  directions <- lambda * directions
  cp <- colSums(directions) / (n_dirs * 0.5)

  return(cp)
}

scaleDiagrams <- function(diagrams, n_dirs, lambda, scaleConstant) {
  scale_denom <- (-1) * sum(lambda)

  if (scale_denom <= 0) {
    stop("Incorrect lambda's computed. Cannot scale by a negative value.")
  }

  scale_value <- scaleConstant / scale_denom

  for (i in 1:n_dirs) {
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
