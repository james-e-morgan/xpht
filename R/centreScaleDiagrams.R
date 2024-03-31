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
#' @inheritParams extractBoundary
#' @return A list of extended persistence diagrams centred at the origin and
#' with birth and death times scaled by an amount proportional to the sum of the
#' minimum birth times.
#' @export
centreScaleDiagrams <- function(diagrams,
                                scale = TRUE,
                                scaleConstant = 1,
                                saveOutput = FALSE,
                                outputDir = NULL,
                                fName = NULL,
                                verbose = TRUE) {
  if (saveOutput) {
    if (!dir.exists(outputDir)) {
      outputDir <- getwd()
      cat("Output directory doesn't exist. Saving to working directory.\n")
    }
    if (is.null(fName)) {
      fName <- readline(prompt = "Please provide a filename for save: ")
    }
    out_file <- paste(outputDir, "/", fName, ".RDS", sep = "")
  }

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

  if (scale) {
    adjusted_diagrams <- scaleDiagrams(
      diagrams = adjusted_diagrams,
      n_dirs = n_dirs,
      lambda = lambda,
      scaleConstant = scaleConstant
    )
  }

  if (saveOutput) {
    saveRDS(adjusted_diagrams, file = out_file)
    if (verbose) {
      cat("Successfully saved ", out_file, "\n", sep = "")
    }
  } else {
    if (verbose) {
      cat("Diagrams successfully centred.")
      if (scale) {
        cat("Diagrams successfully scaled.")
      }
    }
    return(adjusted_diagrams)
  }
}

#' Centre and Scale Extended Persistence Diagrams for multiple XPHTs
#'
#' Runs [centreScaleDiagrams()] on all .RDS files in a given directory.
#'
#' @inheritParams centreScaleDiagrams
#' @inheritParams multiExtractBoundary
#' @return If `saveOutput = TRUE`, then the output of [centreScaleDiagrams()]
#'  for each object will be saved in a separate file. If `saveOutput = FALSE`,
#'  then multiCentreScaleDiagrams()] returns a list containing the extended
#'  persistence diagrames centred at the origin and with birth and death times
#'  scaled by an amount proportional to the sum of the minimum birth times for
#'  each object in `inputDir`. The \eqn{i}-th entry in the list is the centred
#'  extended persistence diagrams for the \eqn{i}-th file in `inputDir`.
#' @seealso centreScaleDiagrams()
#' #' @importFrom utils tail
#' @export
multiCentreScaleDiagrams <- function(
  inputDir,
  scale = TRUE,
  scaleConstant = 1,
  saveOutput = FALSE,
  outputDir = NULL,
  verbose = TRUE
) {
  files <- list.files(
    path = inputDir,
    pattern = "*.RDS",
    full.names = TRUE,
    recursive = FALSE
  )

  if (saveOutput) {
    if (!dir.exists(outputDir)) {
      outputDir <- getwd()
      cat("Output directory doesn't exist. Saving to working directory.\n")
    }
  } else {
    centred_diagrams <- vector(mode = "list", length = length(files))
  }

  for (i in seq_along(files)) {
    f <- files[[i]]
    if (verbose) {
      cat("Commencing,", f, "\n", sep = " ")
    }

    fName <- tail(strsplit(strsplit(f, ".", fixed = TRUE)[[1]][1], "/",
                    fixed = TRUE
                  )[[1]], n = 1)

    diagrams <- readRDS(f)

    if (saveOutput) {
      centreScaleDiagrams(
        diagrams = diagrams,
        scale = scale,
        scaleConstant = scaleConstant,
        saveOutput = TRUE,
        outputDir = outputDir,
        fName = fName,
        verbose = verbose
      )
    } else {
      centred_diagrams[[i]] <- centreScaleDiagrams(
        diagrams = diagrams,
        scale = scale,
        scaleConstant = scaleConstant,
        saveOutput = FALSE,
        verbose = verbose
      )
    }
  }

  if (saveOutput) {
    cat("All centred extended persistence diagrams saved in:\n",
        outputDir, "\n", sep = "")
  } else {
    cat("All centred extended persistence diagrams computed.\n")
    return(centred_diagrams)
  }
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
  K <- sum(directions[, 1]^2)
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
