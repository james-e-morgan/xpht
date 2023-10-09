#' Compute the distance matrix between Extended Persistent Homology Transforms
#'
#' The function [computeDistanceMatrix()] computes the q-Wasserstein distance
#'  between the extended persistent homology transforms of multiple images.
#'
#' Given a collection of objects \eqn{O_1,\dots,O_n} and directions
#'  \eqn{v_1,\dots,v_K}, compute the distance matrix \eqn{D}$ between all pairs
#'  of objects \eqn{O_i,O_j} for \eqn{i\neq j}. The objects \eqn{O_i} and 
#'  \eqn{O_j} have extended persistence diagrams \eqn{X_{i1},\dots,X_{iK}} and
#'  \eqn{X_{j1},\dots,X_{jK}}, respectively. The distance between $O_i$ and
#'  $O_j$ is the average distance between their diagrams.
#'
#' @param diagrams `list` of extended persistence diagrams of multiple objects.
#'  Each object must have their extended persistence diagrams computed in the
#'  same number of directions. If \eqn{K} directions were used, then the first
#'  \eqn{K} entries must correspond to the first object, the next \eqn{K}
#'  entries to the next object and so on.
#' @param nObjects `int` giving the total number of objects. Used to determine
#'   the number of directions.
#' @param q `numeric` determining the \eqn{L^q} distance to use. Must have
#'  \eqn{1\leq q < \infty}. The default value is 1.
#' @param aligned `bool` to indicate whether the original images have the same
#'  alignment. The default value is `FALSE`.
#' @param verbose `bool` If `TRUE` will print when the distance between each
#' pair of objects has been computed.
#' @return A square, symmetric `matrix` where entry \eqn{(i,j)} contains the
#' \eqn{q}-Wasserstein distance between the XPHTs of objects i and j. 
#' @export
computeDistanceMatrix <- function(diagrams,
                                  nObjects,
                                  q = 1,
                                  aligned = FALSE,
                                  verbose = TRUE) {
  n_dirs <- length(diagrams) / nObjects

  distance_matrix <- matrix(0, nrow = nObjects, ncol = nObjects)

  # Distance matrix is symmetric, so only compute upper triangular part.
  # Diagonal entries are zero.
  for (i in 1:(nObjects - 1)) {
    idx_1 <- ((i - 1) * n_dirs + 1):(i * n_dirs)
    for (j in (i + 1):(nObjects)) {
      idx_2 <- ((j - 1) * n_dirs + 1):(j * n_dirs)

      object_1 <- lapply(idx_1, function(k) diagrams[[k]])
      object_2 <- lapply(idx_2, function(k) diagrams[[k]])

      if (aligned) {
        d_ij <- alignedDistance(
          object_1 = object_1,
          object_2 = object_2,
          q = q,
          n_dirs = n_dirs
        )
      } else {
        d_ij <- unalignedDistance(
          object_1 = object_1,
          object_2 = object_2,
          q = q,
          n_dirs = n_dirs
        )
      }

      d_ij <- (d_ij / n_dirs)^(1 / q)

      distance_matrix[i, j] <- d_ij
      distance_matrix[j, i] <- d_ij
      
      if (verbose) {
        cat("Computed distance (", i, ",", j, ") and (",
          j, ",", i, ")",
          sep = ""
        )
      }
    }
  }
  cat("Successfully computed all distances.")
  return(distance_matrix)
}

#' Stack Multiple XPHTs into a Single List
#' 
#' Given a directory of XPHTs, load each XPHT and store as an element in a
#' list.
#' 
#' @param inputDir The directory containing the XPHTs stored as `.RDS` files.
#' @return A list of XPHTs
#' @export
stackDiagrams <- function(inputDir) {
  files = list.files(inputDir, full.names=TRUE, recursive=FALSE)
  data = vector(mode = "list")
  cx = 1
  for (i in seq_along(files)) {
    dgm = readRDS(files[[i]])
    for (j in seq_along(dgm)) {
      data[[cx]] = dgm[[j]]
      cx = cx + 1
    }
  }
  cat('Successfully loaded XPHTs of', length(files), 'images.\n', sep=" ")
  return(data)
}

unalignedDistance <- function(object_1, object_2, q, n_dirs) {
  ds <- c("Ess0", "Ess1", "Ord0", "Rel1")
  # Aligned distance gives a starting point
  min_dist <- alignedDistance(
    object_1 = object_1,
    object_2 = object_2,
    q = q,
    n_dirs = n_dirs
  )
  if (min_dist == 0) {
    return(min_dist)
  }
  # Only compute distances until min_dist is exceeded
  for (m in 0:(n_dirs - 1)) {
    total <- 0
    for (x in ds) {
      d <- sum(sapply(
        1:n_dirs,
        function(k) {
          pointDistance(
            points_1 = object_1[[k]][[x]],
            points_2 = object_2[[((k + m) %% n_dirs) + 1]][[x]],
            q = q
          )
        }
      ))

      total <- total + d
      if (total >= min_dist) {
        break
      }
    }
    if (total < min_dist) {
      min_dist <- total
    }
  }
  return(min_dist)
}

alignedDistance <- function(object_1, object_2, q, n_dirs) {
  ds <- c("Ord0", "Rel1", "Ess0", "Ess1")
  total <- 0

  for (x in ds) {
    d <- sum(sapply(
      1:n_dirs,
      function(k) {
        pointDistance(
          points_1 = object_1[[k]][[x]],
          points_2 = object_2[[k]][[x]],
          q = q
        )
      }
    ))

    total <- total + d
  }

  return(total)
}

pointDistance <- function(points_1, points_2, q) {
  n_1 <- length(points_1) / 2
  n_2 <- length(points_2) / 2

  if (n_1 + n_2 == 0) {
    # No points
    return(0)
  } else if (n_1 > 0 && n_2 == 0) {
    # Only points_1 has points
    d <- sum(sapply(
      1:n_1,
      function(i) distanceToDiagonal(points_1[i, ], q)
    ))

    return(d)
  } else if (n_1 == 0 && n_2 > 0) {
    # Only points_2 has points
    d <- sum(sapply(
      1:n_2,
      function(i) distanceToDiagonal(points_2[i, ], q)
    ))

    return(d)
  } else {
    # Hungarian Algorithm
    cost_matrix <- vector()

    for (i in 1:n_1) {
      r_1 <- apply(points_2, 1, function(x) distanceLq(points_1[i, ], x, q))

      dist_xy <- distanceToDiagonal(points_1[i, ], q)
      r_2 <- rep(dist_xy, n_1)

      cost_matrix <- rbind(cost_matrix, c(r_1, r_2))
    }

      r_1 <- apply(points_2, 1, function(x) distanceToDiagonal(x, q))
      r_2 <- rep(0, n_1)

      for (j in 1:n_2) {
        cost_matrix <- rbind(cost_matrix, c(r_1, r_2))
      }
      
      pairing <- RcppHungarian::HungarianSolver(cost_matrix)
      idxs <- pairing[["pairs"]]

      pair_vals <- apply(
        idxs, 1,
        function(x) cost_matrix[x[1], x[2]]
      )

      return(sum(pair_vals))
  }
}

distanceToDiagonal <- function(p, q) {
  return(2 * ((abs(p[2] - p[1]) / 2)^q))
}

distanceLq <- function(p1, p2, q) {
  return(abs(p1[1] - p2[1])^q + abs(p1[2] - p2[2])^q)
}
