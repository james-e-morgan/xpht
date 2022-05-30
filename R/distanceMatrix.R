#' @export
computeDistanceMatrix <- function(diagrams,
                                  n_objects,
                                  q = 1,
                                  aligned = FALSE) {
  n_dirs <- length(diagrams) / n_objects

  distance_matrix <- matrix(0, nrow = n_objects, ncol = n_objects)

  # Distance matrix is symmetric, so only compute upper triangular part.
  # Diagonal entries are zero.
  for (i in 1:(n_objects - 1)) {
    idx_1 <- ((i - 1) * n_dirs + 1):(i * n_dirs)
    for (j in (i + 1):(n_objects)) {
      idx_2 <- ((j - 1) * n_dirs + 1):(j * n_dirs)

      object_1 <- lapply(idx_1, function(k) diagrams[[k]])
      object_2 <- lapply(idx_2, function(k) diagrams[[k]])

      if (aligned) {
        d_ij <- computeAlignedDistance(
          object_1 = object_1,
          object_2 = object_2,
          q = q,
          n_dirs = n_dirs
        )
      } else {
        d_ij <- computeUnalignedDistance(
          object_1 = object_1,
          object_2 = object_2,
          q = q,
          n_dirs = n_dirs
        )
      }

      d_ij <- (d_ij / n_dirs)^(1 / q)

      distance_matrix[i, j] <- d_ij
      distance_matrix[j, i] <- d_ij

      print(paste("Computed distance (", i, ",", j, ") and (",
        j, ",", i, ")",
        sep = ""
      ))
    }
  }

  return(distance_matrix)
}

computeUnalignedDistance <- function(object_1, object_2, q, n_dirs) {
  ds <- c("Ext0", "Ext1", "Ord0", "Rel1")
  # Aligned distance gives a starting point
  min_dist <- computeAlignedDistance(
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
          computePointDistance(
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

computeAlignedDistance <- function(object_1, object_2, q, n_dirs) {
  ds <- c("Ord0", "Rel1", "Ext0", "Ext1")
  total <- 0

  for (x in ds) {
    d <- sum(sapply(
      1:n_dirs,
      function(k) {
        computePointDistance(
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

computePointDistance <- function(points_1, points_2, q) {
  n_1 <- length(points_1) / 2
  n_2 <- length(points_2) / 2

  if (n_1 + n_2 == 0) {
    # No points
    return(0)
  } else if (n_1 > 0 && n_2 == 0) {
    # Only points_1 has points
    d <- sum(sapply(
      1:n_1,
      function(i) computeDistanceToDiagonal(points_1[i, ], q)
    ))

    return(d)
  } else if (n_1 == 0 && n_2 > 0) {
    # Only points_2 has points
    d <- sum(sapply(
      1:n_2,
      function(i) computeDistanceToDiagonal(points_2[i, ], q)
    ))

    return(d)
  } else {
    # Hungarian Algorithm
    cost_matrix <- vector()

    for (i in 1:n_1) {
      r_1 <- apply(points_2, 1, function(x) compute_lp_distance(points_1[i, ], x, q))

      dist_xy <- computeDistanceToDiagonal(points_1[i, ], q)
      r_2 <- rep(dist_xy, n_1)

      cost_matrix <- rbind(cost_matrix, c(r_1, r_2))

      r_1 <- apply(points_2, 1, function(x) computeDistanceToDiagonal(x, q))
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
}

computeDistanceToDiagonal <- function(p, q) {
  return(2 * ((abs(p[2] - p[1]) / 2)^q))
}

compute_lp_distance <- function(p1, p2, q) {
  return(abs(p1[1] - p2[1])^q + abs(p1[2] - p2[2])^q)
}
