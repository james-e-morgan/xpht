#' Compute the Extended Persistence Diagrams of a given image
#'
#' Given a disjoint collection of simple closed curves representing the
#' boundaries of a two dimensional shape, compute the zeroth order extended
#' persistent homology transform of the shape in a given number of directions.
#'
#' @param bdry_curves A list containing points along the boundary curves of the
#' image. The list must have length number_of_curves + 1 with the first entry a
#' vector whose entries are the numbers of positively oriented curves and
#' negatively oriented curves, respectively. All of the positively oriented
#' curves must be listed before the negatively oriented curves.
#' @param img_name A string giving the name of the image.
#' @param n_dirs An integer giving the number of directions to use. This
#' must be positive and even.
#' @param tolerance Parameter to control the noise from the image. Any vertices
#' closer than this in a particular direction will be classified as noise and
#' ignored. (Default: 1/sqrt(2))
#' @param save_output If TRUE, will save output to directory specified by
#' output_dir. (Default: FALSE)
#' @param output_dir The directory to save the output. If save_output is TRUE and
#' no directory is specified, saves to working directory. (Default: NULL)
#' @param f_name The name of the output file saved. If save_output is TRUE and no
#' filename specified, prompts user for filename. (Default: NULL)
#' @param verbose If TRUE, prints indictors of progress throughout.
#' (Default: TRUE)
#' @return A list of length 5 with names "Name", "Ord0", "Rel1", "Ext0",
#' "Ext1". Each entry is a list
#' containing the relevant part of the extended persistence diagram. The i-th
#' entry corresponds to the i-th direction.
#' @export
computeExtendedPersistence <- function(bdry_curves,
                                       img_name,
                                       n_dirs,
                                       tolerance = 1 / sqrt(2),
                                       save_output = FALSE,
                                       output_dir = NULL,
                                       f_name = NULL,
                                       verbose = TRUE) {
  if (save_output) {
    if (!dir.exists(output_dir)) {
      output_dir <- getwd()
      print("Output directory doesn't exist.
                  Saving to working directory.")
    }
    if (is.null(f_name)) {
      f_name <- readLine(prompt = "Please provide a filename for save: ")
    }
    out_file <- paste(output_dir, "/", f_name, ".RDS", sep = "")
  }
  if (n_dirs < 2 || n_dirs %% 2 == 1) {
    stop("Number of directions must be at least two and even")
  } else {
    directions <- t(sapply(1:n_dirs, function(r) {
      c(
        cos(2 * pi * r / n_dirs),
        sin(2 * pi * r / n_dirs)
      )
    }))
    midpoint <- (n_dirs / 2) + 1
  }

  x_diagram <- vector(mode = "list")

  skeleton <- parseSkeleton(bdry_curves)

  ac_curves <- bdry_curves[[1]][1]
  c_curves <- bdry_curves[[1]][2]

  for (d in 1:(n_dirs / 2)) {
    diagram_name <- paste(img_name, "-", toString(d), sep = "")

    dir_vector <- directions[d, ]

    n_components <- length(skeleton)
    zeroth_diagram <- vector(mode = "list", length = n_components)

    count <- 1

    for (i in 1:n_components) {
      skeleton_component <- skeleton[i][[1]]

      filtration <- computeHeightFiltration(
        curve = skeleton_component,
        direction = dir_vector
      )

      zeroth_diagram[[i]] <- computeDiagram(
        filtration = filtration,
        direction = dir_vector,
        tolerance = tolerance
      )
    }
    extendedDiagram <- computeExtendedDiagram(
      zerothDiagram = zeroth_diagram,
      nComponents = n_components,
      diagramName = diagramName
    )
    class(x_diagram) <- "extended_diagram"
    x_diagram[[d]] <- x_diagram
  }

  for (d in midpoint:n_dirs) {
    diagram_name <- paste(img_name, "-", toString(d), sep = "")
    neg_diagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
    negDiagram <- sapply(neg_diagram.names, function(x) NULL)

    neg_diagram[["Name"]] <- diagram_name
    k <- d - midpoint + 1

    neg_diagram[["Ord0"]] <- -x_diagram[[k]][["Rel1"]]
    neg_diagram[["Rel1"]] <- -x_diagram[[k]][["Ord0"]]
    neg_diagram[["Ext0"]] <- -matrix(x_diagram[[k]][["Ext0"]][, c(2, 1)],
      ncol = 2
    )

    if (length(x_diagram[[k]][["Ext1"]]) > 0) {
      neg_diagram[["Ext1"]] <- -matrix(x_diagram[[k]][["Ext1"]][, c(2, 1)],
        ncol = 2
      )
    } else {
      neg_diagram[["Ext1"]] <- vector()
    }
    class(neg_diagram) <- "extended_diagram"
    x_diagram[[d]] <- neg_diagram
  }

  class(x_diagram) <- "extended_diagram"

  if (save_output) {
    saveRDS(x_diagram, file = out_file)
    if (verbose) {
      cat("Successfully saved ", out_file, "\n", sep = "")
    }
  } else {
    if (verbose) {
      cat("Extended persistence diagrams successfully computed for",
        n_dirs, "directions.\n",
        sep = " "
      )
    }
    return(x_diagram)
  }
}

#' Compute the Extended PHT for multiple images
#'
#' Runs \code{\link{extractBoundary}} on all .RDS files in a given directory.
#'
#' Each persistence diagram is \emph{saved} in in the specified
#' output directory.
#'
#' For more information, see \code{\link{extractBoundary}}.
#'
#' @param input_dir The directory containing the extracted boundary curves.
#' @param output_dir The directory to save the output. If save_output is TRUE
#' and no directory is specified, saves to working directory. (Default: NULL)
#' @param n_dirs An integer giving the number of directions to use. This
#' must be positive and even.
#' @param tolerance Parameter to control the noise from the image. Any vertices
#' closer than this in a particular direction will be classified as noise and
#' ignored. (Default: 1/sqrt(2))
#' @param verbose If TRUE, prints indictors of progress throughout.
#' (Default: TRUE)
#' @export
multiExtendedPersistence <- function(input_dir,
                                     output_dir,
                                     n_dirs,
                                     tolerance = 1 / sqrt(2),
                                     verbose = TRUE) {
  if (!dir.exists(output_dir)) {
    output_dir <- getwd()
    print("Output directory doesn't exist. Saving to working directory.")
  }

  files <- list.files(
    path = input_dir,
    pattern = "*.RDS",
    full.names = TRUE,
    recursive = FALSE
  )

  for (i in seq_along(files)) {
    f <- files[[i]]
    if (verbose) {
      cat("Commencing,", f, "\n", sep = " ")
    }

    f_name <- tail(strsplit(strsplit(f, ".", fixed = TRUE)[[1]][1],
      "/",
      fixed = TRUE
    )[[1]], n = 1)

    bdry_curves <- readRDS(f)

    extendedPersistence(
      bdry_curves = bdry_curves,
      img_name = f_name,
      n_dirs = n_dirs,
      tolerance = tolerance,
      save_output = TRUE,
      output_dir = output_dir,
      f_name = f_name,
      verbose = verbose
    )
  }
}

parseSkeleton <- function(bdry_curves) {
  complex <- vector(mode = "list")

  skeleton.names <- c("vertex", "edge", "coords")

  ac_curves <- bdry_curves[[1]][1]
  c_curves <- bdry_curves[[1]][2]
  n_curves <- ac_curves + c_curves

  for (i in 2:(n_curves + 1)) {
    curve <- bdry_curves[[i]]
    skeleton <- sapply(skeleton.names, function(x) NULL)

    if (all(curve[1, 1:2] == curve[nrow(curve), 1:2])) {
      np <- nrow(curve) - 1
    } else {
      np <- nrow(curve)
    }

    skeleton[["vertex"]] <- 1:np
    skeleton[["edge"]] <- vector()
    skeleton[["coords"]] <- vector()

    for (j in 1:(np - 1)) {
      skeleton[["coords"]] <- rbind(
        skeleton[["coords"]],
        curve[j, 1:2]
      )
      skeleton[["edge"]] <- rbind(
        skeleton[["edge"]],
        c(j, j + 1)
      )
    }
    skeleton[["edge"]] <- rbind(
      skeleton[["edge"]],
      c(np, 1)
    )
    skeleton[["coords"]] <- rbind(
      skeleton[["coords"]],
      curve[np, 1:2]
    )

    complex[[i - 1]] <- skeleton
  }

  return(complex)
}

computeHeightFiltration <- function(curve, direction) {
  filtration.names <- c("height", "lowerNbrs", "coords")
  filtration <- sapply(filtration.names, function(x) NULL)

  filtration[["coords"]] <- curve[["coords"]]

  filtration[["height"]] <- apply(
    curve[["coords"]], 1,
    function(v) dotProduct(v, direction)
  )

  n_vertices <- max(curve[["vertex"]])

  filtration[["lowerNbrs"]] <- vector(mode = "list", length = n_vertices)

  for (i in 1:n_vertices) {
    e <- curve[["edge"]][i, ]

    h_1 <- filtration[["height"]][e[1]]
    h_2 <- filtration[["height"]][e[2]]

    if (h_1 < h_2) {
      filtration[["lowerNbrs"]][[e[2]]] <- append(
        filtration[["lowerNbrs"]][[e[2]]],
        e[1]
      )
    } else if (h_1 == h_2) {
      # If two adjacent vertices have the same height, then the one with the
      # lower index is the lower neighbour
      if (e[1] < e[2]) {
        filtration[["lowerNbrs"]][[e[2]]] <- append(
          filtration[["lowerNbrs"]][[e[2]]],
          e[1]
        )
      } else {
        filtration[["lowerNbrs"]][[e[1]]] <- append(
          filtration[["lowerNbrs"]][[e[1]]],
          e[2]
        )
      }
    } else {
      filtration[["lowerNbrs"]][[e[1]]] <- append(
        filtration[["lowerNbrs"]][[e[1]]],
        e[2]
      )
    }
  }
  return(filtration)
}

dotProduct <- function(v_1, v_2) {
  return(sum(v_1 * v_2))
}

computeDiagram <- function(filtration, direction, tolerance) {
  sorted_heights <- unique(sort(filtration[["height"]]))
  n_vertices <- length(filtration[["height"]])

  diagram.names <- c("finite", "finite_minimal", "extended", "extended_minimal")
  diagram <- sapply(diagram.names, function(x) NULL)

  parents <- vector(
    mode = "list",
    length = length(filtration[["height"]])
  )

  diagram[["extended"]] <- c(
    head(sorted_heights, 1),
    tail(sorted_heights, 1)
  )

  for (h in sorted_heights) {
    h_vertices <- which(filtration[["height"]] == h)

    for (v in h_vertices) {
      if (is.null(filtration[["lowerNbrs"]][[v]])) {
        # The vertex v has no lower neighbours.
        parents[[v]] <- v
      } else {
        components <- sapply(
          filtration[["lowerNbrs"]][[v]],
          function(x) findParent(x, parents)
        )
        # ignore multiple paths to same parent
        components <- unique(components)

        if (length(components) == 1) {
          parents[[v]] <- components
        } else {
          birth_times <- sapply(
            components,
            function(x) filtration[["height"]][x]
          )
          min_birth_time <- min(birth_times)

          components <- sort(components)

          # If there are two candidates for the birth of a component,
          # we take the one with the lower index.
          # This is for consistency.
          count <- 0

          for (x in components) {
            h_x <- filtration[["height"]][x]
            if (h_x > min_birth_time) {
              if (h_x < h) {
                # Component born at height of x dies at
                # the current height h.
                if (abs(h_x - h) > tolerance) {
                  diagram[["finite"]] <- rbind(
                    diagram[["finite"]],
                    c(h_x, h)
                  )
                  # Check if vertex x is minimal
                  idxs <- (c(x - 2, x - 1, x) %% n_vertices) + 1
                  vertices <- c(
                    filtration[["coords"]][idxs[1], ],
                    filtration[["coords"]][idxs[2], ],
                    filtration[["coords"]][idxs[3], ]
                  )
                  heights <- c(
                    filtration[["height"]][idxs[1]],
                    filtration[["height"]][idxs[2]],
                    filtration[["height"]][idxs[3]]
                  )
                  diagram[["finite_minimal"]] <- append(
                    diagram[["finite_minimal"]],
                    testMinimality(
                      vertices = vertices,
                      heights = heights,
                      direction = direction
                    )
                  )
                }
              }
            } else if (h_x == min_birth_time) {
              if (count == 0) {
                new_component <- x
                parents[[v]] <- new_component
                count <- 1
              } else {
                if (h_x < h) {
                  if (abs(h_x - h) > tolerance) {
                    diagram[["finite"]] <- rbind(
                      diagram[["finite"]],
                      c(h_x, h)
                    )
                    # Check if vertex x is minimal
                    idxs <- (c(x - 2, x - 1, x) %% n_vertices) + 1
                    vertices <- c(
                      filtration[["coords"]][idxs[1], ],
                      filtration[["coords"]][idxs[2], ],
                      filtration[["coords"]][idxs[3], ]
                    )
                    heights <- c(
                      filtration[["height"]][idxs[1]],
                      filtration[["height"]][idxs[2]],
                      filtration[["height"]][idxs[3]]
                    )
                    diagram[["finite_minimal"]] <- append(
                      diagram[["finite_minimal"]],
                      testMinimality(
                        vertices = vertices,
                        heights = heights,
                        direction = direction
                      )
                    )
                  }
                }
              }
            }
          }
          # All components found are part of same connected component.
          # Update this information
          for (x in components) {
            parents[[x]] <- new_component
          }
        }
      }
    }
  }
  birth_point <- unique(sapply(
    unique(parents),
    function(x) findParent(x, parents)
  ))

  if (length(birth_point) > 1) {
    stop("Simple closed curve has more than one essential class.")
  }

  # Check if vertex x is minimal
  x <- birth_point[1]
  idxs <- (c(x - 2, x - 1, x) %% n_vertices) + 1
  vertices <- c(
    filtration[["coords"]][idxs[1], ],
    filtration[["coords"]][idxs[2], ],
    filtration[["coords"]][idxs[3], ]
  )
  heights <- c(
    filtration[["height"]][idxs[1]],
    filtration[["height"]][idxs[2]],
    filtration[["height"]][idxs[3]]
  )

  diagram[["extended_minimal"]] <- testMinimality(
    vertices = vertices,
    heights = heights,
    direction = direction
  )

  return(diagram)
}

findParent <- function(x,
                       parents) {
  if (parents[[x]] == x) {
    return(x)
  } else {
    px <- findParent(parents[[x]], parents)
    return(px)
  }
}

computeExtendedDiagram <- function(zeroth_diagram,
                                   n_components,
                                   diagram_name) {
  ex_diagram.names <- c("Name", "Ord0", "Rel1", "Ext0", "Ext1")
  ex_diagram <- sapply(ex_diagram.names, function(x) NULL)

  ord_0 <- vector()
  rel_1 <- vector()
  ext_0 <- vector()
  ext_1 <- vector()

  ex_diagram[["Name"]] <- diagram_name
  for (i in 1:n_components) {
    diagram <- zeroth_diagram[[i]]

    if (length(diagram[["finite"]]) > 0) {
      nFinite <- nrow(diagram[["finite"]])

      for (j in 1:nFinite) {
        if (diagram[["finite_minimal"]][j]) {
          ord_0 <- rbind(
            ord_0,
            diagram[["finite"]][j, ]
          )
        } else {
          rel_1 <- rbind(
            rel_1,
            rev(diagram[["finite"]][j, ])
          )
        }
      }
    }

    if (diagram[["extended_minimal"]]) {
      ext_0 <- rbind(
        ext_0,
        diagram[["extended"]]
      )
    } else {
      ext_1 <- rbind(
        ext_1,
        rev(diagram[["extended"]])
      )
    }
  }
  ex_diagram[["Ord0"]] <- ord_0
  ex_diagram[["Rel1"]] <- rel_1
  ex_diagram[["Ext0"]] <- ext_0
  ex_diagram[["Ext1"]] <- ext_1

  return(ex_diagram)
}

testMinimality <- function(vertices,
                           heights,
                           direction,
                           colinear_cond = 1e-8) {
  v_k <- vertices[2, ]
  v_prev <- vertices[1, ]
  v_next <- vertices[3, ]

  h_k <- heights[2]
  h_prev <- heights[1]
  h_next <- heights[3]

  if (h_k > h_prev || h_k > h_next) {
    return(FALSE)
  }

  if (abs(h_k - h_prev) <= colinear_cond &&
    abs(h_k - h_next) <= colinear_cond) {
    return(testNormalVector(v_k, v_next, direction))
  } else {
    P <- c(
      0.5 * (v_next[1] + v_prev[1]),
      0.5 * (v_next[2] + v_prev[2])
    )

    delta <- (v_next[1] - v_k[1]) * (P[2] - v_k[2]) -
      (v_next[2] - v_k[2]) * (P[1] - v_k[1])

    if (delta != 0) {
      return(delta > 0)
    } else {
      stop("Numerical error in minimality test.
                 Could not determine minimality.")
    }
  }
}

testNormalVector <- function(v_1,
                             v_2,
                             direction) {
  normal_vect <- c(v_2[2] - v_1[2], v_2[1] - v_1[1]) * c(-1, 1)

  return(dotProduct(direction, normal_vect) > 0)
}
