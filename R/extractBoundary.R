#' Extract the boundary curves of a binary image
#'
#' The function [extractBoundary()] extracts and stores each boundary curve
#'  of a binary image.
#'
#' Given a binary image (pixel values 0 and 1) \code{extractBoundary}
#'  will trace all boundary curves along the midpoints of pixel edges between
#'  foreground and background regions. Curves which bound the outside of a
#'  foreground region are traced in the anticlockwise direction while curves
#'  which bound the inside of a foreground region are traced in the clockwise
#'  direction. The image is padded to ensure outermost pixels are part of the
#'  background.
#'
#' @param img A binary image of class `cimg' read in using
#'  [imager::load.image()].
#' @param background The value of the background, 0 for black and 1 for white.
#'  The default value is \code{1}.
#' @param saveOutput If TRUE, will save output as a \code{.RDS} file to
#'  directory specified by outputDir. The default values is `FALSE`.
#' @param outputDir The directory to save the output. If saveOutput is TRUE
#'  and no directory is specified, saves to working directory. The default
#'  is `NULL`.
#' @param fName The name of the output file saved. If saveOutput is TRUE and
#'  no filename specified, prompts user for filename. The default value is
#'  `NULL`.
#' @param verbose If TRUE, prints indicators of progress throughout. The default
#'  value is `TRUE`.
#' @return A list containing the points on the boundary curves structured as:
#'    * The first entry is a vector \code{(m,n)} containing the number of
#'      anticlockwise (m) and clockwise (n) curves.
#'    * Entries \code{2} to \code{m+1}  contain the anticlockwise curves.
#'    * Entries \code{m+2} to \code{n+1} contain the clockwise curves.
#' @seealso [multiExtractBoundary()], [imager::load.image()]
#' @export
extractBoundary <- function(img,
                            background = 1,
                            saveOutput = FALSE,
                            outputDir = NULL,
                            fName = NULL,
                            verbose = TRUE) {
  if (saveOutput) {
    if (is.null(outputDir)) {
      outputDir <- getwd()
    } else if (!dir.exists(outputDir)) {
      outputDir <- getwd()
    }
    if (is.null(fName)) {
      fName <- readline(prompt = "Please provide a filename for save: ")
    }
    out_file <- paste(outputDir, "/", fName, ".RDS", sep = "")
  }
  # Pad image
  img <- imager::pad(img, nPix = 1, axes = "x", pos = -1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "x", pos = 1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "y", pos = -1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "y", pos = 1, val = background)

  img_matrix <- componentLabelling(img, background, verbose)

  boundary <- boundaryTrace(img_matrix, verbose)

  if (saveOutput) {
    saveRDS(boundary, file = out_file)
    if (verbose) {
      cat("Successfully saved ", out_file, "\n", sep = "")
    }
  } else {
    if (verbose) {
      cat("Boundary extraction successful.\n")
    }
    return(boundary)
  }
}

#' Extract the boundary curves of binary images in a given directory
#'
#' Runs [extractBoundary()] on all images of a specified type in a
#' given directory. Images must be binary (pixel values 0 and 1).
#'
#' @param inputDir The directory containing the images.
#' @param imgType The image file type for input files, specified without a
#' ".". The default value is `pn g`.
#' @inheritParams extractBoundary
#' @return A list for each image containing the points of the boundary curves,
#' structured as in [extractBoundary()].
#' @seealso [extractBoundary()]
#' @export
multiExtractBoundary <- function(inputDir,
                                 imgType = "png",
                                 background = 1,
                                 saveOutput = FALSE,
                                 outputDir = NULL,
                                 verbose = TRUE) {
  f_type <- paste("*.", imgType, sep = "")
  files <- list.files(
    path = inputDir,
    pattern = f_type,
    full.names = TRUE,
    recursive = FALSE
  )

  if (saveOutput) {
    if (is.null(outputDir)) {
      outputDIr <- getwd()
    } else if (!dir.exists(outputDir)) {
      outputDir <- getwd()
    }
  } else {
    boundaries <- vector(mode = "list", length = length(files))
  }

  for (i in seq_along(files)) {
    if (verbose) {
      cat("### Commencing", files[[i]], "\n", sep = " ")
    }

    img <- imager::load.image(files[[i]])

    if (saveOutput) {
      f_string <- strsplit(files[[i]], "/", fixed = TRUE)
      fName <- strsplit(f_string[[1]][length(f_string[[1]])], ".",
        fixed = TRUE
      )[[1]][1]
      extractBoundary(
        img = img,
        background = background,
        saveOutput = TRUE,
        outputDir = outputDir,
        fName = fName,
        verbose = verbose
      )
    } else {
      boundaries[[i]] <- extractBoundary(
        img = img,
        background = background,
        verbose = verbose
      )
    }
  }

  if (saveOutput) {
    cat("All boundaries successfully saved in:\n",
        outputDir, "\n", sep = "")
  } else {
    return(boundaries)
  }
}

componentLabelling <- function(img,
                               background,
                               verbose) {
  w <- imager::width(img)
  h <- imager::height(img)

  img_matrix <- matrix(img, nrow = h, ncol = w, byrow = TRUE)
  img_matrix <- apply(img_matrix, 2, rev)

  if (any(img_matrix > 1) || any(img_matrix < 0)) {
    stop("Image not binary. Pixel values must be either 0 or 1.")
  }

  labels_fg <- vector(mode = "list")
  labels_bg <- vector(mode = "list")

  current_fg <- 1
  current_bg <- -1

  foreground <- 1 - background
  # image is scanned left-to-right, top-to-bottom
  for (i in 1:h) {
    for (j in 1:w) {
      if (img_matrix[i, j] == foreground) {
        # Hit a foreground pixel. Assess the mask (X's):
        # X X X
        # X * o
        # o o o
        mask <- c(
          img_matrix[i, j - 1],
          img_matrix[i - 1, j - 1],
          img_matrix[i - 1, j],
          img_matrix[i - 1, j + 1]
        )

        if (mask[3] > 0) {
          img_matrix[i, j] <- mask[3]
        } else if (mask[1] > 0) {
          img_matrix[i, j] <- mask[1]
          if (mask[4] > 0) {
            labels_fg <- resolveLabels(mask[4],
              mask[1],
              labels_fg,
              bg = FALSE
            )
          }
        } else if (mask[2] > 0) {
          img_matrix[i, j] <- mask[2]
          if (mask[4] > 0) {
            labels_fg <- resolveLabels(mask[2],
              mask[4],
              labels_fg,
              bg = FALSE
            )
          }
        } else if (mask[4] > 0) {
          img_matrix[i, j] <- mask[4]
        } else {
          img_matrix[i, j] <- current_fg
          labels_fg[[current_fg]] <- current_fg
          current_fg <- current_fg + 1
        }
      } else {
        # Found a background pixel
        mask <- c(0, 0)
        if (i > 1) {
          mask[2] <- img_matrix[i - 1, j]
        }
        if (j > 1) {
          mask[1] <- img_matrix[i, j - 1]
        }

        if (mask[2] < 0) {
          img_matrix[i, j] <- mask[2]
          if (mask[1] < 0) {
            labels_bg <- resolveLabels(-mask[1],
              -mask[2],
              labels_bg,
              bg = TRUE
            )
          }
        } else if (mask[1] < 0) {
          img_matrix[i, j] <- mask[1]
        } else {
          img_matrix[i, j] <- current_bg
          x <- -current_bg
          labels_bg[[x]] <- current_bg
          current_bg <- current_bg - 1
        }
      }
    }
  }
  if (verbose) {
    cat("Completed first raster for component labelling.\n",
        "Now constructing representative tables.\n", sep = "")
  }

  for (i in seq_len(length(labels_fg))) {
    for (lab in labels_fg[[i]]) {
      labels_fg[[lab]] <- unique(c(
        labels_fg[[lab]],
        labels_fg[[i]]
      ))
    }
  }

  for (i in seq_len(length(labels_bg))) {
    for (lab in labels_bg[[i]]) {
      x <- -lab
      labels_bg[[x]] <- unique(c(
        labels_bg[[x]],
        labels_bg[[i]]
      ))
    }
  }

  rep_table_fg <- vector()
  rep_table_bg <- vector()

  for (i in seq_len(length(labels_fg))) {
    rep_table_fg <- append(rep_table_fg, min(labels_fg[[i]]))
  }
  for (i in seq_len(length(labels_bg))) {
    rep_table_bg <- append(rep_table_bg, max(labels_bg[[i]]))
  }

  lab <- 1
  for (i in 1:(current_fg - 1)) {
    if (rep_table_fg[i] == i) {
      rep_table_fg[i] <- lab
      lab <- lab + 1
    } else {
      rep_table_fg[i] <- rep_table_fg[rep_table_fg[i]]
    }
  }

  lab <- -1
  for (i in 1:(abs(current_bg) - 1)) {
    if (rep_table_bg[i] == -i) {
      rep_table_bg[i] <- lab
      lab <- lab - 1
    } else {
      x <- -rep_table_bg[i]
      rep_table_bg[i] <- rep_table_bg[x]
    }
  }

  if (verbose) {
    cat("Representative tables constructed. Now relabelling components.\n")
  }

  for (i in 1:h) {
    for (j in 1:w) {
      if (img_matrix[i, j] > 0) {
        img_matrix[i, j] <- rep_table_fg[img_matrix[i, j]]
      } else {
        x <- -img_matrix[i, j]
        img_matrix[i, j] <- rep_table_bg[x]
      }
    }
  }

  return(img_matrix)
}

resolveLabels <- function(l_1,
                          l_2,
                          labs,
                          bg) {
  if (bg) {
    labs[[l_1]] <- unique(c(labs[[l_1]], labs[[l_2]], -l_2))
    labs[[l_2]] <- unique(c(labs[[l_1]], labs[[l_2]], -l_1))
  } else {
    labs[[l_1]] <- unique(c(
      labs[[l_1]],
      labs[[l_2]],
      l_1,
      l_2
    ))
    labs[[l_2]] <- unique(c(
      labs[[l_1]],
      labs[[l_2]],
      l_1,
      l_2
    ))
  }
  return(labs)
}

boundaryTrace <- function(img_matrix,
                          verbose) {
  w <- ncol(img_matrix) - 1
  h <- nrow(img_matrix) - 1

  regions <- vector()

  extreme_labels <- c(min(img_matrix), max(img_matrix))
  final <- FALSE

  c_curves <- vector(mode = "list") # clockwise curve, i.e. inner
  ac_curves <- vector(mode = "list") # anticlockwise curve, i.e. outer.

  ac_count <- 0
  c_count <- 0

  for (i in 2:h) {
    for (j in 2:w) {
      if (img_matrix[i, j] != -1) {
        if (!any(regions == img_matrix[i, j])) {
          lab <- img_matrix[i, j]
          if (lab > 0) {
            ac_count <- ac_count + 1
            ac_curves[[ac_count]] <- traceCurve(img_matrix, i, j, lab)
          } else {
            c_count <- c_count + 1
            c_curves[[c_count]] <- traceCurve(img_matrix, i, j, lab)
          }
          regions <- append(
            regions,
            img_matrix[i, j]
          )
        }
        # if (any(regions == extreme_labels[1]) &&
        #  any(regions == extreme_labels[2])) {
        if (ac_count == extreme_labels[2] &&
          -c_count == extreme_labels[1] + 1) {
          final <- TRUE
          break
        }
      }
    }
    if (final) {
      break
    }
  }

  total_curves <- length(c_curves) + length(ac_curves)

  if (verbose) {
    cat("Successfully extracted", total_curves, "curves:\n",
      length(ac_curves),
      "anticlockwise curves\n",
      length(c_curves),
      "clockwise curves\n",
      sep = " "
    )
  }
  bdry_curves <- vector(mode = "list")
  bdry_curves[[1]] <- c(length(ac_curves), length(c_curves))
  bdry_curves <- do.call(c, list(bdry_curves, ac_curves, c_curves))

  return(bdry_curves)
}

traceCurve <- function(img_matrix, i, j, lab) {
  # Start trace from north edge of pixel
  p0 <- c(i - 0.5, j)

  p1 <- c(i, j)
  pc <- c(0, 0)

  # Keep track of which edge we are on
  # * 1 *
  # 2 * 4
  # * 3 *
  px_loc <- 1

  curve_pts <- matrix(data = p0, nrow = 1, ncol = 2)

  # Define directions (mod 8)
  # 3 2 1
  # 4 * 0
  # 5 6 7
  dir <- 7

  while (any(pc != p1)) {
    nbhd <- matrix(c(
      i, i - 1, i - 1, i - 1, i, i + 1, i + 1, i + 1,
      j + 1, j + 1, j, j - 1, j - 1, j - 1, j, j + 1
    ),
    nrow = 8,
    ncol = 2
    )

    if (dir %% 2 == 1) {
      dir <- dir - 1
    }

    for (k in ((dir + 7):(dir + 14) %% 8)) {
      nbr <- FALSE
      # adjust for R indexing
      x <- nbhd[k + 1, ][1]
      y <- nbhd[k + 1, ][2]

      if (img_matrix[x, y] == lab && any(c(x, y) != pc)) {
        nbr <- TRUE
        prev <- switch(px_loc,
          switch(k + 1,
            rbind(
              c(i, j - 0.5),
              c(i + 0.5, j),
              c(x + 0.5, y)
            ),
            rbind(
              c(i, j - 0.5),
              c(i + 0.5, j),
              c(i, j + 0.5),
              c(x + 0.5, y)
            ),
            c(),
            c(x, y + 0.5),
            c(x - 0.5, y),
            rbind(
              c(i, j - 0.5),
              c(x - 0.5, y)
            ),
            rbind(
              c(i, j - 0.5),
              c(x, y - 0.5)
            ),
            rbind(
              c(i, j - 0.5),
              c(i + 0.5, j),
              c(x, y - 0.5)
            )
          ),
          switch(k + 1,
            rbind(
              c(i + 0.5, j),
              c(x + 0.5, y)
            ),
            rbind(
              c(i + 0.5, j),
              c(i, j + 0.5),
              c(x + 0.5, y)
            ),
            rbind(
              c(i + 0.5, j),
              c(i, j + 0.5),
              c(x, j + 0.5)
            ),
            rbind(
              c(i + 0.5, j),
              c(i, j + 0.5),
              c(i - 0.5, j),
              c(x, y + 0.5)
            ),
            c(),
            c(x - 0.5, y),
            c(x, y - 0.5),
            rbind(
              c(i + 0.5, j),
              c(x, y - 0.5)
            )
          ),
          switch(k + 1,
            c(x + 0.5, y),
            rbind(
              c(i, j + 0.5),
              c(x + 0.5, y)
            ),
            rbind(
              c(i, j + 0.5),
              c(x, y + 0.5)
            ),
            rbind(
              c(i, j + 0.5),
              c(i - 0.5, j),
              c(x, y + 0.5)
            ),
            rbind(
              c(i, j + 0.5),
              c(i - 0.5, j),
              c(x - 0.5, y)
            ),
            rbind(
              c(i, j + 0.5),
              c(i - 0.5, j),
              c(i, j - 0.5),
              c(x - 0.5, y)
            ),
            c(),
            c(x, y - 0.5)
          ),
          switch(k + 1,
            c(),
            c(x + 0.5, y),
            c(x, y + 0.5),
            rbind(
              c(i - 0.5, j),
              c(x, y + 0.5)
            ),
            rbind(
              c(i - 0.5, j),
              c(x - 0.5, y)
            ),
            rbind(
              c(i - 0.5, j),
              c(i, j - 0.5),
              c(x - 0.5, y)
            ),
            rbind(
              c(i - 0.5, j),
              c(i, j - 0.5),
              c(x, y - 0.5)
            ),
            rbind(
              c(i - 0.5, j),
              c(i, j - 0.5),
              c(i + 0.5, j),
              c(x, y - 0.5)
            )
          )
        )

        curve_pts <- unname(rbind(curve_pts, prev))
        px_loc <- switch(px_loc,
          switch(k + 1,
            3,
            3,
            NULL,
            4,
            1,
            1,
            2,
            2
          ),
          switch(k + 1,
            3,
            3,
            4,
            4,
            NULL,
            1,
            2,
            2
          ),
          switch(k + 1,
            3,
            3,
            4,
            4,
            1,
            1,
            NULL,
            2
          ),
          switch(k + 1,
            NULL,
            3,
            4,
            4,
            1,
            1,
            2,
            2
          )
        )
        dir <- k
        i <- x
        j <- y
        pc <- c(i, j)
        break
      }
    }
    if (!nbr) {
      # The region is a single pixel
      curve_pts <- rbind(
        curve_pts,
        c(i, j - 0.5),
        c(i + 0.5, j),
        c(i, j + 0.5)
      )
      pc <- p1
    }
  }

  curve_pts[, c(1, 2)] <- curve_pts[, c(2, 1)]
  if (lab > 0) {
    curve_pts[, 1] <- rev(curve_pts[, 1])
    curve_pts[, 2] <- rev(curve_pts[, 2])
  } 
  
  return(curve_pts)
}
