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
                            background = 0,
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

  img_matrix <- prepImage(img)

  h <- dim(img_matrix)[1]
  w <- dim(img_matrix)[2]

  indices <- expand.grid(1:(w - 1), 1:(h - 1))

  boundaries <- apply(indices, 1, getEdge, img_matrix = img_matrix)
  boundaries[sapply(boundaries, is.null)] <- NULL

  b <- matrix(unlist(boundaries), ncol = 4, byrow = TRUE)

  a <- vector(mode = "logical", length = nrow(b))
  remaining_ind <- seq_len(nrow(b))

  for (i in seq_len(nrow(b))) {
    n <- findNext(i, b, remaining_ind)
    a[i] <- remaining_ind[n]
    remaining_ind <- remaining_ind[-n]
  }

  c_curves <- vector(mode = "list") # clockwise curve, i.e. inner
  ac_curves <- vector(mode = "list") # anticlockwise curve, i.e. outer.

  ac_count <- 0
  c_count <- 0

  assigned_check <- vector(mode = "logical", length = length(a))

  for (i in seq_along(a)) {
    if (!assigned_check[i]) {
      ordering <- vector()
      ordering <- addNext(i, a)
      assigned_check <- replace(assigned_check, ordering, TRUE)
      bound <- b[ordering, 1:2]
      if (orientation(bound)) {
        c_count <- c_count + 1
        c_curves[[c_count]] <- bound
      } else {
        ac_count <- ac_count + 1
        ac_curves[[ac_count]] <- bound
      }
    }
  }
  total_curves <- length(c_curves) + length(ac_curves)

  if (TRUE) {
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

  if (saveOutput) {
    saveRDS(bdry_curves, file = out_file)
    if (verbose) {
      cat("Successfully saved ", out_file, "\n", sep = "")
    }
  } else {
    if (verbose) {
      cat("Boundary extraction successful.\n")
    }
    return(bdry_curves)
  }
}

prepImage <- function(img,
                      background = 1) {

  if (any(img > 1) || any(img < 0)) {
    stop("Image not binary. Pixel values must be either 0 or 1.")
  }


  # Pad image
  img <- imager::pad(img, nPix = 1, axes = "x", pos = -1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "x", pos = 1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "y", pos = -1, val = background)
  img <- imager::pad(img, nPix = 1, axes = "y", pos = 1, val = background)

  w <- imager::width(img)
  h <- imager::height(img)

  img_matrix <- 1 - matrix(img, nrow = h, ncol = w, byrow = TRUE)
  img_matrix <- apply(img_matrix, 2, rev)


  if (any(img_matrix > 1) || any(img_matrix < 0)) {
    stop("Image not binary. Pixel values must be either 0 or 1.")
  }

  return(img_matrix)
}



getEdge <- function(img_matrix, index, background = 0, verbose) {
  width_value <- index[[1]]
  height_value <- index[[2]]

  block <- c(
    img_matrix[height_value, width_value],
    img_matrix[height_value, width_value + 1],
    img_matrix[height_value + 1, width_value],
    img_matrix[height_value + 1, width_value + 1]
  )

  edge_mat <- addEdge(block, height_value, width_value, background, verbose)
  return(edge_mat)
}

addEdge <- function(block, i, j,
                    background = 0,
                    verbose) {
  table <- list(
    c(),                                                            # 0000
    c(0.5, 1, 1, 0.5) + c(j, i, j, i),                              # 0001
    c(0, 0.5, 0.5, 1) + c(j, i, j, i),                              # 0010
    c(0, 0.5, 1, 0.5) + c(j, i, j, i),                              # 0011
    c(1, 0.5, 0.5, 0) + c(j, i, j, i),                              # 0100
    c(0.5, 1, 0.5, 0) + c(j, i, j, i),                              # 0101
    c(0, 0.5, 0.5, 0, 1, 0.5, 0.5, 1) + c(j, i, j, i, j, i, j, i),  # 0110
    c(0, 0.5, 0.5, 0) + c(j, i, j, i),                              # 0111
    c(0.5, 0, 0, 0.5) + c(j, i, j, i),                              # 1000
    c(0.5, 1, 0, 0.5, 0.5, 0, 1, 0.5) + c(j, i, j, i, j, i, j, i),  # 1001
    c(0.5, 0, 0.5, 1) + c(j, i, j, i),                              # 1010
    c(0.5, 0, 1, 0.5) + c(j, i, j, i),                              # 1011
    c(1, 0.5, 0, 0.5) + c(j, i, j, i),                              # 1100
    c(0.5, 1, 0, 0.5) + c(j, i, j, i),                              # 1101
    c(1, 0.5, 0.5, 1) + c(j, i, j, i),                              # 1110
    c()                                                             # 1111
  )

  index <- strtoi(paste(block, collapse = ""), base = 2) + 1

  return(table[[index]])
}


findNext <- function(e, edges, index_table) {
  end <- edges[e, 3:4]
  for (k in seq_along(index_table)) {
    if (identical(end, edges[index_table[k], 1:2])) {
      return(k)
    }
  }
  stop("Unable to findNext")
}

addNext <- function(i, a) {
  bound_indices <- i
  new <- a[i]
  while (new != i) {
    bound_indices <- c(bound_indices, new)
    new <- a[new]
  }

  return(c(bound_indices, i))
}

orientation <- function(bound) {
  xa <- bound[nrow(bound) - 1, 1]
  ya <- bound[nrow(bound) - 1, 2]
  xb <- bound[1, 1]
  yb <- bound[1, 2]
  xc <- bound[2, 1]
  yc <- bound[2, 2]

  det <- (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)
  return(det < 0)
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
      outputDir <- getwd()
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
      fName <- strsplit(f_string[[1]][length(f_string[[1]])], ".", fixed = TRUE
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
