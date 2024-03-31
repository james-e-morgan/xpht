## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  cache = FALSE,
  verbose = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----dev='jpeg'---------------------------------------------------------------
library(imager)
gridA <- system.file("extdata/gridA.png", package = "XPHT")
gridg <- system.file("extdata/gridg.png", package = "XPHT")
imgA <- load.image(gridA)
imgg <- load.image(gridg)
par(mfrow = c(1, 2))
plot(imgA, axes = FALSE)
plot(imgg, axes = FALSE)

## -----------------------------------------------------------------------------
fpath <- system.file("extdata/letterA/courier-new.png", package = "XPHT")
img <- load.image(fpath)
plot(img)

## -----------------------------------------------------------------------------
library(XPHT)
boundaryA <- extractBoundary(img)

## -----------------------------------------------------------------------------
class(boundaryA)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(boundaryA[[2]], axes = FALSE, type = "l")
arrows(x0 = boundaryA[[2]][1, 1], y0 = boundaryA[[2]][1, 2],
       x1 = boundaryA[[2]][50, 1], y1 = boundaryA[[2]][50, 2],
       col = "blue")
plot(boundaryA[[3]], axes = FALSE, type = "l")
arrows(x0 = boundaryA[[3]][1, 1], y0 = boundaryA[[3]][1, 2],
       x1 = boundaryA[[3]][20, 1], y1 = boundaryA[[3]][20, 2],
       col = "blue")

## -----------------------------------------------------------------------------
fpath <- system.file("extdata/letterg/courier-new.png", package = "XPHT")
img <- load.image(fpath)
plot(img)

## -----------------------------------------------------------------------------
boundaryg <- extractBoundary(img, verbose = FALSE)
par(mfrow = c(1, 2))
plot(boundaryg[[2]], axes = FALSE, type = "l")
arrows(x0 = boundaryg[[2]][1, 1], y0 = boundaryg[[2]][1, 2],
       x1 = boundaryg[[2]][50, 1], y1 = boundaryg[[2]][50, 2],
       col = "blue")
plot(boundaryg[[3]], axes = FALSE, type = "l")
arrows(x0 = boundaryg[[3]][1, 1], y0 = boundaryg[[3]][1, 2],
       x1 = boundaryg[[3]][20, 1], y1 = boundaryg[[3]][20, 2],
       col = "blue")

## -----------------------------------------------------------------------------
fpath <- system.file("extdata/letterA", package = "XPHT")
ABoundaries <- multiExtractBoundary(fpath, imgType = "png", verbose = FALSE)
length(ABoundaries[1:10])
class(ABoundaries)

## -----------------------------------------------------------------------------
testBoundary <- ABoundaries[[10]]
par(mfrow = c(1, 2))
plot(testBoundary[[2]], axes = FALSE, type = "l")
arrows(x0 = testBoundary[[2]][1, 1], y0 = testBoundary[[2]][1, 2],
       x1 = testBoundary[[2]][50, 1], y1 = testBoundary[[2]][50, 2],
       col = "blue")
plot(testBoundary[[3]], axes = FALSE, type = "l")
arrows(x0 = testBoundary[[3]][1, 1], y0 = testBoundary[[3]][1, 2],
       x1 = testBoundary[[3]][20, 1], y1 = testBoundary[[3]][20, 2],
       col = "blue")

## -----------------------------------------------------------------------------
xphtA <- extendedPersistence(boundaryA, "A-courier-new", 8)

## -----------------------------------------------------------------------------
length(xphtA)

## -----------------------------------------------------------------------------
class(xphtA[[1]])
names(xphtA[[1]])
xphtA[[1]]

## -----------------------------------------------------------------------------
xphtg <- extendedPersistence(boundaryg, "g-courier-new", 8)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(xphtA[[1]], main = "Extended Persistence of A\nin Direction pi/4")
legend("bottomright", legend = c("Ord0", "Rel1", "Ess0", "Ess1"),
       col = getDefaultColours(), pch = 15:18)
plot(xphtg[[1]], main = "Extended Persistence of g\nin Direction pi/4")
legend("topright", legend = c("Ord0", "Rel1", "Ess0", "Ess1"),
       col = getDefaultColours(), pch = 15:18)

## -----------------------------------------------------------------------------
par(mfrow = c(2, 1), mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
plot(xphtA[[1]], barcode = TRUE,
     main = "Extended Persistence Barcode of A\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)
plot(xphtg[[1]], barcode = TRUE,
     main = "Extended Persistence Barcode of g\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)

## -----------------------------------------------------------------------------
centredA <- centreScaleDiagrams(xphtA, scale = FALSE)
centredg <- centreScaleDiagrams(xphtg, scale = FALSE)
par(mfrow = c(2, 1), mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
plot(centredA[[1]], barcode = TRUE,
     main = "Centred Extended Persistence Barcode of A\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)
plot(centredg[[1]], barcode = TRUE,
     main = "Centred Extended Persistence Barcode of g\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)

## -----------------------------------------------------------------------------
scaledA <- centreScaleDiagrams(xphtA, scale = TRUE, scaleConstant = 5)
scaledg <- centreScaleDiagrams(xphtg, scale = TRUE, scaleConstant = 5)
par(mfrow = c(2, 1), mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
plot(scaledA[[1]], barcode = TRUE, main = "Centred and Scaled Extended
     Persistence Barcode of A\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)
plot(scaledg[[1]], barcode = TRUE, main = "Centred and Scaled Extended
     Persistence Barcode of g\nin Direction pi/4")
legend("topright", inset = c(-0.3, 0),
       legend = c("Ord0", "Rel1", "Ess0", "Ess1"), lty = c(1, 1, 1, 1),
       col = getDefaultColours(), lwd = 2)

## -----------------------------------------------------------------------------
fpathA <- system.file("extdata/xpht32A/", package = "XPHT")
fpathg <- system.file("extdata/xpht32g/", package = "XPHT")
diagramsA <- stackDiagrams(fpathA)
diagramsg <- stackDiagrams(fpathg)

## -----------------------------------------------------------------------------
distMatrixA.aligned <- computeDistanceMatrix(diagramsA, 13, q = 2,
                                             aligned = TRUE, verbose = FALSE)
distMatrixg.unaligned <- computeDistanceMatrix(diagramsg, 13, q = 2,
                                               verbose = FALSE)

## -----------------------------------------------------------------------------
distMatrixA.aligned[1:6, 1:6]

## -----------------------------------------------------------------------------
plot.mds <- function(x, labels, main = NULL,
                     colour = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     legend.pos = "topleft", cex = 1.3, pch = 19,
                     cex.main = 1.5, cex.legend = 0.9, ...) {
  points.mds <- x$points

  limxy <- range(points.mds)
  limxy <- limxy + ((limxy[2] - limxy[1]) * 0.15) * c(-0.5, 0.5)

  par(mar = c(0.2, 0.8, 1.1, 0.8), ps = 10)
  plot(limxy, limxy, type = "n", axes = FALSE, frame = FALSE)
  rect(limxy[1], limxy[1], limxy[2], limxy[2], border = "#999999", lwd = 0.3)

  points(points.mds[, 1], points.mds[, 2], col = colour[as.integer(labels)],
         cex = cex, pch = pch)
  mtext(side = 3, main, cex = cex.main)

  labels.u <- unique(labels)
  legend.text <- as.character(labels.u)

  legend(legend.pos, legend = legend.text, inset = 0.03,
         col = colour[as.integer(labels.u)], bty = "n", pch = pch,
         cex = cex.legend)
}

## -----------------------------------------------------------------------------
fpath.matrix <- system.file("extdata/distMatg.RDS", package = "XPHT")
fpath.labels <- system.file("extdata/storey.txt", package = "XPHT")
distMat <- readRDS(fpath.matrix)
labels <- factor(scan(fpath.labels, character(), quote = ""))

fonts.mds <- cmdscale(distMat, eig = TRUE, k = 2)
plot.mds(fonts.mds, labels = labels,
         main = "MDS of Letter g Images using 2-Wasserstein Distance")

