#' @title Plot taper curve of a tree
#' @description creating a plot of the taper curve of a tree, over or under bark
#' @param x an object of class 'datBDAT'
#' @param bark either NULL or logical; if TRUE taper curve over bark is plotted,
#' if FALSE taper curve under bark is plotted; if NULL, both are plotted
#' @param col.bark color to be used for plot of bark, if plot of taper curve
#' over and under bark is requested
#' @param legend logical, if legend should be added
#' @param assort assortments produced by \code{getAssortment(, value="merge")}
#' @param ... further arguments for \code{plot} and \code{points}
#' @details Creates graphics of the taper curve of trees. Either over bark or
#' under bark, or both. Elements design can partly be chosen. If assortments are
#' given, these are added to the plot. Doing that, the assortment bottom and top
#' position is indicated by a vertical line and mid-diameter is shown as a point
#' with vertical dashed line. N.B. the mid-diameter shown is under bark and
#' rounded downwards for 0.5 cm if mid-diameter < 20 and for 0.75 cm if bigger.
#' Assortment volume is calculated using this diameter according to the legal
#' rules for roundwood assortments (formerly german Forst-HKS and now RVR).
#' Additionally, assortment names are indicated.
#' One can provide assortment names in a column of \code{assort} named
#' 'assortname', which will be used if available, otherwise the 'Sort'-column
#'  will be used. See Examples.
#' @return No return value, called for side effects
#' @examples
#' ## plotting the taper curve of a tree
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 1))
#' t <- data.frame(spp = 1, D1 = 40, H = 35)
#' tree <- buildTree(tree = t)
#' plot(tree, type = "l", las = 1, legend = TRUE)
#' plot(tree, bark = TRUE, las = 1)
#' plot(tree, bark = FALSE, las = 1)
#' t <- data.frame(spp = c(1, 1), D1 = c(40, 35), H = c(35, 30))
#' tree <- buildTree(tree = t)
#' plot(tree, bark = FALSE, las = 1, legend = TRUE)
#' plot(tree, bark = TRUE, las = 1, legend = TRUE)
#'
#' t <- data.frame(spp = c(1, 8), D1 = 40, H = 35)
#' tree <- buildTree(tree = t)
#' plot(tree, bark = NULL, las = 1, col.bark = "blue", legend = TRUE)
#' plot(tree[1, ], main = getSpeciesCode(tree[1, ]$spp, out = "long"))
#' plot(tree[2, ], main = getSpeciesCode(tree[2, ]$spp, out = "scientific"))
#' par(mfrow = c(2, 1))
#' plot(tree, bark = TRUE, las = 1)
#'
#' ## now add assortments into taper curve
#' par(mfrow = c(1, 1))
#' ass <- getAssortment(tree, sort = list(lX = 1, fixN = 2, fixL = 4, fixA = 10))
#' plot(tree, assort = ass)
#' plot(tree, bark = FALSE, assort = ass)
#' plot(tree, bark = FALSE, assort = ass, legend = TRUE)
#' plot(tree[1, ], assort = ass[ass$tree == 1, ], main = "first tree in subset")
#' plot(tree[2, ], assort = ass[ass$tree == 2, ], main = "second tree in subset")
#'
#' ## adding own assortment labels using column 'assortname'
#' ass$assortname <- ifelse(grepl("Fix", ass$Sort), paste0("Fix:", ass$length), ass$Sort)
#' plot(tree, assort = ass)
#' par(oldpar)
#' @importFrom graphics abline lines plot points text
#' @export

plot.datBDAT <- function(x, bark = NULL, col.bark = NULL, legend = FALSE,
                         assort = NULL, ...) {
  ## catch call
  mc <- match.call()
  mc[[1]] <- as.name("plot")
  mc$bark <- NULL # remove from call to plot and points
  mc$col.bark <- NULL
  mc$legend <- NULL
  mc$assort <- NULL

  ## for each tree do...
  for (i in seq(along = x$D1)) {
    mci <- mc # for each tree a fresh call object to be manipulated
    if (is.null(mci$main)) {
      mci$main <- paste("spp=", x[i, ]$spp,
        ", D1=", round(x[i, ]$D1, 1),
        ", H=", round(x[i, ]$H, 1),
        sep = ""
      )
    }
    col <- ifelse(is.null(mci$col), "black", mci$col)
    if (is.null(mci$col)) mci$col <- col
    if (is.null(mci$xlab)) mci$xlab <- "H [m]"
    if (is.null(mci$ylab)) mci$ylab <- "D [cm]"
    if (is.null(mci$type)) mci$type <- "l"
    ## plotting
    plotx <- seq(0, x[i, ]$H, 0.1)
    mci$x <- plotx
    ## if bark is NULL, use bark==TRUE
    ploty <- getDiameter(x[i, ], Hx = plotx, bark = ifelse(is.null(bark), TRUE, bark))
    mci$y <- ploty
    eval(mci, parent.frame()) # plot the taper curve
    abline(v = 0, h = 0)

    if (is.null(bark)) {
      col.bark <- ifelse(is.null(col.bark), "grey", col.bark)
      mci[[1]] <- as.name("points")
      mci$col <- col.bark
      mci$y <- getDiameter(x[i, ], Hx = plotx, bark = FALSE)
      eval(mci, parent.frame()) # add lines for bark taper curve via 'points()'
      if (legend == TRUE) {
        legend("topright",
          legend = c("taper curve over bark", "taper curve under bark"),
          col = c(col, col.bark), lty = 1
        )
      }
    } else {
      if (legend == TRUE) {
        legend("topright",
          legend = ifelse(bark == TRUE, "taper curve over bark",
            "taper curve under bark"
          ),
          lty = 1, col = col
        )
      }
    }
    if (!is.null(assort)) {
      j <- unique(assort$tree)[i] # select correct assortments for actual tree
      asstmp <- assort[assort$tree == j & assort$length != 0, ]
      if (nrow(asstmp) > 0) {
        htmp <- unique(round(c(asstmp$height, asstmp$height + asstmp$length), 4))
        dtmp <- getDiameter(x[i, ],
          Hx = htmp,
          bark = ifelse(is.null(bark), TRUE, bark)
        )
        invisible(sapply(seq(length(htmp)), function(a) {
          lines(x = rep(htmp[a], 2), y = c(0, dtmp[a]))
        }))
        invisible(sapply(seq(nrow(asstmp)), function(a) {
          lines(
            x = rep(asstmp$height[a] + asstmp$length[a] / 2, 2),
            y = c(0, asstmp$midD[a]),
            lty = 2, col = "light grey"
          )
        }))
        points(x = asstmp$height + asstmp$length / 2, y = asstmp$midD, pch = 16)
        if (!is.null(asstmp$assortname)) {
          asslabel <- asstmp$assortname
        } else {
          asslabel <- asstmp$Sort
        }
        text(x = asstmp$height + asstmp$length / 2, y = 0, pos = 3, labels = asslabel)
      }
    }
  }
}
