#' @title Get double bark thickness of tree at given height Hx
#' @description this function calculates double bark thickness in given height
#' for a given tree
#' @param tree either a data.frame or a list containing the variables needed to
#' decribe a tree, i.e. spp, D1, H, and optionally H1, D2, H2. Additionally,
#' parameter \code{Hx} might be directly given via \code{tree}. See
#' \code{\link{buildTree}} for more details and parameter \code{mapping} for
#'  mapping of variable names.
#' @param Hx height in tree for which double bark thickness is required
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names. See details.
#' @param ... passing arguments to methods.
#' @details
#' This function returns double bark thickness in given height \code{Hx} in stem
#' taper (hence, it depends on the diameter in given height). This can be added
#' to an diameter under bark to get the diameter over bark.
#' Parameter \code{tree} is able to take either a data.frame with correct
#' variables names or arbitrary names if \code{mapping} is provided to map the
#' data.frame names to the required names by \code{c("df-colname" = "var-name")}
#' or to take a named list.
#'
#' @return vector of double bark thickness given height \code{Hx} inside stem taper
#' @examples
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 25), H = c(25, 25), Hx = c(1.3, 22.248))
#' getBark(tree)
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 25), Hx = c(1.3, 22.248))
#' getBark(tree, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 25))
#' Hx <- list(Hx = c(1.3, 22.248))
#' getBark(tree = tree, Hx = Hx, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getBark <- function (tree, ...) {
  UseMethod("getBark", tree)
}

#' @describeIn getBark transforming \code{data.frame} before calling
#' \code{getBark} using \code{buildTree}
#' @export
getBark.data.frame <- function(tree, Hx = NULL, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "bark",
                    vars = Hx,
                    mapping = mapping)
  getBark(tree, Hx = Hx, mapping = mapping)
}

#' @describeIn getBark transforming \code{list} before calling
#' \code{getBark} using \code{buildTree}
#' @export
getBark.list <- function(tree, Hx = NULL, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "bark",
                    vars = Hx,
                    mapping = mapping)
  getBark(tree, Hx = Hx, mapping = mapping)
}

#' @describeIn getBark class method for class \code{datBDAT}
#' @export
getBark.datBDAT <- function(tree, Hx = NULL, mapping = NULL, ...) {
  if (!("datBDAT.bark" %in% class(tree)) | !is.null(Hx)) {
    tree <- buildTree(tree, check = "bark", vars = Hx, mapping = mapping)
  }

  nc <- ifelse(!is.null(Hx), length(Hx), 1) # number of columns for final matrix
  nr <- nrow(tree) / nc # number of rows for final matrix

  if ("datBDAT.bark" %in% class(tree)) {

    ## get double bark thickness
    n <- nrow(tree)
    res <- as.vector(
      .Fortran("vBDATRinde2Hx",
        n = as.integer(n),
        vBDATBArtNr = as.integer(tree$spp),
        vD1 = as.single(tree$D1),
        vH1 = as.single(tree$H1),
        vD2 = as.single(tree$D2),
        vH2 = as.single(tree$H2),
        vHges = as.single(tree$H),
        vHx = as.single(tree$Hx),
        vIErr = as.integer(rep(0, n)),
        vRINDE2Hx = as.single(rep(0, n)),
        PACKAGE = "rBDAT"
      )$vRINDE2Hx
    )
    res <- ifelse(tree$D1 <= 0, 0, res)
  }

  return(matrix(res, nrow = nr, ncol = nc)[, , drop = TRUE])
}
