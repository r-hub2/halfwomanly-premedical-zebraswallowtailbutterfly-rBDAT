#' @title Get height of given diameter inside tree taper
#' @description this function calculates the height of a given diameter inside
#' or outside bark for a given tree
#' @param tree either a data.frame or a list containing the variables needed to
#' decribe a tree, i.e. spp, D1, H, and optionally H1, D2, H2. See
#' \code{\link{buildTree}} for details and parameter \code{mapping} for
#'  mapping of variable names
#' @param Dx diameter of tree for which height is required; defaults to NULL
#' @param bark logical, if TRUE given diameter \code{Dx} is considered over
#' bark, if FALSE diameter is considered to be under bark. Coerced to logical by
#' \code{\link{as.logical}(bark[1])}.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names.
#' @param ... passing arguments to methods.
#' @details see \code{\link{buildTree}} for how to specify a tree object
#' @return a matrix with one row for each tree and one column for each \code{Dx}
#' given, holding the height of provided diameter \code{Dx} inside stem taper.
#' The matrix is simplified by \code{[,,drop=TRUE]}, especially if
#' \code{Dx=NULL}.
#' @examples
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 25), H = c(25, 20), Dx = c(7, 7))
#' getHeight(tree, bark = TRUE)
#' getHeight(tree, bark = FALSE)
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 20), Dx = c(7, 7))
#' getHeight(tree, bark = TRUE, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 20))
#' Dx <- c(7, 5)
#' getHeight(tree, Dx = Dx, bark = TRUE, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#'
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 25), H = c(25, 20), Dx = c(7, 7))
#' getHeight(tree, Dx = c(1:5), bark = TRUE)
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getHeight <- function (tree, ...) {
  UseMethod("getHeight", tree)
}

#' @describeIn getHeight transforming \code{data.frame} before calling
#' \code{getHeight} using \code{buildTree}
#' @export
getHeight.data.frame <- function(tree, Dx = NULL, bark = TRUE,
                                   mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "height",
                    vars = Dx,
                    mapping = mapping)
  getHeight(tree, Dx = Dx, bark = bark, mapping = mapping)
}

#' @describeIn getHeight transforming \code{list} before calling
#' \code{getHeight} using \code{buildTree}
#' @export
getHeight.list <- function(tree, Dx = NULL, bark = TRUE, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "height",
                    vars = Dx,
                    mapping = mapping)
  getHeight(tree, Dx = Dx, bark = bark, mapping = mapping)
}

#' @describeIn getHeight class method for class \code{datBDAT}
#' @export
getHeight.datBDAT <- function(tree, Dx = NULL, bark = TRUE,
                              mapping = NULL, ...) {
  if (!("datBDAT.height" %in% class(tree)) | !is.null(Dx)) {
    tree <- buildTree(tree, check = "height", vars = Dx, mapping = mapping)
  }

  bark <- as.logical(bark[1])

  nc <- ifelse(!is.null(Dx), length(Dx), 1) # number of columns for final matrix
  nr <- nrow(tree) / nc # number of rows for final matrix

  if ("datBDAT.height" %in% class(tree)) {

    ## get height of diameter inside stem
    n <- nrow(tree)
    if (isTRUE(bark)) {
      res <- as.vector(
        .Fortran("vbdathxdx",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vH = as.single(tree$H),
          vDx = as.single(tree$Dx),
          vHx = as.single(rep(0, n)),
          vIErr = as.integer(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vHx
      )
      res <- ifelse(tree$D1 <= 0, 0, res)
    }
    else {
      res <- as.vector(
        .Fortran("vbdathxdxor",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vH = as.single(tree$H),
          vDx = as.single(tree$Dx),
          vHxor = as.single(rep(0, n)),
          vIErr = as.integer(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vHxor
      )
      res <- ifelse(tree$D1 <= 0, 0, res)
    }
  }
  return(matrix(res, nrow = nr, ncol = nc)[, , drop = TRUE])
}
