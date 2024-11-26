#' @title Get diameter in given height inside tree taper
#' @description this function calculates the diameter inside or outside bark of
#' in given height for a given tree
#' @param tree either a data.frame, a list or an object of class \code{datBDAT}
#' containing the variables needed to decribe a tree, i.e. spp, D1, H, and
#' optionally H1, D2, H2. See \code{\link{buildTree}} for details and parameter
#' \code{mapping} for mapping of variable names
#' @param Hx height in tree for which diameter over or under bark is required;
#' defaults to NULL
#' @param bark logical, if TRUE returned diameter \code{Dx} is over bark,
#' if FALSE returned diameter is under bark. Coerced to logical by
#' \code{\link{as.logical}(bark[1])}.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names. See details.
#' @param ... passing arguments to methods.
#' @details if tree does not includes variable Hx, a full outer join is
#' generated between both
#'
#' @return a matrix with one row for each tree and one column for each \code{Hx}
#' given, holding the diameter over or under bark of provided height \code{Hx}
#' inside stem taper. The matrix is simplified by \code{[,,drop=TRUE]},
#' especially if \code{Hx=NULL}.
#' @examples
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 30), H = c(25, 25), Hx = c(1.3, 22.248))
#' getDiameter(tree, bark = TRUE)
#' getDiameter(tree, bark = FALSE)
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 30), h = c(25, 25), Hx = c(1.3, 22.248))
#' getDiameter(tree, bark = TRUE, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 30), h = c(25, 25))
#' Hx <- c(1.3, 22.248)
#' getDiameter(tree, Hx = Hx, bark = TRUE, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getDiameter <- function (tree, ...) {
  UseMethod("getDiameter", tree)
}

#' @describeIn getDiameter transforming \code{data.frame} before calling
#' \code{getDiameter} using \code{buildTree}
#' @export
getDiameter.data.frame <- function(tree, Hx = NULL, bark = TRUE,
                                   mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "diameter",
                    vars = Hx,
                    mapping = mapping)
  getDiameter(tree, Hx = Hx, bark = bark, mapping = mapping)
}

#' @describeIn getDiameter transforming \code{list} before calling
#' \code{getDiameter} using \code{buildTree}
#' @export
getDiameter.list <- function(tree, Hx = NULL, bark = TRUE, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "diameter",
                    vars = Hx,
                    mapping = mapping)
  getDiameter(tree, Hx = Hx, bark = bark, mapping = mapping)
}

#' @describeIn getDiameter class method for class \code{datBDAT}
#' @export
getDiameter.datBDAT <- function(tree, Hx = NULL, bark = TRUE,
                                mapping = NULL, ...) {
  if (!("datBDAT.diameter" %in% class(tree)) | !is.null(Hx)) {
    tree <- buildTree(tree, check = "diameter", vars = Hx, mapping = mapping)
  }

  bark <- as.logical(bark[1])

  nc <- ifelse(!is.null(Hx), length(Hx), 1) # number of columns for final matrix
  nr <- nrow(tree) / nc # number of rows for final matrix

  if ("datBDAT.diameter" %in% class(tree)) {

    ## calculate diameter
    n <- nrow(tree)
    if (isTRUE(bark)) {
      res <- as.vector(
        .Fortran("vbdatdmrhx",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vHges = as.single(tree$H),
          vHx = as.single(tree$Hx),
          vIErr = as.integer(rep(0, n)),
          vdmrhx = as.single(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vdmrhx
      )
      res <- ifelse(tree$D1 <= 0, 0, res)
    }
    else {
      res <- as.vector(
        .Fortran("vbdatdorhx",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vHges = as.single(tree$H),
          vHx = as.single(tree$Hx),
          vIErr = as.integer(rep(0, n)),
          vdorhx = as.single(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vdorhx
      )
      res <- ifelse(tree$D1 <= 0, 0, res)
    }
  }
  else {
    stop("data not of appropriate class!")
  }

  return(matrix(res, nrow = nr, ncol = nc)[, , drop = TRUE])
}
