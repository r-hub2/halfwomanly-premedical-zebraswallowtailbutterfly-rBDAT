#' @title Get estimated mean of form factor q03
#' @description this function returns for a given dbh-height-class and species
#' the estimated mean of the q03-distribution, which can be considered as a
#' form factor for different inventories (i.e. NFI 1, NSoG, IS08, NFI3).
#' @param tree either a data.frame or a list containing the variables needed to
#' decribe a tree, i.e. spp, D1, H, and optionally H1, D2, H2. See
#' \code{\link{buildTree}} for details and parameter \code{mapping} for
#'  mapping of variable names
#' @param inv integer, indicator for which inventory the form factor should be
#' returned; defaults to 1 (0=volume tables of Grundner-Schwappach 1921,
#' 1=NFI1, 2=NSoG, 3=IS08, 4=NFI3).
#' Any other value is silently set back to \code{inv=1}.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names.
#' @param ... passing arguments to methods.
#' @details
#' This function returns the estimated mean of the q03-distribution for a given
#' dbh-height-class and species. The q03-distribution is the ratio between the
#' diameter in 30\% of tree height and 5\% of tree height and describes the form
#' of the taper curve. This value can be returned for different inventories,
#' during which upper stem diameters have been sampled.
#'
#' Only spp, D1 and H are used to estimate mean q03, if provided and H1 equals
#' 1.3m, otherwise dbh is estimated using BDAT functions (i.e.
#' getDiameter(tree, Hx=1.3)).
#'
#'
#'
#' @return vector of form factor q03.
#' @examples
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 25), H = c(25, 25))
#' tree <- buildTree(tree, check = "form", vars = list(inv = 1))
#' str(tree)
#' getForm(tree) # default is 'inv=1'
#' getForm(tree, inv = 0) # taper form of BDAT from model fit
#' getForm(tree, inv = 1) # NFI 1 (without NSoG)
#' getForm(tree, inv = 2) # NSoG = New States of Germany
#' getForm(tree, inv = 3) # carbon inventory 2008
#' getForm(tree, inv = 4) # NFI 3 (reunificated Germany)

#' @useDynLib rBDAT, .registration=TRUE
#' @export
getForm <- function (tree, ...) {
  UseMethod("getForm", tree)
}

#' @describeIn getForm transforming \code{data.frame} before calling
#' \code{getForm} using \code{buildTree}
#' @export
getForm.data.frame <- function(tree, inv = NULL, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "form",
                    vars = inv,
                    mapping = mapping)
  getForm(tree, inv = inv, mapping = mapping)
}

#' @describeIn getForm transforming \code{list} before calling
#' \code{getForm} using \code{buildTree}
#' @export
getForm.list <- function(tree, inv = NULL, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "form",
                    vars = inv,
                    mapping = mapping)
  getForm(tree, inv = inv, mapping = mapping)
}

#' @describeIn getForm class method for class \code{datBDAT}
#' @export
getForm.datBDAT <- function(tree, inv = NULL, mapping = NULL, ...) {
  if (!("datBDAT.form" %in% class(tree)) | !is.null(inv)) {
    tree <- buildTree(tree, check = "form", vars = inv, mapping = mapping)
  }

  if ("datBDAT.form" %in% class(tree)) {

    ## get q03
    n <- nrow(tree)
    if (n >= 1) {
      res <- as.vector(
        .Fortran("vbdatformtarif",
          n = as.integer(n),
          vTarif = as.integer(tree$inv),
          vBDATBArtNr = as.integer(tree$spp),
          vD = as.single(tree$D13),
          vH = as.single(tree$H),
          vMwQ03BWI = as.single(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vMwQ03BWI
      )
      res <- ifelse(tree$D13 <= 0, 0, res)
    } else {
      res <- 0
    }
  }

  return(res)
}
