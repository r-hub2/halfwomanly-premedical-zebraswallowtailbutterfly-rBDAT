#' @title Get total aboveground biomass of tree
#' @description this function calculates total aboveground biomass
#' for a given tree
#' @param tree either a data.frame or a list containing the variables needed to
#' decribe a tree, i.e. spp, D1, H, and optionally H1, D2, H2. See
#' \code{\link{buildTree}} for details and parameter \code{mapping} for
#'  mapping of variable names
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names. See details.
#' @param ... passing arguments to methods.
#' @details
#' This function returns total aboveground biomass according to the biomass
#' functions developed for the german NFI3 (BWI3) on the basis of a field survey
#' covering whole Germany from 2007 to 2010. Hence, the base data differs from
#' those used to fit the taper functions. Neverteless and although the original
#' fortran code does not provide an interface to pass H2 or negative D2 values
#' to the functions (they were fitted with abolute D03-values, i.e. diameter
#' in 30% of tree height) this R function does allow for such use of D2 and H2
#' for reasons of consistency.
#'
#' The functions themselves are fitted for four species (Norway spruce, Scots
#' pine, European beech and Oak spp.) directly, and for another fourteen species
#' by generating pseudo-observations and fitting the functions to these. See
#' the enclosed report for details.
#'
#' The biomass functions integrate four different ranges of validity: (i) tree
#' height below 1.3m, (ii) dbh below 10cm, (iii) dbh between 10cm and the
#' 99th-quantile of the data and (iv) dbh above.
#'
#' @return vector of total aboveground biomass
#' @references Riedel, T. and G. Kaendler (2017). "Nationale
#' Treibhausgasberichterstattung: Neue Funktionen zur Schätzung der
#' oberirdischen Biomasse am Einzelbaum." Forstarchiv 88(2): 31-38.
#' @examples
#' tree <- data.frame(spp = c(1, 1), D1 = c(30, 25), H = c(25, 25))
#' getBiomass(tree)
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 25))
#' getBiomass(tree, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#'
#' tree <- data.frame(BDATCode = c(1, 1), dbh = c(30, 25), h = c(25, 25))
#' getBiomass(tree = tree, mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H"))
#'
#' ## this is standard usage of the fortran code
#' ## now changing the taper form via second diameter D2
#' tree <- data.frame(spp = 1, D1 = 30, D2 = c(26, 24), H = 25)
#' getBiomass(tree)
#'
#' ## the following usage of D2 and H2 w.r.t. the biomass functions are new
#' ## and unique to R
#' ## now changing the taper form via quantile of second diameter in H2
#' tree <- data.frame(spp = 1, D1 = 30, H2 = c(25, 50, 75), H = 25)
#' getBiomass(tree)
#'
#' ## now changing the taper form via form quotient (i.e. negative D2)
#' fq <- getForm(list(spp = 1, D1 = 30, H = 25), inv = 1:4)
#' tree <- data.frame(spp = 1, D1 = 30, D2 = -fq, H = 25)
#' getBiomass(tree) # biomass of an average tree according to different inventories
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getBiomass <- function (tree, ...) {
  UseMethod("getBiomass", tree)
}

#' @describeIn getBiomass transforming \code{data.frame} before calling
#' \code{getBiomass} using \code{buildTree}
#' @export
getBiomass.data.frame <- function(tree, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "biomass",
                    vars = NULL,
                    mapping = mapping)
  getBiomass(tree, mapping = mapping)
}

#' @describeIn getBiomass transforming \code{list} before calling
#' \code{getBiomass} using \code{buildTree}
#' @export
getBiomass.list <- function(tree, mapping = NULL, ...){
  tree <- buildTree(tree,
                    check = "biomass",
                    vars = NULL,
                    mapping = mapping)
  getBiomass(tree, mapping = mapping)
}

#' @describeIn getBiomass class method for class \code{datBDAT}
#' @export

getBiomass.datBDAT <- function(tree, mapping = NULL, ...) {
  if (!("datBDAT.biomass" %in% class(tree)) | is.null(tree$D13) | is.null(tree$D03)) {
    tree <- buildTree(tree,
      check = "biomass", vars = NULL,
      mapping = mapping
    )
  }

  ## since the fortran function does not evaluate negative D2 or any H2 values
  ## here the real D03 value in the correct height is calculated
  # tree$Hx <- 0.3 * tree$H
  # tree$D03 <- getDiameter(tree, bark = TRUE)

  if ("datBDAT.biomass" %in% class(tree)) {

    ## get total above-ground biomass
    res <- as.vector(
      .Fortran("vbiomasse",
        n = as.integer(nrow(tree)), # length of data
        as.integer(tree$spp), # BDAT-species-code
        as.single(tree$D13), # real D13
        as.single(tree$D03), # real D03
        as.single(0.3 * tree$H), # real Height of D03
        as.single(tree$H), # H
        vBiom = as.single(rep(0, nrow(tree))), # Biomasse
        PACKAGE = "rBDAT"
      )$vBiom
    )
  }
  return(res)
}
