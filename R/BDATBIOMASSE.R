#' @title Get total aboveground biomasse
#'
#' @description BDAT-Function to get total aboveground biomass.
#' @param BDATArtNr numeric vector of species code, see
#' \code{\link{getSpeciesCode}}
#' @param D1 first measured diameter of tree [cm], e.g. diameter in breast
#' height
#' @param H1 measurement height of \code{D1} [m]
#' @param D2 second measured diameter of tree, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H2 measurement height of D2, see \code{\link{buildTree}} for details
#' on how to specify different taper forms
#' @param H total tree height [m]
#' @details This function returns total aboveground biomass for given tree/s,
#' based on the biomass functions developed for the german NFI 3. See the
#' additional material for some german reference.
#'
#' @return vector of biomass for given trees
#' @examples
#' ## simple call of function, with all parameters given
#' BDATBIOMASSE(1, 30, 1.3, 0, 0, 25)
#'
#' ## same with variables
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' BDATBIOMASSE(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H)
#'
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATBIOMASSE(getSpeciesCode(c("Fi", "Bu")), 30, H = 25)
#' ## english abbreviated
#' BDATBIOMASSE(getSpeciesCode(c("NS", "BE")), 30, H = 25)
#' @seealso \code{\link{getBiomass}} for a function with more convenient english
#' name
#' @export
BDATBIOMASSE <- function(BDATArtNr, D1, H1 = 0, D2 = 0, H2 = 0, H) {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H
  )
  res <- getBiomass(tree = df, mapping = NULL)
  return(res)
}
