#' @title Calculate double bark thickness
#'
#' @description BDAT-Function to get double bark thickness at given height Hx
#' @param BDATArtNr numeric vector of species code; see
#' \code{\link{getSpeciesCode}}.
#' @param D1 first measured diameter of tree [cm] e.g. diameter in breast height
#' @param H1 measurement height of \code{D1} [m]
#' @param D2 second measured diameter of tree, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H2 measurement height of D2, see \code{\link{buildTree}} for details
#' on how to specify different taper forms
#' @param H total tree height [m]
#' @param Hx height for which double bark thickness is required [m]
#' @details This function returns double bark thickness in given height
#' \code{Hx} in stem
#' taper (hence, it depends on the diameter in given height). This can be added
#' onto a diameter under bark to receive diameter over bark.
#'
#' @return vector of double bark thickness given height \code{Hx} inside stem
#' taper.
#' @examples
#' ## simple call of function, with all parameters given
#' BDATRINDE2HX(1, 30, 1.3, 0, 0, 25, Hx = 1.3)
#'
#' ## same with variables
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' Hx <- 1.3
#' BDATRINDE2HX(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H, Hx = Hx)
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATRINDE2HX(getSpeciesCode(c("Fi", "Bu")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' ## english abbreviated
#' BDATRINDE2HX(getSpeciesCode(c("NS", "BE")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' @seealso \code{\link{getBark}} for a function with a convenient english name
#' @export
BDATRINDE2HX <- function(BDATArtNr, D1, H1 = 1.3, D2 = 0, H2 = 0, H, Hx) {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H, Hx = Hx
  )
  res <- getBark(tree = df, Hx = NULL, mapping = NULL)
  return(res)
}
