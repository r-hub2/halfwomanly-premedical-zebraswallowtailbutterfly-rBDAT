#' @title Calculate diameter in height Hx under bark.
#'
#' @description Function to call BDAT Fortran subroutine to calculate diameter
#' under bark in height Hx for specified tree/s.
#' @param BDATArtNr numeric vector of species code; see
#' \code{\link{getSpeciesCode}}.
#' @param D1 first measured diameter of tree [cm], e.g. diameter in breast
#' height
#' @param H1 measurement height of \code{D1} [m]
#' @param D2 second measured diameter of tree, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H2 measurement height of D2, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H total tree height [m]
#' @param Hx height in tree for which diameter under bark is required
#' @details conventional function interface for BDATDORHX. See
#' \code{\link{getDiameter}} for more details.
#' @return vector of diameters under bark
#' @examples
#' # simple call of function, with all parameters
#' BDATDORHX(1, 30, 1.3, 0, 0, 25, Hx = 1.3)
#' # same with variables
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' Hx <- 1.3
#' BDATDORHX(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H, Hx = Hx)
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATDORHX(getSpeciesCode(c("Fi", "Bu")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' ## english species names abbreviated
#' BDATDORHX(getSpeciesCode(c("NS", "BE")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' @seealso \code{\link{BDATDMRHX}} for BDAT routine calculating diameter over
#' bark, \code{\link{getDiameter}} for a function with a more convenient english
#' name, more options and a bark switch.
#' @export
BDATDORHX <- function(BDATArtNr, D1, H1, D2, H2, H, Hx) {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H, Hx = Hx
  )
  res <- getDiameter(tree = df, Hx = NULL, bark = F, mapping = NULL)
  return(res)
}
