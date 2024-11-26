#' @title Calculate diameter in height Hx over bark.
#'
#' @description Function to call BDAT Fortran subroutine to calculate diameter
#' over bark in height Hx for specified tree/s.
#' @param BDATArtNr numeric vector of species code; see
#' \code{\link{getSpeciesCode}}.
#' @param D1 first measured diameter of tree [cm], e.g. diameter in breast
#' height.
#' @param H1 measurement height of \code{D1} [m]
#' @param D2 second measured diameter of tree, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H2 measurement height of D2, see \code{\link{buildTree}} for details
#' on how to specify different taper forms
#' @param H total tree height [m]
#' @param Hx height in tree for which diameter over bark is required
#' @details conventional function interface for Fortran function BDATDMRHX. See
#' \code{\link{getDiameter}} for more details.
#' @return vector of diameters over bark
#' @examples
#' # simple call of function, with all parameters
#' BDATDMRHX(1, 30, 1.3, 0, 0, 25, Hx = 1.3)
#' # same with variables
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' Hx <- 1.3
#' BDATDMRHX(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H, Hx = Hx)
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATDMRHX(getSpeciesCode(c("Fi", "Bu")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' ## english species names abbreviated
#' BDATDMRHX(getSpeciesCode(c("NS", "BE")), 30, 0, 0, 0, H = 25, Hx = 1.3)
#' @seealso \code{\link{BDATDORHX}} for BDAT routine calculating diameter under
#' bark and \code{\link{getDiameter}} for a function with a more convenient
#' english name, more options and including a bark switch.
#' @export
BDATDMRHX <- function(BDATArtNr, D1, H1, D2, H2, H, Hx) {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H, Hx = Hx
  )
  res <- getDiameter(tree = df, Hx = NULL, bark = T, mapping = NULL)
  return(res)
}
