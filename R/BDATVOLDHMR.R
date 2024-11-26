#' @title Calculate (coarse) wood volume over bark of a tree up to given diameter
#'
#' @description BDAT-Function to get (coarse) wood volume over bark up to a
#' given diameter of one or many trees.
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
#' @param DHGrz diameter over bark up to which volume should be calculated
#' @param SekLng length of section over which the integral of taper form should
#' be applied, defaults to 2.0m.
#' @details Volume is calculated using BDAT Fortran routines. In particular,
#' \code{BDATVOLDHMR} internally calls \code{\link{BDATVOLABMR}} with parameter
#' \code{A = 0}, i.e. volume is calculated from forest floor up to given
#' diameter, which itself is transformed into height \code{B}.
#'
#' @return vector of same length as input variables transformed into a
#' data.frame, returning the required volume in cubic meter.
#' @examples
#' ## simple call of function, with all parameters given
#' BDATVOLDHMR(1, 30, 1.3, 0, 0, 25, 7, 2.0)
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' DHGrz <- 7
#'
#' ## same using variables
#' BDATVOLDHMR(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H, DHGrz = DHGrz)
#'
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATVOLDHMR(getSpeciesCode(c("Fi", "Bu")), 30, 0, 0, 0, H = 25, DHGrz = 7)
#' ## english abbreviation
#' BDATVOLDHMR(getSpeciesCode(c("NS", "BE")), 30, 0, 0, 0, H = 25, DHGrz = 7)
#' @seealso \code{\link{BDATVOLDHOR}} for BDAT routine calculating volume under
#' bark, \code{\link{getVolume}} for a function with a convenient english name,
#' more options and including a bark switch.
#' @export
BDATVOLDHMR <- function(BDATArtNr, D1, H1 = 1.3, D2 = 0, H2 = 0, H, DHGrz = 7,
                        SekLng = 2.0) {
  HDHGrz <- getHeight(tree = list(
    spp = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H,
    Dx = DHGrz
  ), bark = TRUE)
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H, A = 0, B = HDHGrz, sl = SekLng
  )
  res <- getVolume(tree = df, iAB = "H", bark = TRUE, mapping = NULL)
  return(res)
}
