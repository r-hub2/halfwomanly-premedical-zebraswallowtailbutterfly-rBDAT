#' @title Calculate volume under bark of a tree between height A and B
#'
#' @description BDAT-Function to get wood volume under bark of one or many trees
#' of a section between height A and height B
#' @param BDATArtNr numeric vector of species code; see \code{\link{getSpeciesCode}}
#' @param D1 first measured diameter of tree [cm], e.g. diameter in breast height
#' @param H1 measurement height of \code{D1} [m]
#' @param D2 second measured diameter of tree, see \code{\link{buildTree}} for
#' details on how to specify different taper forms
#' @param H2 measurement height of D2, see \code{\link{buildTree}} for details
#' on how to specify different taper forms
#' @param H total tree height [m]
#' @param A lower height of section for which volume is required [m]
#' @param B upper height of section for which volume is required [m]
#' @param SekLng length of section over which the integral of taper form should
#' be applied, defaults to 2.0m
#' @details wood volume is calculated using BDAT Fortran routines.
#'
#' @return vector of same length as input variables transformed into a
#' data.frame, returning the required wood volume in cubic meter.
#' @examples
#' ## simple call of function, with all parameters
#' BDATVOLABOR(1, 30, 1.3, 0, 0, 25, .25, 5.25, 2.0)
#' BDATArtNr <- 1
#' D1 <- 30
#' H1 <- 1.3
#' D2 <- 0
#' H2 <- 0
#' H <- 25
#' A <- 1
#' B <- 10
#'
#' ## same with variables
#' BDATVOLABOR(BDATArtNr = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H, A = A, B = B)
#'
#' ## calling with a subset of tree characteristics
#' ## german species names, abbreviated
#' BDATVOLABOR(getSpeciesCode(c("Fi", "Bu")), 30, 0, 0, 0, H = 25, A = 0, B = 25)
#' ## english abbreviation
#' BDATVOLABOR(getSpeciesCode(c("NS", "BE")), 30, 0, 0, 0, H = 25, A = 0, B = 25)
#' @seealso \code{\link{BDATVOLABMR}} for BDAT routine calculating volume over
#' bark, \code{\link{getVolume}} for a function with a convenient english name,
#' more options and including a bark switch.
#' @export
BDATVOLABOR <- function(BDATArtNr, D1, H1 = 1.3, D2 = 0, H2 = 0, H, A, B,
                        SekLng = 2.0) {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1,
    D2 = D2, H2 = H2, H = H, A = A, B = B, sl = SekLng
  )
  res <- getVolume(tree = df, iAB = "H", bark = F, mapping = NULL)
  return(res)
}
