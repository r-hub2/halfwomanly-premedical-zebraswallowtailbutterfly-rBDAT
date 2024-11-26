#' @title BDAT 2.0 assortment function
#'
#' @description
#' Calculates volumes and assortments for given tree/s.
#' @param BDATArtNr numeric, 1 <= spp <= 36, see \code{\link{getSpeciesCode}}
#' @param D1 numeric, first measured diameter [cm], usually in 1.3m
#' @param H1 numeric, height of first measured diameter [m]
#' @param D2 numeric, second measured diameter [cm], or form parameter,
#'   see \code{\link{buildTree}}.
#' @param H2 H2: numeric, height of second measured diameter [m], or form
#'   parameter, see \code{\link{buildTree}}.
#' @param H numeric, tree height [m]
#' @param lX length of unusable wood at stem foot [m], defaults to 0
#' @param Hkz indicator for tree top, 0 - normal (default), 1 - Wipfelbruch,
#' 2 - Gipfelbruch
#' @param Skz indicator for stem type, defaults to 0, see \code{\link{buildTree}}
#' @param Az minimum cutting diameter over bark [cm], defaults to 0, using
#' tabulated data depending on DBH (not documented)
#' @param Hsh usable stem height, defaults to 0, i.e. 0.7*H
#' @param Zsh minimum cutting diameter under bark for stem wood [cm],
#' defaults to 0, using tabulated data depending on DBH (not documented)
#' @param Zab minimum cutting diameter under bark for top segment [cm],
#' defaults to 0, i.e. 14cm under bark.
#' @param Sokz  type assortment calculation, 0 - no assortment,
#' 1 - mid diameter (Mittenstärke), 2 - Heilbronner Sortierung, defaults to 1
#' @param NMaxFixLng number of fixed length assortments at stem foot, defaults
#' to 0 (no fixed length assortments, irrespective of \code{FixLngDef})
#' @param FixLngDef matrix of 4 * length(spp), having minimum cutting diameter,
#' required assortment length, absolute and relative add-on
#' @param result indicator about what information should be returned
#' @return Using default value of \code{result}, which is 'raw', a list is
#' returned keeping information about the input value and produced assortments -
#' the unprocessed returns from the Fortran code.
#'
#' See \code{\link{getAssortment}} for more details, as this function is
#' internally called.
#' @seealso  \code{\link{getAssortment}} for a more flexible function with a
#' more convenient english name.
#' @examples
#' BDAT20(BDATArtNr = c(1, 1), D1 = c(30, 25), H = c(25, 20)) # returns long list
#' BDAT20(BDATArtNr = c(1, 1), D1 = c(30, 25), H = c(25, 20), result = "Vol")
#' # size class
#' BDAT20(BDATArtNr = c(1, 1), D1 = c(30, 25), H = c(25, 20), result = "Skl")
#' BDAT20(
#'   BDATArtNr = c(1, 1), D1 = c(30, 25), H = c(25, 20), NMaxFixLng = 1,
#'   result = "Fix"
#' )
#' @export
BDAT20 <- function(BDATArtNr, D1, H1 = 0, D2 = 0, H2 = 0, H, lX = 0, Hkz = 0, Skz = 0,
                   Az = 0, Hsh = 0, Zsh = 0, Zab = 0, Sokz = 1, NMaxFixLng = 0, FixLngDef = matrix(rep(
                     0,
                     length(BDATArtNr) * 4
                   ), ncol = 4), result = "raw") {
  df <- data.frame(
    spp = BDATArtNr, D1 = D1, H1 = H1, D2 = D2, H2 = H2, H = H,
    lX = lX, Hkz = Hkz, Skz = Skz, Az = Az, Hsh = Hsh, Zsh = Zsh, Zab = Zab,
    Sokz = Sokz, fixN = NMaxFixLng, fixZ = FixLngDef[, 1], fixL = FixLngDef[
      ,
      2
    ], fixA = FixLngDef[, 3], fixR = FixLngDef[, 4]
  )

  res <- getAssortment(tree = df, value = result)
  return(res)
}
