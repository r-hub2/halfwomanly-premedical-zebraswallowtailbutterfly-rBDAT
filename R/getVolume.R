#' @title Get segment volume for one or many trees
#' @description Function to get segment volume over or under bark for a tree of
#' dimension D1, possible D2, and H. Volume is calculated for a given segment
#' defined by A and B, which might be specified as heights or diameters.
#' @param tree either a data.frame or a list containing the variables needed to
#' decribe a tree, i.e. spp, D1, H, and optionally H1, D2, H2. See
#' \code{\link{buildTree}} for details and parameter \code{mapping} for
#'  mapping of variable names
#' @param AB list with heights or diameters A and B of section for which volume
#' over or under bark should be calculated. Additionally, add in \code{sl} for
#' the segment length over which the integral should be calculated. See details.
#' @param iAB character indicating how to interpret given A and B values. Either
#' "H" (the default), "Dob" (diameter over bark) or "Dub" (diameter under bark).
#'  Could be of length one or two, depending on whether A and B are both height
#'  or diameter variables or not. See examples.
#' @param bark boolean of length one indicator for whether required volume should
#' include bark volume or not. Defaults to TRUE. Coerced to logical by
#' \code{\link{as.logical}(bark[1])}.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names. See details.
#' @param ... passing arguments to methods.
#' @details \code{iAB} can be a vector of length two, indicating how to interpret
#' \code{A} and \code{B}. Hence, one can calculate volume between a given
#' height and a given diameter, either over or under bark. If of length one,
#' it is assumed the indicator applies to both A and B.
#'
#' Internally, provided diameters in \code{A} or \code{B} are tranformed to
#' heights by \code{BDATHXDX} if iAB="dob" or by \code{BDATHXDXOR} if
#' iAB="dub".
#' @return vector of same length as tree, returning the required volume in
#' cubic meter.
#' @examples
#' ## default return with just one or several trees given is total coarse wood
#' ## volume over bark, i.e. stock volume (in german: Vfm m.R.)
#' getVolume(list(spp = 1, D1 = 30, H = 25))
#'
#' ## first using a data.frame and height information
#' tree <- data.frame(spp = 1, D1 = 30, H = 25, A = 0.1, B = 10)
#' getVolume(tree) # iAB = "H", bark = TRUE, mapping = NULL as default values
#' getVolume(tree, iAB = "H", bark = FALSE, mapping = NULL)
#'
#' ## now, use diameter information
#' tree <- data.frame(spp = 1, D1 = 30, H = 25, A = 30, B = 7)
#' getVolume(tree, iAB = "Dob", bark = TRUE, mapping = NULL)
#'
#' ## now use both diameter and height information
#' tree <- data.frame(spp = 1, D1 = 30, H = 25, A = 0, B = 7)
#' getVolume(tree, iAB = c("H", "Dob"), bark = TRUE, mapping = NULL)
#' ## is equivalent to the original BDAT-function:
#' BDATVOLDHMR(1, 30, 0, 0, 0, H = 25)
#' ## and
#' getVolume(tree, iAB = c("H", "Dob"), bark = FALSE, mapping = NULL)
#' ## is equivalent to:
#' BDATVOLDHOR(1, 30, 0, 0, 0, H = 25)
#'
#' ## use a misnamed data.frame and mapping argument
#' tree <- data.frame(BDATCode = 1, dbh = 30, h = 25, l = 30, u = 7)
#' getVolume(tree,
#'   iAB = "Dob", bark = TRUE,
#'   mapping = c("BDATCode" = "spp", "dbh" = "D1", "h" = "H", "l" = "A", "u" = "B")
#' )
#'
#' ## using a list to provide the data
#' tree <- list(
#'   spp = c(1, 1), D1 = c(30, 25), H = c(25, 30), A = c(30, 25),
#'   B = c(7, 7)
#' )
#' getVolume(tree, iAB = "Dob") #' defaults: bark = TRUE, mapping = NULL
#' getVolume(tree, iAB = "Dob", bark = FALSE) #' defaults: mapping = NULL
#'
#' ## using a misnamed list and mapping argument
#' tree <- list(
#'   BDATCode = c(1, 1), dbh = c(30, 25), H = c(25, 30),
#'   A = c(0.1, 0.1), B = c(1, 1)
#' )
#' getVolume(tree, mapping = c("BDATCode" = "spp", "dbh" = "D1"))
#'
#' ## using parameter AB to provide the segment data
#' ## in this case a cross join between tree and AB is produced and evaluated
#' tree <- list(BDATCode = c(1, 1), dbh = c(30, 25), H = c(25, 30))
#' AB <- list(A = c(0.1, 1), B = c(1, 2))
#' getVolume(tree, AB, iAB = "H", bark = TRUE, mapping = c("BDATCode" = "spp", "dbh" = "D1"))
#'
#' ## effect of adjusting 'sl'
#' tree <- list(spp = 1, D1 = 30, H = 25)
#' getVolume(tree = tree, AB = list(A = 0.1, B = 5.1, sl = 2)) # default
#' getVolume(tree = tree, AB = list(A = 0.1, B = 5.1, sl = 1))
#' getVolume(tree = tree, AB = list(A = 0.1, B = 5.1, sl = 0.1))
#' getVolume(tree = tree, AB = list(A = 0.1, B = 5.1, sl = 0.01))
#'
#' ## compare Smalian formula to mid-diameter volume for sl=5
#' R1 <- getDiameter(tree, Hx = 1.3) / 100 / 2 # radius in m
#' R2 <- getDiameter(tree, Hx = 6.3) / 100 / 2 # radius in m
#' h <- 5
#' pi * ((R1 + R2) / 2)^2 * h # average of bottom and top diameter
#' (h * pi) / 3 * (R1^2 + R1 * R2 + R2^2) # truncated-cone volume
#' getVolume(tree, AB = list(A = 1.3, B = 6.3, sl = 5)) # mid-diameter volume
#' getVolume(tree, AB = list(A = 1.3, B = 6.3, sl = 0.01)) # physical volume
#'
#' @seealso \code{\link{BDATVOLABMR}} and \code{\link{BDATVOLABOR}} for a
#' BDAT-type wrapper and function to keep back-compatability. To get total
#' coarse wood volume (>=7m diameter over bark) BDAT-type functions are
#' \code{\link{BDATVOLDHMR}} and \code{\link{BDATVOLDHOR}}.
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getVolume <- function(tree, ...){
  UseMethod("getVolume", tree)
}

#' @describeIn getVolume transforming \code{data.frame} before calling
#' \code{getVolume} using \code{buildTree}
#' @export
getVolume.data.frame <- function(tree, AB = NULL, iAB = "H", bark = TRUE,
                                 mapping = NULL, ...){
  ## check if AB is (not) given in any way (i.e. explicitly or implicitly)
  ## if not given, return total coarse wood incl. bark
  if (is.null(AB) & (!all(c("A", "B") %in% names(tree)) &
                     !all(c("A", "B") %in% mapping))) {
    # define parameters for coarse wood volume incl. bark up to 7cm dob
    # german: Vorratsvolumen Vfm m.R.
    AB <- list(A = 0, B = 7)
    iAB <- c("H", "Dob")
  }

  tree <- buildTree(tree,
                    check = "volume",
                    vars = AB,
                    mapping = mapping)
  getVolume(tree, AB = AB, iAB = iAB, bark = bark, mapping = mapping)
}

#' @describeIn getVolume transforming \code{data.frame} before calling
#' \code{getVolume} using \code{buildTree}
#' @export
getVolume.list<- function(tree, AB = NULL, iAB = "H", bark = TRUE,
                          mapping = NULL, ...){
  ## check if AB is (not) given in any way (i.e. explicitly or implicitly)
  ## if not given, return total coarse wood incl. bark
  if (is.null(AB) & (!all(c("A", "B") %in% names(tree)) &
                     !all(c("A", "B") %in% mapping))) {
    # define parameters for coarse wood volume incl. bark up to 7cm dob
    # german: Vorratsvolumen Vfm m.R.
    AB <- list(A = 0, B = 7)
    iAB <- c("H", "Dob")
  }

  tree <- buildTree(tree,
                    check = "volume",
                    vars = AB,
                    mapping = mapping)
  getVolume(tree, AB = AB, iAB = iAB, bark = bark, mapping = mapping)
}

#' @describeIn getVolume class method for class \code{datBDAT}
#' @export
getVolume.datBDAT <- function(tree, AB = NULL, iAB = "H", bark = TRUE,
                      mapping = NULL, ...) {

  ## check if AB is (not) given in any way (i.e. explicitly or implicitly)
  ## if not given, return total coarse wood incl. bark
  if (is.null(AB) & (!all(c("A", "B") %in% names(tree)) &
    !all(c("A", "B") %in% mapping))) {
    # define parameters for coarse wood volume incl. bark up to 7cm dob
    # german: Vorratsvolumen Vfm m.R.
    AB <- list(A = 0, B = 7)
    iAB <- c("H", "Dob")
  }

  if (!("datBDAT.volume" %in% class(tree)) | !is.null(AB)) {
    tree <- buildTree(tree, check = "volume", vars = AB, mapping = mapping)
  }

  bark <- as.logical(bark[1])

  nc <- ifelse(!is.null(AB), max(c(length(AB$A), length(AB$B))), 1) # number of columns for final matrix
  nr <- nrow(tree) / nc # number of rows for final matrix

  if ("datBDAT.volume" %in% class(tree)) {
    n <- nrow(tree) # number of trees to process

    ## transform diameter into height
    if (length(iAB) == 1) {
      iAB <- rep(iAB, 2) ## make indicator available for both A and B
    }

    tmpAB <- c("A", "B")
    for (i in seq(along = tmpAB)) {
      # i <- 2
      if (identical(substr(tolower(iAB[i]), 1, 1), "d")) {
        tree$tmp <- tree[, tmpAB[i]]
        if (identical(tolower(iAB[i]), "dob")) {
          tree$tmp <- as.vector(
            .Fortran("vbdathxdx",
              n = as.integer(n),
              vBDATBArtNr = as.integer(tree$spp),
              vD1 = as.single(tree$D1),
              vH1 = as.single(tree$H1),
              vD2 = as.single(tree$D2),
              vH2 = as.single(tree$H2),
              vH = as.single(tree$H),
              vDx = as.single(tree$tmp),
              vHx = as.single(rep(0, n)),
              vIErr = as.integer(rep(0, n)),
              PACKAGE = "rBDAT"
            )$vHx
          )
        } else if (identical(tolower(iAB[i]), "dub")) {
          tree$tmp <- as.vector(
            .Fortran("vbdathxdxor",
              n = as.integer(n),
              vBDATBArtNr = as.integer(tree$spp),
              vD1 = as.single(tree$D1),
              vH1 = as.single(tree$H1),
              vD2 = as.single(tree$D2),
              vH2 = as.single(tree$H2),
              vH = as.single(tree$H),
              vDx = as.single(tree$tmp),
              vHxor = as.single(rep(0, n)),
              vIErr = as.integer(rep(0, n)),
              PACKAGE = "rBDAT"
            )$vHxor
          )
        } else {
          stop("don't understand iAB for diameter")
        }

        tree[, tmpAB[i]] <- round(tree$tmp, 6)
      }
      tree$tmp <- NULL
    }

    ## correct A and B for D1 == 0 (D1 < 0 caught in buildTree)
    tree$A <- ifelse(tree$D1 <= 0, 0, tree$A)
    tree$B <- ifelse(tree$D1 <= 0, 0, tree$B)

    ## check that A <= B, in case: switch
    tree$tmp <- ifelse(tree$A > tree$B, tree$B, tree$A)
    tree$B <- ifelse(tree$A > tree$B, tree$A, tree$B)
    tree$A <- ifelse(tree$B > tree$tmp, tree$tmp, tree$A)
    tree$tmp <- NULL # remove


    ## choose function depending on bark status
    if (identical(bark, TRUE)) {
      res <- as.vector(
        .Fortran("vbdatvolabmr",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vHges = as.single(tree$H),
          vA = as.single(tree$A),
          vB = as.single(tree$B),
          vSekLng = as.single(tree$sl),
          vIErr = as.integer(rep(0, n)),
          vVolABmr = as.single(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vVolABmr
      )
    }
    else {
      res <- as.vector(
        .Fortran("vbdatvolabor",
          n = as.integer(n),
          vBDATBArtNr = as.integer(tree$spp),
          vD1 = as.single(tree$D1),
          vH1 = as.single(tree$H1),
          vD2 = as.single(tree$D2),
          vH2 = as.single(tree$H2),
          vHges = as.single(tree$H),
          vA = as.single(tree$A),
          vB = as.single(tree$B),
          vSekLng = as.single(tree$sl),
          vIErr = as.integer(rep(0, n)),
          vVolABor = as.single(rep(0, n)),
          PACKAGE = "rBDAT"
        )$vVolABor
      )
    }
  }

  return(matrix(res, nrow = nr, ncol = nc)[, , drop = TRUE])
}
