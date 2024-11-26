#' @title Get assortments for one or many trees
#' @description Function to get assortments given harvest specifications for a
#' tree of dimension D1, possible D2, and H. Assortments could be diameter
#' or length constrained.
#' @param tree either an object of class 'datBDAT.assortment' or a data.frame or
#'  list containing the variables needed, i.e. spp, D1, H, and optionally
#'  H1, D2, H2 for specifying the tree. See \code{\link{buildTree}} for more
#'  information and parameter \code{mapping} for mapping of variable names.
#'  Indeed, \code{tree} could additionally take the variables specified in
#'  \code{sort}; otherwise a full outer join is produced between \code{tree}
#'  and \code{sort}
#' @param sort named list with variables specifying assortments, see
#' \code{\link{buildTree}}.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} between colnames(\code{tree}) and required parameter
#' names. See details.
#' @param value character vector indicating return type: either "Vol", "Skl",
#' "Fix", "LDSort", "merge" (default) or "raw". See section Value.
#' @param ... passing arguments to methods.
#' @details Parameter 'tree' is able to take either a data.frame with correct
#' variables names or arbitrary names if \code{mapping} is provided to map the
#' data.frame names to the required names by \code{c("df-colname" = "var-name")}
#' or to take a named list.
#'
#' Assortments are calculated using BDAT2.0 Fortran routines. Slightly extended,
#' it now is possible to return length and diameter information also about
#' standard assortments (which are held in "Vol"-element of return list,
#' if value="raw" or value="Vol").
#'
#' The standard assortment names are:
#' \itemize{
#'   \item X = non-usable wood at stem foot (X-Holz)
#'   \item Sth = stem wood
#'   \item Ab = upper part of stem wood, second length after transport cut
#'   \item Ind = industrial wood
#'   \item nvDh = non-usable coarse wood
#' }
#'
#' @return depending on \code{value} either information about 'Skl,' 'Vol' or
#' 'Fix' (these are elements of standard output of Fortran BDAT20 subroutine)
#' and, additionally, information about the 'Vol' elements is retrieved by
#' \code{value="LDSort"} or
#' if \code{value} = 'merge' than a combined information of all produced
#' assortments with name, base position, assortment length, mid-diameter,
#' top-diameter and volume is produced. This is most likely what you want,
#' hence the default. Standard output of BDAT is provided by \code{value="raw"}.
#' Since v0.4.0 a vectorized BDAT-fortran-function is implemented, so each
#' element of the returned list in case of value="raw" keeps information of all
#' trees given. Make sure to appropriately split and manipulate this data.
#'
#' @examples
#' tree <- data.frame(spp = 1, D1 = 30, H = 25)
#' sort <- list(Az = 7, Sokz = 1)
#' getAssortment(tree, sort, value = "Vol")
#' getAssortment(tree, sort, value = "Skl")
#' sort <- list(Az = 7, Sokz = 1, fixN = 1, fixZ = 10, fixL = 5,
#'              fixA = 10, fixR = 1)
#' getAssortment(tree, sort, value = "Vol")
#' getAssortment(tree, sort, value = "Skl")
#' getAssortment(tree, sort, value = "Fix")
#' getAssortment(tree, sort, value = "LDSort")
#' getAssortment(tree, sort, value = "merge")
#'
#' ## prepare data for repeated sorting
#' ## (get rid of preparating data handling)
#' n <- 3
#' tree <- data.frame(
#'   spp = rep(1, n), D1 = seq(20, 50, length.out = n),
#'   H = seq(15, 40, length.out = n)
#' )
#' sort <- list(lX = 0, Sokz = 1, Az = 7,
#'              fixN = 2, fixZ = 10, fixL = 5,
#'              fixA = 10, fixR = 1)
#' tree <- buildTree(tree = tree, check = "assortment", vars = sort)
#' getAssortment(tree, value = "Vol")
#' getAssortment(tree, value = "Skl")
#' getAssortment(tree, value = "Fix")
#' getAssortment(tree, value = "LDSort")
#' getAssortment(tree, value = "merge")
#'
#' ## to get bare BDAT-Output, use value='raw'
#' # very long list, each element keeping all trees
#' getAssortment(tree, value="raw")
#' # bonus: it returns the calling parameters as well
#'
#' @seealso \code{\link{BDATVOLABMR}} and \code{\link{BDATVOLABOR}} for a
#' BDAT-type wrapper and function to keep back-compatability. To get coarse wood
#' volume over bark (>=7m diameter over bark) BDAT-type functions are
#' \code{\link{BDATVOLDHMR}} and \code{\link{BDATVOLDHOR}}.
#' @useDynLib rBDAT, .registration=TRUE
#' @export
getAssortment <- function (tree, ...) {
  UseMethod("getAssortment", tree)
}

#' @describeIn getAssortment transforming \code{data.frame} before calling
#' \code{getAssortment} using \code{buildTree}
#' @export
getAssortment.data.frame <- function(tree, sort = NULL, mapping = NULL,
                                     value = "merge", ...){
  tree <- buildTree(tree,
                    check = "assortment",
                    vars = sort,
                    mapping = mapping)
  getAssortment(tree, sort = sort, mapping = mapping, value = value)
}

#' @describeIn getAssortment transforming \code{list} before calling
#' \code{getAssortment} using \code{buildTree}
#' @export
getAssortment.list <- function(tree, sort = NULL, mapping = NULL,
                               value = "merge", ...){
  tree <- buildTree(tree,
                    check = "assortment",
                    vars = sort,
                    mapping = mapping)
  getAssortment(tree, sort = sort, mapping = mapping, value = value)
}

#' @describeIn getAssortment class method for class \code{datBDAT}
#' @export
getAssortment.datBDAT <- function(tree, sort = NULL, mapping = NULL,
                                  value = "merge", ...) {
  value <- ifelse(identical(value, "Fix"), "FixLng", value)
  if (!(value %in% c("Vol", "Skl", "LDSort", "FixLng", "merge", "raw"))) {
    stop("value should be one of 'Vol', 'Skl' 'LDSort', 'Fix', 'merge' or
         'raw', exactly!")
  }

  if (!("datBDAT.assortment" %in% class(tree))) {
    tree <- buildTree(tree, check = "assortment", vars = sort, mapping = mapping)
  }

  # nc=number of assortment specifications per tree
  nas <- ifelse(!is.null(sort), max(sapply(X = sort, FUN = length)), 1)
  nr <- nrow(tree) / nas # number of rows for final matrix


  n <- nrow(tree)
  res <- .Fortran("vbdat20",
    n = as.integer(n),
    BDATArtNr = as.integer(tree$spp),
    D1 = as.single(tree$D1),
    H1 = as.single(tree$H1),
    D2 = as.single(tree$D2),
    H2 = as.single(tree$H2),
    H = as.single(tree$H),
    lX = as.single(tree$lX), # Laenge X-Holz
    Hkz = as.integer(tree$Hkz), # Hoehenkennziffer
    Skz = as.integer(tree$Skz), # Stammkennziffer
    Az = as.single(tree$Az), # Aufarbeitungzopf
    Hsh = as.single(tree$Hsh), # Stammhöhe
    Zsh = as.single(tree$Zsh), # Stammholzzopf
    Zab = as.single(tree$Zab), # Abschnittszopf
    Sokz = as.integer(tree$Sokz), # Sortierkennziffer
    Skl = as.integer(rep(0, 6 * n)), # Stärkeklasse /AUS/
    Vol = as.single(rep(0, 7 * n)), # Volumenangabe /AUS/
    # COMMON /glLDSort/ : Länge und Durchmesser der Sortimente /AUS/
    LDSort = as.single(rep(0, 20 * n)),
    BHD = as.single(rep(0, n)), # geschätzer/gemessener BHD /AUS/
    Ifeh = as.integer(rep(0, n)), # error indicator
    # c(zopf, länge, absZ, relZ%)
    FixLngDef = as.single(t(matrix(c(
      tree$fixZ, tree$fixL,
      tree$fixA, tree$fixR
    ),
    ncol = 4, byrow = F
    ))), # FIX
    NMaxFixLng = as.integer(tree$fixN),
    FixLng = as.single(rep(0, 180 * n)), # Fixlaengen /AUS/
    NFixLng = as.integer(rep(0, n)), # Anzahl Fixlaengen /AUS/
    PACKAGE = "rBDAT"
  )

  ## error code interception
  errInd <- which(res$Ifeh > 0)
  em <- errormessage(res$Ifeh)
  res <- clearError(res)

  if (value %in% c("Skl", "Vol", "LDSort", "FixLng")) {
    res <- matrix(res[[value]], nrow = n, byrow = T)
    res <- transformBDAT20(assort = res, value = value)

  } else if (identical(value, "merge")) {
    r1 <- matrix(res[["Vol"]], nrow = n, byrow = T)
    r1 <- transformBDAT20(assort = r1, value = "Vol")

    r2 <- matrix(res[["LDSort"]], nrow = n, byrow = T)
    r2 <- transformBDAT20(assort = r2, value = "LDSort")

    r2$Vol <- as.vector(t(r1[-c(1, 2, 8)]))

    res <- matrix(res[["FixLng"]], nrow = n, byrow = T)
    res <- transformBDAT20(assort = res, value = "FixLng")

    res$No <- res$No + 5 # increase index by count of conventional assortments
    res <- rbind(r2, res)

    ## order output
    res <- res[order(res$tree), c(
      "tree", "No", "Sort", "height", "length",
      "midD", "topD", "Vol"
    )]

  } else if (identical(value, "raw")) {

    ## do nothing in case value is 'raw'
  }

  if (nas > 1) {
    tmp <- paste0(rep(1:nr, each = nas), ".", rep(1:nas, nr))
    attr(res, "tree-sort-mapping") <- tmp
  }


  if(length(errInd)>0 & !identical(value, "raw")){
    warning("error indication in subroutine BDAT20:\n",
            paste0("tree ", errInd, ": ", em[errInd], "\n"))
  }

  return(res)
}
