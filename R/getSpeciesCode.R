#' @title Get BDAT species code or transform it to a name.
#' @description Function to get BDAT species code, or transform it to a
#' german or english name, possibly an abbreviated version or even a scientific
#' name
#' @param inSp species information given, either numeric or character
#' @param outSp character vector of names, for which information should be returned
#' @param version numeric vector giving which look-up table version is required.
#' Defaults to \code{2}, i.e. the latest version (currently from june 2025).
#' Alternatively, one can use \code{version=1} to get the previous version. See
#' details.
#' @details The function matches inSp to outSp. Depending on inSp, being either
#' a numeric vector of values between 1 and 36 or a character vector of species
#' names. Possible names are those which could be return values. One can get all
#' names and the respective species code by calling the function with
#' \code{inSP=NULL} and \code{outSP=NULL} (the default).
#'
#' English species names and codes are now (June 2025) taken from the German
#' National Forest Inventory (\code{version=2}, default), in the previous
#' version (\code{version=1}) the English species names were taken from
#' https://www.forestry.gov.uk/pdf/PF2011_Tree_Species.pdf/$FILE/PF2011_Tree_Species.pdf
#' while slightly adjusting the codes to be unique compared to the german codes
#' (e.g. European larch is now ELA instead of EL).
#'
#' The version change only refers to column \code{long}, not to column \code{short}.
#'
#' Any given species code outside the interval [1, 36] is set back to the code 1
#' (i.e. European spruce (version=2) or Norway spruce (version=1)), while
#' throwing a warning. If any inSp - name is invalid, i.e. not in species list,
#' this throws an error.
#'
#' All elements of outSp, which are not column names of the default returned
#' data.frame, are silently dropped.
#'
#'
#' @return vector or data.frame, depending on length of 'outSp'.
#' @examples
#' getSpeciesCode(inSp = NULL, outSp = NULL) ## the default
#' getSpeciesCode() ## the same
#' getSpeciesCode(outSp = "scientific")
#' getSpeciesCode(inSp = c(1, 2)) ## giving codes
#' getSpeciesCode(inSp = c(1, 2, -1, 37)) ## values outside [1, 36] are given code 1
#' getSpeciesCode(inSp = c(1, 2), outSp = c("scientific")) ## output a vector
#' getSpeciesCode(inSp = c("Bu", "Fi")) ## asking for codes of abbreviated german names
#' getSpeciesCode(inSp = c("Bu", "Fi", "Bu")) ## order is preserved
#' getSpeciesCode(inSp = c("Buche", "Fichte")) ## asking for codes of german names
#' getSpeciesCode(inSp = c("BE", "NS")) ## ... abbreviated english names
#' getSpeciesCode(inSp = c("Common beech", "European spruce")) ## ... english names
#' getSpeciesCode(inSp = c("Fagus sylvatica", "Picea abies")) ### ... scientific names
#'
#' # get previous species look-up table
#' getSpeciesCode(version=1)
#' # as compared to the current version
#' getSpeciesCode() # or explicitly
#' getSpeciesCode(version=2)
#'
#' @export

getSpeciesCode <- function(inSp = NULL, outSp = NULL, version=c(2, 1)) {

  stopifnot("'version' must be in c(1,2)" = is.numeric(version))
  stopifnot("'version' must be in c(1,2)" = version %in% c(1,2))

  ## look-up table
  df <-
    structure(list(
      ID = 1:36,
      kurz = c(
        "Fi", "SF", "Ta", "KT", "Kie", "SK",
        "WK", "DG", "La", "EL", "JL", "Th", "Ts", "SN", "Bu", "HB", "Ei", "RE",
        "Pa", "BP", "Es", "Ah", "BA", "SA", "FA", "Bi", "Li", "Er", "Kir", "Ul",
        "Ro", "El", "Ka", "We", "LB", "VB"
      ),
      lang = c(
        "Fichte", "Sitka-Fichte", "Tanne", "Kuestentanne",
        "Kiefer", "Schwarzkiefer", "Weymouthskiefer", "Douglasie", "Laerche",
        "Europ. Laerche", "Jap. Laerche", "Thuja", "Tsuga", "sNB", "Buche",
        "Hainbuche", "Eiche", "Roteiche", "Pappel", "Balsampappel", "Esche",
        "Ahorn", "Bergahorn", "Spitzahorn", "Feldahorn", "Birke", "Linde",
        "Erle", "Kirsche", "Ulme", "Robinie", "Elsbeere", "Kastanie",
        "Weide", "sLB", "Vogelbeere"
      ),
      short = c(
        "NS", "SS", "ESF", "GF", "SP", "AUP", "WEP",
        "DF", "XLA", "ELA", "JLA", "RC", "WH", "XC",
        "BE", "HBM", "OK", "ROK", "XPO", "BPO", "AH", "XAH", "SY",
        "NOM", "FM", "XBI", "LI", "AR", "WCH", "EM", "BL", "WST",
        "SC", "XWL", "XB", "ROW"
      ),
      long = c(
        "European spruce", "Sitka spruce", "European silver fir",
        "Grand fir", "Common pine", "European black pine", "Weymouth pine",
        "Douglas fir", "larch", "European larch", "Japanese larch",
        "Arbor vitae", "Canada hemlock", "other conifers",
        "Common beech",
        "hornbeam", "Common oak/sessile oak", "Red oak", "poplar",
        "Balsam poplar", "Common ash", "maple", "sycamore maple", "plane maple",
        "field maple", "birch", "lime tree", "alder", "cherry",
        "elm", "Common robinia", "Wild service tree", "European chestnut",
        "willow", "other broadleaves", "rowan berry"
      ),
      long_v1 = c(
        "Norway spruce", "Sitka spruce", "European silver fir",
        "Grand fir", "Scotch pine", "Austrian pine", "Weymouth pine",
        "Douglas fir", "larch", "European larch", "Japanese larch",
        "Western red cedar", "Western hemlock", "other conifers",
        "beech",
        "hornbeam", "oak (robur/petraea)", "Red oak", "poplar",
        "Balsam poplar", "ash", "maple", "sycamore", "Norway maple",
        "Field maple", "birch", "lime tree", "alder", "Wild cherry",
        "elm", "black locust", "Wild service tree", "Sweet chestnut",
        "willow", "other broadleaves", "rowan"
      ),
      scientific = c(
        "Picea abies", " Picea sitchensis", "Abies alba",
        "Abies grandis", "Pinus sylvestris", "Pinus nigra",
        "Pinus strobus", "Pseudotsuga menziesii", "Larix spp.",
        "Larix decidua", "Larix kaempferi", "Thuja plicata",
        "Tsuga heterophylla", "Coniferales trees",
        "Fagus sylvatica",
        "Carpinus betulus", "Quercus spp.", "Quercus rubra",
        "Populus spp.", "Populus balsamifera", "Fraxinus excelsior",
        "Acer spp.", "Acer pseudoplatanus", "Acer platanoides",
        "Acer campestre", "Betula spp.", "Tilia spp.", "Alnus spp.",
        "Prunus avium", "Ulmus spp.", "Robinia pseudoacacia",
        "Sorbus torminalis", "Castanea sativa", "Salix spp.",
        "Magnoliopsida trees", "Sorbus aucuparia"
      )
    ),
    .Names = c("ID", "kurz", "lang", "short", "long", "long_v1", "scientific"),
    class = "data.frame", row.names = c(NA, -36L)
    )

  ## check which version is required
  if(identical(version, 1)){
    df$long <- df$long_v1
  }
  df$long_v1 <- NULL

  if (!is.null(inSp)) {
    # check outSp
    outSp <- outSp[which(outSp %in% colnames(df))]

    # check inSp
    if (is.numeric(inSp)) {
      if (is.null(outSp)) {
        outSp <- c("ID", "kurz", "lang", "short", "long", "scientific")
      }
      if (any(inSp < 0) | any(inSp > 36)) {
        inSp[which(inSp < 0 | inSp > 36)] <- 1
        warning("some 'inSp' where <0 or >36, set to 1")
      }
      outdf <- merge(df, data.frame(x = inSp, orderx = 1:length(inSp)),
        by.x = "ID", by.y = "x", all.y = T
      )
      outdf <- outdf[order(outdf$orderx), outSp]
    } else if (is.character(inSp)) {
      if (is.null(outSp)) {
        outSp <- c("ID")
      }
      if (all(inSp %in% df$kurz)) {
        scol <- "kurz" ## *s*elected *col*umn
      } else if (all(inSp %in% df$lang)) {
        scol <- "lang"
      } else if (all(inSp %in% df$short)) {
        scol <- "short"
      } else if (all(inSp %in% df$long)) {
        scol <- "long"
      } else if (all(inSp %in% df$scientific)) {
        scol <- "scientific"
      } else {
        print(inSp[which(!(inSp %in% unlist(df[, 2:6])))])
        stop("cannot match given 'inSp' to any column. Check your names.")
      }
      ## match
      outdf <- merge(data.frame(x = inSp, orderx = 1:length(inSp)), df,
        by.x = "x", by.y = scol, all.x = T, sort = F
      )
      outdf <- outdf[order(outdf$orderx), outSp]
    } else {
      warning("inSp of wrong mode, must be numeric or character.")
    }
  } else {
    if (!is.null(outSp)) {
      outdf <- subset(df, select = outSp[which(outSp %in% colnames(df))])
    } else {
      outdf <- df
    }
  }

  return(outdf)
}
