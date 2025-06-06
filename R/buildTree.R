#' @title Build and check tree data for subsequent use in
#' BDAT Fortran subroutines
#' @description this functions takes the data provided and builds a data.frame
#' to be used in other BDAT get*-functions. It discriminates between different
#' type of required output via the \code{check}-parameter. Checks are done on
#' the type and range of the variables given to make sure calls to the
#' Fortran routines do not suffer from type-mismatch (with potential freezing
#' of R).
#' @param tree either a data.frame or a list containing the variables needed,
#' i.e. spp, D1, H and optionally H1, D2, H2. See details for more
#' information and parameter \code{mapping} for mapping of variable names.
#' @param check character vector which indicates the type of required output and
#' determines the checks to be done
#' @param vars named list with additional variables for the specific
#' BDAT-functions; see \code{\link{getDiameter}}, \code{\link{getHeight}},
#' \code{\link{getVolume}}, \code{\link{getBiomass}}, \code{\link{getBark}},
#' \code{\link{getForm}} and \code{\link{getAssortment}}.
#' These variables might be included to \code{tree} as well, see details.
#' @param mapping mapping of variable names in case a data.frame is given into
#' parameter \code{tree} and \code{vars} between final \code{colnames(tree)} and
#' required parameter names. See details.
#' @details Parameter \code{tree} is able to take either a data.frame with
#' correct variables names or arbitrary names if \code{mapping} is provided to
#' map the data.frame names to the required names by
#' \code{c("df-colname" = "var-name")} or to take a named list. If same-named
#' variables are present in both \code{tree} and \code{vars}, priority is put
#' on the ones in \code{vars} since explicitly given.
#'
#' Possible variables are (*=required, depending on function):
#' \itemize{
#'   \item spp*: numeric, 1 <= spp <= 36, see \code{\link{getSpeciesCode}}
#'   \item D1*: numeric, first measured diameter [cm], usually at 1.3m
#'   \item H1: numeric, height of first measured diameter [m], if zero,
#'   internally transformed to 1.3m
#'   \item D2: numeric, second measured diameter [cm], or form parameter: latter
#'   is defined in conjunction with \code{H2}:
#'   \itemize{
#'     \item D2=0 and H2=0 => taper form of volume tables according to
#'     Grundner & Schwappach (1906-1938), the default
#'     \item D2=0 and 0 < H2 < 100 => german NFI1-taper form, with H2 given as
#'     percentile of the NFI1-\eqn{q_{0.30}}{q0.30}-distribution; H2=50
#'     corresponds to mean NFI1 taper form, H2<50 to slenderly and H2>50 to
#'     thicker trees; see \code{\link{getForm}} for more information about
#'     \eqn{q_{0.30}}{q0.30}
#'     \item D2=0 and H2>100 => mean NFI1 taper form
#'     \item D2>0 and H2=0 => D2 is a diameter and H2 is assumed to be 7m
#'     \item D2>0 and H2>0 => D2 and H2 are given as diameter and height
#'     \item -1<D2<0 => abs(D2) is interpreted as \eqn{q_{0.30}}{q0.30}
#'     \item -1>D2 => mean NFI1 taper form
#'   }
#'   \item H2: numeric, height of second measured diameter [m], or in
#'   conjunction with \code{D2}, see there.
#'   \item H*: numeric, tree height [m]
#'   \item A*: numeric, lower diameter [cm] or height [m] of section for which
#'   volume should be calculated, interpretation depends on \code{iAB}, see
#'   \code{\link{getVolume}}
#'   \item B*: numeric, upper diameter [cm] or height [m] of section for which
#'   volume should be calculated, interpretation depends on \code{iAB}, see
#'   \code{\link{getVolume}}
#'   \item sl: numeric, length of section over which should be integrated,
#'   defaults to 2.0m
#'   \item Dx*: diameter for which height or bark thickness is required
#'   \item Hx*: height for which diameter is required
#'   \item inv: inventory for which mean q03 is required, defaults to 1, see
#'   \code{\link{getForm}}
#' }
#'
#' For deriving assortments, the following variables are optional (if not given,
#' default values are used):
#' \itemize{
#'   \item lX: length of unusable wood at stem foot [m], defaults to 0 (X-Holz)
#'   \item Hkz: indicator for tree top, 0 - normal (default), 1 - Wipfelbruch, 2 - Gipfelbruch
#'   \itemize{
#'     \item 0 => H=H
#'     \item 1 => H=H+2
#'     \item 2 => DBH < 30 => H=DBH; dbh > 30 => H = 30 + (DBH-30) * 0.3
#'   }
#'   \item Skz: indicator for stem type, defaults to 0
#'   \itemize{
#'     \item 0 => conifer trees => no restriction; deciduous trees => no assortments
#'     \item 1 => monopodial deciduous trees => Hsh = 0.7*H
#'     \item 2 => branching between dbh and 7m => Hsh = 5m
#'     \item 3 => crown base < 3m => Hsh=0.1
#'     \item 4 => dead or broken stem => Az = H*0.7
#'     \item 5 => dead tree => non-usable wood
#'   }
#'   \item Hsh: usable stem height, defaults to 0, i.e. 0.7*H
#'   \item Az: minimum cutting diameter over bark [cm], defaults to 0,
#'   using an exponential function given DBH
#'   \item Zsh: minimum cutting diameter under bark for stem wood [cm], defaults
#'   to 0, using parameter \code{Az} if estimated length < maximum length (i.e. 20m)
#'   \item Zab: minimum cutting diameter under bark for top segment [cm], defaults
#'   to 0, i.e. 14cm under bark
#'   \item Sokz: type assortment calculation, 0 - no assortment,
#'   1 - Mid diameter (Mittenstärke), 2 - Heilbronner Sortierung, defaults to 1
#'   \item fixN: number of fixed length assortments at stem foot, defaults to 0
#'   (no fixed length assortments, irrespective of other fix* parameters)
#'   \item fixZ: mininum diameter under bark for fixed length assortment at
#'   stem foot, defaults to 0
#'   \item fixL: length of fixed length assortment at stem foot, defaults to 0
#'   \item fixA: fixed length assortement add-on in [cm], defaults to 0
#'   \item fixR: fixed length assortement add-on in [\%], defaults to 0
#' }
#'
#' If parameter \code{tree} is used to hand over all tree data in form of
#' a data.frame, at least the parameter spp, D1, H must be provided, eventually
#' mapped via \code{mapping}.
#' Parameter \code{Hx} and \code{Dx}, which specify height and diameter for
#' which a diameter or height is requested, respectively, can either be included
#' to the definition of the tree data or alternatively given separately using
#' the vars parameter. In that case, vars is used in priority to a identically
#' named variable in \code{tree}. Additionally, \code{tree} and \code{vars} are
#' merged via a full outer join.
#' The add-on in fixed length assortments can be given in absolute and relative
#' units at the same time, but the higher value will be used.
#' @return a data.frame of class datBDAT.<check> having all variables needed in
#' specific functions. If \code{check} is NULL, only a basic tree-data.frame of
#' class "datBDAT" is returned.
#' @examples
#' ## example for only tree data
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' res <- buildTree(tree = tree)
#' head(res)
#' class(res)
#'
#' tree <- list(species = c(1, 1), dbh = c(30, 25), h = c(25, 30))
#' mapping <- c("species" = "spp", "dbh" = "D1", "h" = "H")
#' res <- buildTree(tree = tree, mapping = mapping)
#' head(res)
#' class(res)
#'
#' ## example for diameter calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' vars <- list(Hx = c(1.3, 1.3))
#' mapping <- NULL
#' res <- buildTree(tree = tree, check = "diameter", vars = vars)
#' head(res)
#' class(res)
#' tree <- list(Art = c(1, 1), Bhd = c(30, 25), H = c(25, 30))
#' vars <- list(X = c(1.3, 1.3))
#' mapping <- c("Art" = "spp", "Bhd" = "D1", "X" = "Hx")
#' res <- buildTree(tree = tree, check = "diameter", vars = vars, mapping = mapping)
#' head(res)
#' class(res)
#'
#' ## example with many diameters for one tree
#' tree <- list(spp = c(1), D1 = c(30), H = c(25))
#' vars <- list(Hx = seq(0, 25, 0.1))
#' mapping <- NULL
#' res <- buildTree(tree = tree, check = "diameter", vars = vars)
#'
#' tree <- data.frame(s = 1, d = 30, h = 25, hx = 1.3)
#' mapping <- c("s" = "spp", "d" = "D1", "h" = "H", "hx" = "Hx")
#' res <- buildTree(tree, check = "diameter", mapping = mapping)
#' head(res)
#' class(res)
#'
#' ## example for height calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' vars <- list(Dx = c(30, 25))
#' res <- buildTree(tree = tree, check = "height", vars = vars)
#' head(res)
#' class(res)
#'
#' ## example for volume calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' check <- "volume"
#' vars <- list(A = c(30, 25), B = c(7, 7), sl = 0.1)
#' mapping <- NULL
#' res <- buildTree(tree = tree, check = "volume", vars = vars)
#' head(res)
#' class(res)
#'
#' ## example for bark calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' vars <- list(Hx = c(1.3, 1.3))
#' res <- buildTree(tree = tree, check = "bark", vars = vars)
#' head(res)
#' class(res)
#'
#' ## example for assortment calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' vars <- list(fixN = 1, fixZ = 10, fixL = 5, fixA = 10, fixR = 0.1)
#' res <- buildTree(tree = tree, check = "assortment", vars = vars)
#' head(res)
#' class(res)
#'
#' ## for cases where 'vars' could be a vector (i.e. getBark, getDiameter and
#' ## getHeight), the following is also possible
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' vars <- c(1.3, 1.3)
#' res <- buildTree(tree = tree, check = "bark", vars = vars)
#' head(res)
#' class(res)
#'
#' res <- buildTree(tree = tree, check = "height", vars = vars)
#' head(res)
#' class(res)
#'
#' ## but it is not possible in case of getVolume or getAssortment
#' ## instead, use a named list to achieve a cross join / cartesian product
#' vars <- list(A = rep(1, 3), B = 5:7)
#' res <- buildTree(tree = tree, check = "volume", vars = vars)
#' head(res)
#' class(res)
#'
#' ## example for 'biomass' calculation
#' tree <- list(spp = c(1, 1), D1 = c(30, 25), H = c(25, 30))
#' res <- buildTree(tree = tree, check = "biomass")
#' head(res)
#' class(res)
#'
#' ## example with H1 != 1.3m
#' tree <- list(
#'   spp = c(1, 1), D1 = c(30, 25), H1 = c(2, 2), H = c(25, 30)
#' )
#' res <- buildTree(tree = tree, check = "biomass")
#' head(res)
#' class(res)
#' getBiomass(res)
#' @export

buildTree <- function(tree, check = NULL, vars = NULL, mapping = NULL) {
  if (!is.null(vars) & is.null(check)) {
    stop("'vars' is given, but 'check' is NULL! Bad combination!")
  }
  if(any(sapply(vars, is.null))){
    stop("'vars' is named but some elements are NULL")
  }

  if ("datBDAT" %in% class(tree)) {
    if (!all(c("spp", "D1", "H1", "D2", "H2", "H") %in% colnames(tree))) {
      class(tree) <- class(tree)[-which(class(tree) == "datBDAT")]
      warning("'tree' pretended to be of class 'datBDAT' but was not! Removed class and proceed.")
    }
  }

  if (!("datBDAT" %in% class(tree))) {
    if (is.null(tree)) stop("'tree' must not be NULL!")

    ## a list or data.frame (which is a list as well)
    if (identical(class(tree), "list")) {
      if (!is.null(names(tree))) {
        tree <- as.data.frame(tree)
      } else {
        stop("'tree' must be provided as *named* list!")
      }
    }
    ## a data.frame
    if (is.data.frame(tree)) {
      ## add variables from parameter 'vars'
      if (!is.null(vars) & identical(class(vars), "list")) {
        if (!is.null(names(vars))) {
          slct <- which(!(names(tree) %in% names(vars)))
          # tree <- cbind(tree[slct], as.data.frame(vars))
          tree <- merge(tree[slct], as.data.frame(vars))
          vars <- NULL ## destroy, so its clear it has been added into 'tree'
        } else {
          stop("'vars' must be provided as *named* list!")
        }
      }

      ## adjust data.frame to contain all necessary variables
      if (!("H1" %in% colnames(tree)) & !("H1" %in% mapping)) tree$H1 <- 1.3
      if (!("D2" %in% colnames(tree)) & !("D2" %in% mapping)) tree$D2 <- 0
      if (!("H2" %in% colnames(tree)) & !("H2" %in% mapping)) tree$H2 <- 0

      if (is.null(mapping)) {
        if (!all(c("spp", "D1", "H1", "D2", "H2", "H") %in% colnames(tree))) {
          stop("either provide 'mapping' or name parameter 'tree' appropriately!")
        }
      } else {
        if (!(is.character(mapping) && length(mapping) > 0)) {
          stop("'mapping' is not of type character or is not of appropriate length!")
        }
        for (i in seq(along = mapping)) {
          n <- which(colnames(tree) == names(mapping)[i])
          colnames(tree)[n] <- mapping[i]
        }
        if (!all(c("spp", "D1", "H1", "D2", "H2", "H") %in% colnames(tree))) {
          stop("provided name mapping not sufficient!")
        }
      }
    }


    ## check tree data for correct types and range
    if (!is.numeric(tree$spp)) stop("spp must be numeric!")
    if (any(tree$spp < 1) | any(tree$spp > 36)) stop("spp must be >=1 and <= 36!")
    if (!is.numeric(tree$D1)) stop("D1 must be numeric!")
    if (any(tree$D1 < 0)) stop("D1 must be >= 0!")
    if (!is.numeric(tree$H1)) stop("H1 must be numeric!")
    if (any(tree$H1 < 0) | any(tree$H1 > 2.5)) stop("H1 must be >= 0 and <= 2.5!")
    if (!is.numeric(tree$D2)) stop("D2 must be numeric!")
    if (!is.numeric(tree$H2)) stop("H2 must be numeric!")
    if (any(tree$H2 < 0)) stop("H2 must be >= 0!")
    if (!is.numeric(tree$H)) stop("H must be numeric!")
    if (any(tree$H <= 0)) stop("H must be > 0!")
  } else {

  }

  ## check additional data for correct types and range
  if (identical(check, "diameter")) {
    if (!is.list(vars) & is.numeric(vars)) {
      vars <- list(Hx = vars)
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }

    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      # tree <- cbind(tree[slct], as.data.frame(vars))
      tree <- merge(tree[slct], as.data.frame(vars))
    }
    ## check if variable exists in data.frame
    if (!all(c("Hx") %in% colnames(tree))) {
      stop("variable 'Hx' missing! Maybe mapping not sufficient?")
    }
    ## additionally need to check Hx
    if (!is.numeric(tree$Hx)) stop("Hx must be numeric!")
    if (any(tree$Hx < 0)) stop("Hx must be >=0")
  }
  else if (identical(check, "height")) {
    if (!is.list(vars) & is.numeric(vars)) {
      vars <- list(Dx = vars)
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }
    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      # tree <- cbind(tree[slct], as.data.frame(vars))
      tree <- merge(tree[slct], as.data.frame(vars))
    }
    ## check if variable exists in data.frame
    if (!all(c("Dx") %in% colnames(tree))) {
      stop("variable 'Dx' missing! Maybe mapping not sufficient?")
    }
    ## additionally need to check Hx
    if (!is.numeric(tree$Dx)) stop("Dx must be numeric!")
    if (any(tree$Dx < 0)) stop("Dx must be >=0")
    # maxD <- getDiameter(cbind(tree, Hx=0), bark = bark)
    # if(any(tree$Dx > maxD)) stop("Dx must be < D at height=0!")
  }
  else if (identical(check, "volume")) {
    if (!is.list(vars) & is.numeric(vars)) {
      stop("provided 'vars' make no sense, should be data.frame or *named* list!")
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }
    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      # tree <- cbind(tree[slct], as.data.frame(vars))
      tree <- merge(tree[slct], as.data.frame(vars))
    }
    ## adjust data.frame to contain all necessary variables
    ## A and B must be provided! No defaults!
    if (!("sl" %in% colnames(tree))) tree$sl <- 2.0

    ## check if variable exists in data.frame
    if (!all(c("A", "B", "sl") %in% colnames(tree))) {
      stop("at least one variable of 'A', 'B', 'sl' is missing!  Maybe mapping not sufficient?")
    }
    ## additionally need to check A and B
    if (!is.numeric(tree$A)) stop("A must be numeric!")
    if (any(tree$A < 0)) stop("A must be >=0")
    if (!is.numeric(tree$B)) stop("B must be numeric!")
    if (any(tree$B < 0)) stop("B must be >=0")
    if (!is.numeric(tree$sl)) tree$sl <- rep(2.0, nrow(tree))
    if (any(tree$sl <= 0)) stop("sl should be > 0!")
  }
  else if (identical(check, "assortment")) {
    if (!is.list(vars) & is.numeric(vars)) {
      stop("provided 'vars' make no sense, should be data.frame or *named* list!")
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }
    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      # tree <- cbind(tree[slct], as.data.frame(vars))
      tree <- merge(tree[slct], as.data.frame(vars))
    }
    ## adjust data.frame to contain all necessary variables
    nv <- c(
      "lX", "Hkz", "Skz", "Hsh", "Az", "Zsh", "Zab", "Sokz", "fixN",
      "fixZ", "fixL", "fixA", "fixR"
    )
    for (v in nv) {
      # v <- "lX"
      if (!(v %in% colnames(tree))) {
        ## default values
        if (identical(v, "Sokz")) {
          tree$Sokz <- 1
        }
        else {
          tree[, v] <- 0
        }
      }
      else {
        if (identical(v, "Hkz")) {
          if (!all(tree$Hkz %in% c(0, 1, 2))) stop("Hkz must be 0, 1 or 2!")
        }
        else if (identical(v, "Skz")) {
          if (!all(tree$Skz %in% 0:5)) stop("Skz must be in as.integer(0:5)!")
        }
        else if (identical(v, "Sokz")) {
          if (!all(tree$Sokz %in% 0:2)) stop("Sokz must be 0, 1 or 2!")
        }
        # else if (identical(v, "lX")) {
        #   if (!all(tree$lX <= tree$H)) stop("lX [m] must be below tree height H!")
        # }
        else {
          if (any(tree[, v] < 0)) stop(paste0(v, " must be numeric and >=0!"))
        }
      }
    }
  }
  else if (identical(check, "bark")) {
    if (!is.list(vars) & is.numeric(vars)) {
      vars <- list(Hx = vars)
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }
    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      # tree <- cbind(tree[slct], as.data.frame(vars))
      tree <- merge(tree[slct], as.data.frame(vars))
    }
    ## check if variable exists in data.frame
    if (!all(c("Hx") %in% colnames(tree))) {
      stop("variable 'Hx' missing! Maybe mapping not sufficient?")
    }
    ## additionally need to check Hx
    if (!is.numeric(tree$Hx)) stop("Hx must be numeric!")
    if (any(tree$Hx < 0)) stop("Hx must be >=0")
  }
  else if (identical(check, "biomass")) {

    ## no additional 'vars' necessary (no parameter in function)
    if (!("D13" %in% colnames(tree))) {
      ## if D13 is NOT already available
      if (any(tree$H1 != 0 & tree$H1 != 1.3)) {
        ## H1 = 0 is equal to H1=1.3m i.e. diameter in breast height
        tree$D13 <- getDiameter(tree, Hx = 1.3, bark = TRUE, mapping = NULL)
      } else {
        # if H1 =0 | 1.3
        tree$D13 <- tree$D1
      }
    }
    else {
      ## nothing, yet. (D13 is given)
    }
    if (!("D03" %in% colnames(tree))) {
      ## if D03 is NOT already available
      if (any(abs((tree$H2) - (0.3 * tree$H)) > 0.01)) {
        ## if difference in given H2 and 30% of tree height is bigger than 1cm
        tree$Hx <- 0.3*tree$H
        tree$D03 <- getDiameter(tree, bark = TRUE, mapping = NULL)
        tree$Hx <- NULL
        # in case of non-merchantable trees, see Riedel & Kändler 2017, p. 35
        # D03 not necessary and eventually wrong.
        tree$D03 <- ifelse(tree$D13 < 10, 0, tree$D03)
      } else {
        # if H2 = 0.3 * H
        tree$D03 <- tree$D2
      }
    }
    else {
      ## nothing, yet. (D03 is given)
    }
  }
  else if (identical(check, "form")) {

    ## check if vars is NULL
    if(is.null(vars)){
      vars <- 1 # NULL needs to be replaced
    }
    ## check if vars can further be used
    if (!is.list(vars) & is.numeric(vars)) {
      vars <- list(inv = vars)
    } else if (identical(class(vars), "data.frame")) {
      vars <- as.list(vars)
    }

    ## add 'vars' if provided and not used already (then NULL)
    if (!is.null(vars)) {
      slct <- which(!(names(tree) %in% names(vars)))
      tree <- merge(tree[slct], as.data.frame(vars))
    }

    ## check if given values for 'inv' are inside [0, 4]
    tree$inv[tree$inv < 0 | tree$inv > 4] <- 1

    ## check for D13
    if (!("D13" %in% colnames(tree))) {
      ## if D13 is NOT already available
      if (any(tree$H1 != 0 & tree$H1 != 1.3)) {
        ## H1 = 0 is equal to H1=1.3m i.e. diameter in breast height
        tree$D13 <- getDiameter(tree, Hx = 1.3, bark = TRUE, mapping = NULL)
      } else {
        # if H1 =0 | 1.3
        tree$D13 <- tree$D1
      }
    }
    else {
      ## nothing, yet. (D13 is given.)
    }
  }

  if (is.null(check)) {
    class(tree) <- c("datBDAT", class(tree))
  } else {
    class(tree) <- unique(c(
      paste0(c(
        paste0(c("datBDAT", check), collapse = "."),
        "datBDAT"
      )),
      class(tree)
    ))
  }

  return(tree)
}
